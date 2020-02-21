"""
Functions to build and train neural networks and autoencoder
"""

import argparse
import glob
from scipy import interp
import shutil  # remove nonempty dirs
from datetime import datetime
from sklearn.preprocessing import StandardScaler, MaxAbsScaler, LabelEncoder
from sklearn.metrics import roc_curve, auc, precision_recall_curve, f1_score
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.model_selection import KFold, StratifiedKFold
from keras.layers import Input, Dense, Dropout, Flatten, Activation
from keras.layers.normalization import BatchNormalization
from keras.models import Model, load_model
from keras.optimizers import Adam, SGD, Nadam, Adadelta, Adagrad
from keras.constraints import nonneg
from keras import backend as K
from keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping, CSVLogger

from setup import get_param_str
from custom_callbacks import *
from sparse_x_generator import *
from bootstrap import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def none_or_str(value):
    if value == 'None':
        return None
    return value


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataType', help='output data type', choices=['vasari', 'verhaak', 'radiomic'], default='vasari')
    parser.add_argument('--predType', help='prediction type', choices=['binaryClass', 'multiClass', 'regression'], default='binaryClass')
    parser.add_argument('--exp', help='experiment name')
    parser.add_argument('--dir', help='full path to experiment directory, where all data will be saved')
    parser.add_argument('--data', help='full path to data directory')
    parser.add_argument('--label', help='name of supervised label, must be at least a partial match, ignored if autoencoder')

    parser.add_argument('--folds', help='number of cross validation folds', type=int, default=10)
    parser.add_argument('--retrain', help='retrain model, 1 = True/yes 0 = False/no', type=int, default=0)
    parser.add_argument('--pretrain', help='path to pretrained model', type=none_or_str, default=None)
    parser.add_argument('--freeze',
                        help='in pretraining, freeze nl layers?, 1 = True/yes freeze 0 = False/no/just set the weights',
                        type=int, default=0)
    parser.add_argument('--num_ae_layers', help='if using ae, the num of layers to set or freeze nn weights', type=int, default=None)
    parser.add_argument('--gpu', help='gpu num', type=int, default=0)
    parser.add_argument('--seed', help='cv split seed', type=int, default=425)
    parser.add_argument('--bs_iter', help='boostrap ith iteration', type=int, default=None)
    parser.add_argument('--bs_method', help='boostrap method', type=int, choices=[1, 2], default=1)

    parser.add_argument('--opt', help='optimizer name')
    parser.add_argument('--loss', help='loss function', default='binary_crossentropy')
    parser.add_argument('--act', help='activation function')
    parser.add_argument('--drop', help='dropout', type=float, default=0)
    parser.add_argument('--decay', help='decay', type=float, default=0)
    parser.add_argument('--batch', help='batch size', type=int, default=10)
    parser.add_argument('--epoch', help='number of epochs', type=int, default=500)
    parser.add_argument('--patience', help='number of epochs to wait for improvement', type=int, default=200)
    parser.add_argument('--learn', help='default or specify value', type=float, default=0)
    parser.add_argument('--h1', help='nodes in hidden layer 1', type=int, default=0)
    parser.add_argument('--h2', help='nodes in hidden layer 2', type=int, default=0)
    parser.add_argument('--h3', help='nodes in hidden layer 3', type=int, default=0)
    parser.add_argument('--h4', help='nodes in hidden layer 4', type=int, default=0)
    parser.add_argument('--h5', help='nodes in hidden layer 5', type=int, default=0)

    args = parser.parse_args()
    args = vars(args)  # convert to dictionary
    return args


def setup_experiment(args, overwrite=True):
    """ Create directories to save experiment data.

    Args:
        args: dictionary containing experiment parameters, see get_args()
        overwrite: boolean to overwrite files if already exist in experiment directory

    Returns:
        strings:
        experiment directory
        tensorboard directory
        experiment log filename
        input scaler filename
        output scaler filename
        parameter string to name neural networks in experiment - based on hyperparameters
    """

    # experiment folder structure
    exp_dir = os.path.join(args['dir'], args['exp'], args['label'], 'neuralnets')
    tb_dir = os.path.join(args['dir'], args['exp'], args['label'], 'tensorboard')

    # parameter specific, will remove data in order to overwrite previous experiment parameters
    param_str = get_param_str(args)
    exp_dir = os.path.join(exp_dir, param_str)
    tb_dir = os.path.join(tb_dir, param_str)

    if args['bs_iter'] is None:
        d = [exp_dir, tb_dir]
    else:
        d = [exp_dir]  # don't create tensorboard dirs for bootstrap

    if overwrite:
        for dirName in d:
            if os.path.exists(dirName):
                shutil.rmtree(dirName)
            os.makedirs(dirName)

    log_file = exp_dir + '/log_' + datetime.now().strftime("%Y_%m_%d") + '.txt'
    scale_in = exp_dir + '/geneTrainScalers'
    scale_out = exp_dir + '/imgTrainScalers'
    return exp_dir, tb_dir, log_file, scale_in, scale_out, param_str


def get_opt(opt_name):
    return {
        'Adam': Adam(),
        'SGD': SGD(),
        'Nadam': Nadam(),
        'Adagrad': Adagrad(),
        'Adadelta': Adadelta()
    }[opt_name]


def r2(y_true, y_pred):
    """Calculate the coefficient of determination
     taken from https://jmlb.github.io/ml/2017/03/20/CoeffDetermination_CustomMetric4Keras/
    """

    SS_res = K.sum(K.square(y_true-y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return 1 - SS_res/(SS_tot + K.epsilon())


def my_mape(y_true, y_pred):
    """Calculate mean absolute percent error"""

    diff = np.abs((y_true-y_pred) / abs(y_true))
    return 100. * np.mean(diff, axis=0)


def get_nonneg(y):
    """Check if any outputs are negative

    Args:
        y: outputs in dataset

    Returns: if y is nonneg

    """
    nonneg = True
    if np.min(y) < 0:
        nonneg = False
    return nonneg


def get_callbacks(args, train=[None, None], val=[None, None]):
    """

    Args:
        args: dictionary containing experiment parameters, see get_args()
        train: training data: [train_x, train_y]
        val: validation data: [val_x, val_y]

    Returns: list of callbacks

    """
    mf = os.path.join(args['exp_dir'], args['model_name'])
    lf = os.path.join(args['exp_dir'], args['logger'])

    vm, vmode, min_d, p = get_monitor_info(args['predType'], val=val[0], patience=args['patience'])

    early = EarlyStopping(monitor=vm, min_delta=min_d, patience=p, verbose=2, mode=vmode)

    if args['save']:
        checks = ModelCheckpoint(mf, monitor=vm, verbose=2, save_best_only=True, mode=vmode, period=1)
        callbacks = [early, checks, CSVLogger(lf)]
    else:
        callbacks = [early, CSVLogger(lf)]  # cross-validation, don't save model

    print('min d: ', min_d)

    if args['bs_iter'] is None:
        tb = args['tb_dir']
        if val[0] is not None:
            tb = args['cv_tb_dir']
        callbacks.insert(0, TensorBoard(log_dir=tb))

    if args['predType'] == 'binaryClass':
        callbacks.insert(0, roc_callback(training_data=train, validation_data=val, binary=True))

    if args['predType'] == 'multiClass':
        callbacks.insert(0, roc_callback(training_data=train, validation_data=val, binary=False))
    return callbacks


def encode_y(y):
    # Convert string labels to ints
    y_labels = y.reshape(y.shape[0], )
    le = LabelEncoder()
    y_labels = le.fit_transform(y_labels.ravel())  # strings to ints
    classes = le.classes_
    return y_labels, le, classes


def get_model(args, input_dim, output_dim):
    """Build neural network layers
    If using pretrained autoencoder, assumes the first three hidden layers match

    Args:
        args: dictionary containing experiment parameters, see get_args()
        input_dim: int size of input vector
        output_dim: int size of output vector

    Returns: Keras Model

    """
    kc = None
    bc = None
    opt = get_opt(args['opt'])
    d0 = Dropout(args['drop'], seed=args['h1'])
    d1 = Dropout(args['drop'], seed=args['h2'])
    d2 = Dropout(args['drop'], seed=args['h3'])
    d3 = Dropout(args['drop'], seed=args['h4'])
    d4 = Dropout(args['drop'], seed=args['h5'])
    d5 = Dropout(args['drop'], seed=args['h5'])

    # build network
    input_genes = Input(shape=(input_dim, 1))  # add 3rd dimension for keras-vis
    e = Flatten()(input_genes)
    e = d0(e)
    e = Dense(args['h1'])(e)
    e = BatchNormalization()(e)
    e = Activation(args['act'])(e)
    e = d1(e)
    if args['h2'] > 0:
        e = Dense(args['h2'])(e)
        e = BatchNormalization()(e)
        e = Activation(args['act'])(e)
        e = d2(e)
    if args['h3'] > 0:
        e = Dense(args['h3'])(e)
        e = BatchNormalization()(e)
        e = Activation(args['act'])(e)
        e = d3(e)
    if args['h4'] > 0:
        e = Dense(args['h4'])(e)
        e = BatchNormalization()(e)
        e = Activation(args['act'])(e)
        e = d4(e)
    if args['h5'] > 0:
        e = Dense(args['h5'])(e)
        e = BatchNormalization()(e)
        e = Activation(args['act'])(e)
        e = d5(e)

    if args['predType'] == 'regression':
        metrics = ['mae', 'mape', r2]
        act = 'linear'
        if args['nonneg']:
            kc = nonneg()
            bc = nonneg()
    else:
        metrics = ['acc']
        if args['predType'] == 'multiClass':
            act = 'softmax'
        else:
            act = 'sigmoid'  # binary class

    predictor = Dense(output_dim, name='preds', kernel_constraint=kc, bias_constraint=bc)(e)
    predictor = Activation(act)(predictor)

    model = Model(input_genes, predictor)
    model.compile(loss=args['loss'], optimizer=opt, metrics=metrics)

    model.summary()

    if args['pretrain'] is not None:
        print(' set up ae')
        # set pre-trained weights from ae
        ae = load_model(args['ae_model'], custom_objects={'r2': r2})

        # just set the weights
        model.layers[3].set_weights(ae.layers[1].get_weights())  # 3 dropout layers
        model.layers[7].set_weights(ae.layers[4].get_weights())
        model.layers[11].set_weights(ae.layers[7].get_weights())

        if args['freeze'] == 1:
            # freeze layers
            print(' freeze weights')
            if args['num_ae_layers'] > 0:
                model.layers[3].trainable = False  # 3 dropout layers
            if args['num_ae_layers'] > 1:
                model.layers[7].trainable = False

        del ae

        model = Model(input_genes, predictor)
        model.compile(loss=args['loss'], optimizer=opt, metrics=metrics)
    return model


def get_regression_metrics(y_true, y_pred):
    # metrics for each feature, even in grouped
    mae = mean_absolute_error(y_true=y_true, y_pred=y_pred, multioutput='raw_values')
    mape = my_mape(y_true=y_true, y_pred=y_pred)
    r2_value = r2_score(y_true=y_true, y_pred=y_pred, multioutput='raw_values')
    return mae, mape, r2_value


def get_regression_names():
    res_colnames = ['t_mae', 't_mape', 't_r2', 'v_mae', 'v_mape', 'v_r2']
    res_rename = ['t_mae', 'v_mae', 't_mape', 'v_mape', 't_r2', 'v_r2']
    return res_colnames, res_rename


def fit_model(args, x_train, x_val, y_train, y_val, sample_weights=None):
    """Train model

    Args:
        args: dictionary containing experiment parameters, see get_args()
        x_train: ndarray with shape (sample_size, num_features), training data, inputs
        x_val:   ndarray with shape (sample_size, num_features), validation data, inputs
        y_train: ndarray with shape (sample_size, ), training data, outputs
        y_val:   ndarray with shape (sample_size, ), validation data, outputs
        sample_weights: ndarray with shape (sample_size, )

    Returns: nothing if doing cross validation due to CSVLogger, otherwise, returns model, predictions, and metrics

    """
    model_file = os.path.join(args['exp_dir'], args['model_name'])
    train_size, input_dim = x_train.shape
    if y_train.ndim > 1:
        _, output_dim = y_train.shape
    else:
        output_dim = 1
    print('train size: ', train_size)

    # prep input
    if args['pretrain'] is None:
        print('... preprocessing')
        scaler = StandardScaler()  # std input
        x_train = scaler.fit_transform(x_train)
        if args['save']:
            pickle.dump(scaler, open(args['scale_in']+'.pkl', 'wb'))
    else:
        print('... preprocessing with ae scaler')
        scaler = pickle.load(open(args['ae_scaler'], 'rb'))
        x_train = scaler.transform(x_train)

    x_train = x_train.reshape(len(x_train), input_dim, 1)  # reshape input into 3-D shape

    # prep output
    if args['dataType'] == 'radiomic':   # also scale output
        out_scaler = MaxAbsScaler()
        y_train = out_scaler.fit_transform(y_train)
        if args['save']:
            pickle.dump(out_scaler, open(args['scale_out']+'.pkl', 'wb'))
            pickle.dump(out_scaler, open(args['scale_out']+'.csv', 'wb'))

    # prep val if cross-validation
    train_metrics = None
    val_metrics = None

    if x_val is not None:
        val_size, _ = x_val.shape
        print('val size: ', val_size)
        x_val = scaler.transform(x_val)
        x_val = x_val.reshape(len(x_val), input_dim, 1)

        if args['dataType'] == 'radiomic':  # also scale output
            y_val = out_scaler.transform(y_val)
        val_data = (x_val, y_val)
    else:
        val_data = None
        val_size = None

    callbacks = get_callbacks(args, train=[x_train, y_train], val=[x_val, y_val])

    # train
    model = get_model(args=args, input_dim=input_dim, output_dim=output_dim)
    model.fit(x_train, y_train,
              sample_weight=sample_weights,
              epochs=args['epoch'],
              batch_size=args['batch'],
              shuffle=True,
              validation_data=val_data,
              verbose=2,
              callbacks=callbacks)

    if args['save']:
        # validation metrics
        model = load_model(model_file, custom_objects={'r2': r2})
        train_preds = model.predict(x_train)
        if x_val is None:
            val_preds = None
        else:
            val_preds = model.predict(x_val)

        # metrics for each feature in y
        if args['predType'] == 'regression':
            tmae, tmape, tr2 = get_regression_metrics(y_true=y_train, y_pred=train_preds)
            train_metrics = np.vstack((tmae, tmape, tr2)).T  # rows are features, cols are metrics
            if x_val is not None:
                vmae, vmape, vr2 = get_regression_metrics(y_true=y_val, y_pred=val_preds)
                val_metrics = np.vstack((vmae, vmape, vr2)).T

        return model, train_metrics, val_metrics, train_preds, val_preds


def get_autoencoder(args, input_dim):
    """Build autoencoder

    Args:
        args: dictionary containing experiment parameters, see get_args()
        input_dim: int size of input vector

    Returns: Keras Model

    """
    opt = get_opt(args['opt'])

    # build network
    input_genes = Input(shape=(input_dim,))

    e = Dense(args['h1'])(input_genes)                    # encoding hidden layer 1 (takes genes as inputs)
    e = BatchNormalization()(e)
    e = Activation(args['act'])(e)

    e = Dense(args['h2'])(e)                            # encoding hidden layer 2
    e = BatchNormalization()(e)
    e = Activation(args['act'])(e)

    e = Dense(args['h3'])(e)                         # encoding hidden layer 3 (deep representation)
    e = BatchNormalization()(e)
    e = Activation(args['act'])(e)

    if args['h4'] > 0:
        e = Dense(args['h4'])(e)  # encoding hidden layer 4
        e = BatchNormalization()(e)
        e = Activation(args['act'])(e)

    if args['h5'] > 0:
        e = Dense(args['h5'])(e)  # encoding hidden layer 5
        e = BatchNormalization()(e)
        e = Activation(args['act'])(e)

    if args['h5'] > 0:
        e = Dense(args['h4'])(e)  # decoding hidden layer 1
        e = BatchNormalization()(e)
        e = Activation(args['act'])(e)

    if args['h4'] > 0:
        e = Dense(args['h3'])(e)  # decoding hidden layer 1
        e = BatchNormalization()(e)
        e = Activation(args['act'])(e)

    d = Dense(args['h2'])(e)                        # decoding hidden layer 1
    d = BatchNormalization()(d)
    d = Activation(args['act'])(d)

    d = Dense(args['h1'])(d)                        # decoding hidden layer 2
    d = BatchNormalization()(d)
    d = Activation(args['act'])(d)

    d = Dense(input_dim)(d)   # reconstruction
    d = BatchNormalization()(d)
    decoded = Activation(args['act'])(d)

    # set up model
    autoencoder = Model(input_genes, decoded)
    autoencoder.compile(loss=args['loss'], optimizer=opt, metrics=['mae', 'mse', r2])
    autoencoder.summary()
    return autoencoder


def fit_autoencoder(args, x_train, x_val):
    """Train model

    Args:
        args: dictionary containing experiment parameters, see get_args()
        x_train: ndarray with shape (sample_size, num_features), training data, inputs
        x_val:   ndarray with shape (sample_size, num_features), validation data, inputs

    Returns: nothing if doing cross validation due to CSVLogger, otherwise, returns model, predictions, and metrics

    """

    # setup
    model_file = os.path.join(args['exp_dir'], args['model_name'])
    train_size, input_dim = x_train.shape

    # std input
    scaler = StandardScaler()
    x_train = scaler.fit_transform(x_train)
    if not os.path.isfile(args['scale_in']) & args['save'] is True:
        pickle.dump(scaler, open(args['scale_in'] + '.pkl', 'wb'))

    # train
    ae = get_autoencoder(args=args, input_dim=input_dim)

    # prep val if cross-validation
    val_data = None
    val_metrics = None
    if x_val is not None:
        x_val = scaler.transform(x_val)
        val_data = (x_val, x_val)
        val_size, _ = x_val.shape

    callbacks = get_callbacks(args, train=[x_train, x_train], val=[x_val, x_val])

    ae.fit(x_train, x_train,
           epochs=args['epoch'],
           batch_size=args['batch'],
           shuffle=True,
           verbose=2,
           validation_data=val_data,
           callbacks=callbacks)

    del ae
    K.clear_session()

    if args['save']:
        # validation
        ae = load_model(model_file, custom_objects={'r2': r2})

        # metics
        train_preds = ae.predict(x_train)
        if x_val is None:
            val_preds = None
        else:
            val_preds = ae.predict(x_val)

        # metrics for each feature in y
        tmae, tmape, tr2 = get_regression_metrics(y_true=x_train, y_pred=train_preds)
        train_metrics = np.vstack((tmae, tmape, tr2)).T  # rows are features, cols are metrics
        if x_val is not None:
            vmae, vmape, vr2 = get_regression_metrics(y_true=x_val, y_pred=val_preds)
            val_metrics = np.vstack((vmae, vmape, vr2)).T

        return ae, train_metrics, val_metrics, train_preds, val_preds


def get_split(x, y, y_labels, folds, pred_type, seed=425):
    """Cross-validation fold splits

    Args:
        x: float ndarray (sample_size, feature_size)
        y: int ndarray y class labels
        y_labels: str ndarray (sample_size, ) of original labels
        folds: int folds to create
        pred_type: str type of prediction, for class stratification
        seed: int seed to do split

    Returns: generator, fold indices for splits

    """

    if pred_type == 'regression':
        kf = KFold(n_splits=folds, random_state=seed)
        fold_indices = kf.split(y)
    else:
        skf = StratifiedKFold(n_splits=folds, random_state=seed, shuffle=True)
        if y.ndim > 1:
            fold_indices = skf.split(X=x, y=y_labels)
        else:
            fold_indices = skf.split(X=x, y=y)

    return fold_indices


def cross_validate(args, x, y, cv, model_type, y_labels=None, y_names=None, ids=None, classes=None):
    """

    Args:
        args: dictionary containing experiment parameters, see get_args()
        x: float ndarray (sample_size, feature_size)
        y: int ndarray y class labels
        cv: generator fold indices
        model_type: str, 'nn' for regular feedforward neural net, 'ae' for autoencoder
        y_labels: str ndarray y class names
        y_names: str column names for y labels
        ids: str ndarray id names for rows/samples
        classes: str ndarray of unique classes

    Returns: cross validation results

    """
    args['model_name'] = 'model_cv.h5'
    train_metrics = []
    val_metrics = []

    fold = 0
    for train_index, val_index in cv:

        fold += 1
        args['logger'] = 'train_logger_fold_' + str(fold) + '.csv'
        args['cv_tb_dir'] = os.path.join(args['tb_dir'], 'fold_' + str(fold))
        if os.path.exists(args['cv_tb_dir']):    # overwrite tensorboard events
            shutil.rmtree(args['cv_tb_dir'])

        if args['predType'] == 'regression':
            sample_weights = None
        else:
            sample_weights = compute_sample_weight(class_weight='balanced', y=y_labels[train_index])

        if model_type == 'ae':
            if args['save']:
                nn, t_m, v_m, _, _ = fit_autoencoder(args,
                                                     x_train=x[train_index, ],
                                                     x_val=x[val_index, ])
                train_metrics.append(t_m)
                val_metrics.append(v_m)
                cnames = ['gene']
            else:
                fit_autoencoder(args,
                                x_train=x[train_index,],
                                x_val=x[val_index,])
        else:
            # regular nn
            if args['save']:
                nn, t_m, v_m, t_preds, v_preds = fit_model(args,
                                                           x_train=x[train_index, ],
                                                           y_train=y[train_index, ],
                                                           x_val=x[val_index, ],
                                                           y_val=y[val_index, ],
                                                           sample_weights=sample_weights)
                if args['predType'] == 'regression':
                    train_metrics.append(t_m)
                    val_metrics.append(v_m)
                    cnames = y_names
                else:
                    # get roc plots
                    t_plot = os.path.join(args['exp_dir'], 'train_fold_' + str(fold))
                    v_plot = os.path.join(args['exp_dir'], 'val_fold_' + str(fold))

                    if args['predType'] == 'multiClass':
                        t_roc, _, = multiclass_metrics(classes,
                                                       y_truth=y[train_index, ],
                                                       y_preds=t_preds,
                                                       fn=t_plot,
                                                       t=args['label'],
                                                       plot=args['save'])
                        v_roc, _, = multiclass_metrics(classes,
                                                       y_truth=y[val_index, ],
                                                       y_preds=v_preds,
                                                       fn=v_plot,
                                                       t=args['label'],
                                                       plot=args['save'])
                        cnames = classes
                    else:
                        t_roc, _ = binaryclass_metrics(y_truth=y[train_index, ],
                                                       y_preds=t_preds,
                                                       fn=t_plot,
                                                       t=args['label'],
                                                       plot=args['save'])
                        v_roc, _ = binaryclass_metrics(y_truth=y[val_index, ],
                                                       y_preds=v_preds,
                                                       fn=v_plot,
                                                       t=args['label'],
                                                       plot=args['save'])
                        cnames = [args['label']]

                    # save information
                    train_metrics.append(t_roc)
                    val_metrics.append(v_roc)
                    y_names = None

                if model_type == 'nn':
                    t_preds = pd.DataFrame(data=t_preds, columns=cnames, index=ids[train_index])
                    v_preds = pd.DataFrame(data=v_preds, columns=cnames, index=ids[val_index])
                    t_preds.to_csv(os.path.join(args['exp_dir'], 'train_preds_fold_' + str(fold) + '.csv'))
                    v_preds.to_csv(os.path.join(args['exp_dir'], 'val_preds_fold_' + str(fold) + '.csv'))
            else:
                fit_model(args,
                          x_train=x[train_index,],
                          y_train=y[train_index,],
                          x_val=x[val_index,],
                          y_val=y[val_index,],
                          sample_weights=sample_weights)

            K.clear_session()
            # del nn

    # average results
    if args['save']:
        cv_log, cv_metrics = get_cv_averages(args, train_metrics, val_metrics, y_names)
    else:
        cv_log, cv_metrics = get_cv_averages(args)
    return cv_log, cv_metrics


def get_cv_averages(args, train_metrics=None, val_metrics=None, index_names=None):
    """Calculate average scores based on CSVlogger

    Args:
        args: dictionary containing experiment parameters, see get_args()
        train_metrics: float list of train metric for each fold
        val_metrics: float list of val metric for each fold
        index_names: index names for multiouput regression metric dataframe

    Returns: dataframe of cross validation averages

    """
    # calculate average scores based on logger
    # also return each fold metrics
    cv_metrics = None
    pattern = os.path.join(args['exp_dir'], '*train_logger*.csv')
    m, vmode, _, _ = get_monitor_info(prediction_type=args['predType'], val=1, patience=args['patience'])  # which metric was being used to save models

    cv_log = cv_average(log_file_pattern=pattern,
                        metric=m,
                        vmode=vmode)

    if args['save']:
        # calculate additional scores/metrics, each metric input has different shapes and type
        if args['predType'] == 'multiClass':
            # for classification save end aucs for each class
            res_colnames = ['train_{}'.format(i) for i in list(train_metrics[0].keys())]
            b = ['val_{}'.format(i) for i in list(val_metrics[0].keys())]
            for i in b:
                res_colnames.append(i)

            train_metrics = [list(i.values()) for i in train_metrics]
            val_metrics = [list(i.values()) for i in val_metrics]
            inames = ['fold_' + str(i) for i in range(1, len(train_metrics) + 1)]

            train_metrics = np.vstack((train_metrics, np.mean(np.array(train_metrics), axis=0)))   # average over folds
            val_metrics = np.vstack((val_metrics, np.mean(np.array(val_metrics), axis=0)))
            inames.append('avg')
            cv = pd.DataFrame(data=np.hstack((train_metrics, val_metrics)))
            cv.columns = res_colnames
            cv.index = inames

        elif args['predType'] == 'binaryClass':
            # currently not really necessary, same as logger info
            inames = ['fold_' + str(i) for i in range(1, len(train_metrics) + 1)]
            train_metrics = np.hstack((train_metrics, np.mean(np.array(train_metrics), axis=0)))   # average over folds
            val_metrics = np.hstack((val_metrics, np.mean(np.array(val_metrics), axis=0)))
            inames.append('avg')
            cv = pd.DataFrame(data={'roc': train_metrics,
                                    'val_roc': val_metrics},
                              index=inames)
        else:
            # for regression save individual feature metrics
            # return only fold averages for simplicity
            train_metrics = np.mean(np.array(train_metrics), axis=0)  # average over folds
            val_metrics = np.mean(np.array(val_metrics), axis=0)
            metric_names, metric_renames = get_regression_names()
            cv = pd.DataFrame(data=np.hstack((train_metrics, val_metrics)),
                              columns=metric_names,
                              index=index_names)
            cv_metrics = cv[metric_renames]
    return cv_log, cv_metrics


def cv_average(log_file_pattern, metric, vmode):
    """Calculate cv average based on CSVlogger files for each cross-validation fold

    Args:
        log_file_pattern: path to logger files and pattern of file to search for
        metric: str colname of metric to average by
        vmode: str, get 'min' or 'max' of metric to report

    Returns:

    """
    cv_log = []
    names = []
    for log in glob.glob(log_file_pattern):
        res = pd.read_csv(log)
        if vmode == 'max':
            res = res.iloc[[res[metric].idxmax()]]  # callback was based on val_loss
        else: #min
            res = res.iloc[[res[metric].idxmin()]]  # callback was based on val_loss

        cv_log.append(res)
        names.append(os.path.basename(log))

        # if args['save'] is False:
        #     os.remove(log)  # e.g., don't need to save during bootstrap

    cv_log = pd.concat(cv_log)
    avg_cv = cv_log.mean().to_frame().T
    cv_log = pd.concat([cv_log, avg_cv])
    names.append('avg')
    cv_log.index = names

    return cv_log


def binaryclass_metrics(y_truth, y_preds, fn=None, t=None, plot=True):
    """Get area under the curve for roc and precission-recall, can also plot roc

    Args:
        y_truth: int ndarray true class
        y_preds: float ndarray predictions, class probability
        fn: file name, full path
        t: plot title
        plot: whether to plot roc

    Returns: auc of roc and precission-recall

    """
    fpr, tpr, _ = roc_curve(y_truth, y_preds)
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y_truth, y_preds)
    pr_auc = auc(recall, precision)

    # plot
    if plot:
        color = 'darkorange'
        lw=2
        plt.figure()
        plt.plot(fpr, tpr, color=color, lw=lw,
                 label='ROC area = {0:0.2f}'.format(roc_auc))
        plt.plot([0, 1], [0, 1], 'k--', lw=lw)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(t)
        plt.legend(loc="lower right")
        plt.savefig(fn + '_roc.png')

    return roc_auc, pr_auc


def multiclass_metrics(classes, y_truth, y_preds, fn=None, t=None, plot=True):  # assumes up to 4 classes
    """Get area under the curve for roc and precission-recall, for each class and micro- and macro-averages
    can also plot roc

    Args:
        classes: str ndarray of class names
        y_truth: int ndarray (sample_size, num_classes) of true classes, one-hot encoded
        y_preds: float ndarray (sample_size, num_classes) of true classes
        fn: file name, full path
        t: plot title
        plot: whether to plot roc

    Returns: auc of roc and precission-recall

    """
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    pr_auc = dict()
    precision = dict()
    recall = dict()
    for i in range(len(classes)):
        n = classes[i]
        fpr[n], tpr[n], _ = roc_curve(y_truth[:, i], y_preds[:, i])
        roc_auc[n] = auc(fpr[n], tpr[n])
        precision[n], recall[n], _ = precision_recall_curve(y_truth[:, i], y_preds[:, i])
        pr_auc[n] = auc(recall[n], precision[n])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_truth.ravel(), y_preds.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # A "micro-average": quantifying score on all classes jointly
    precision["micro"], recall["micro"], _ = precision_recall_curve(y_truth.ravel(), y_preds.ravel())
    pr_auc["micro"] = auc(recall["micro"], precision["micro"])

    # tag on f1_score
    roc_auc["f1_micro"] = f1_score(y_true=np.array(y_truth).argmax(axis=-1),
                                   y_pred=np.array(y_preds).argmax(axis=-1),
                                   average='micro')

    # Compute macro-average ROC curve and ROC area
    # macro-average ROC curves (average per class in a 1-vs-all fashion)
    # micro-averaged ROC curves (consider all positives and negatives together as single class)

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[classes[i]] for i in range(len(classes))]))
    all_recall = np.unique(np.concatenate([recall[classes[i]] for i in range(len(classes))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(len(classes)):
        mean_tpr += interp(all_fpr, fpr[classes[i]], tpr[classes[i]])

    mean_precision = np.zeros_like(all_recall)
    for i in range(len(classes)):
        mean_precision += interp(all_recall, recall[classes[i]], precision[classes[i]])

    # Finally average it and compute AUC
    mean_tpr /= len(classes)
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    mean_precision /= len(classes)
    recall["macro"] = all_recall
    precision["macro"] = mean_precision
    pr_auc["macro"] = auc(recall["macro"], precision["macro"])

    if plot:
        colors = ['aqua', 'darkorange', 'cornflowerblue', 'orangered']
        lw=2

        # plot roc
        plt.figure()
        plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["micro"]),
             color='deeppink', linestyle=':', linewidth=4)

        plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["macro"]),
             color='navy', linestyle=':', linewidth=4)

        for i, color, name in zip(range(len(classes)), colors, classes):
            i = classes[i]
            plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                     label='ROC curve of class {0} (area = {1:0.2f})'
                     ''.format(name, roc_auc[i]))

        plt.plot([0, 1], [0, 1], 'k--', lw=lw)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(t)
        plt.legend(loc="lower right")
        plt.savefig(fn + '_roc.png')
    plt.close('all')

    return roc_auc, pr_auc
