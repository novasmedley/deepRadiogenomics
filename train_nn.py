"""
Train neural networks
"""
from neuralnet import *
from setup import *
from bootstrap import *
import time
from keras.utils import np_utils
from tensorflow import set_random_seed


def main():
    start_time = time.time()
    model_type = 'nn'  # neural network, not autoencoder (ae)

    args = get_args()
    args['num_boots'] = 500  # bootstrap seed generation, constant

    np.random.seed(args['seed'])
    set_random_seed(args['seed'])

    if args['bs_iter'] is not None:  # boostrapping?
        print('bootstrapping')
        args['exp'] = args['exp'] + '/boot_' + str(args['bs_iter'])

    # experiment params and setup
    args['exp_dir'], args['tb_dir'], args['log_file'], \
        args['scale_in'], args['scale_out'], args['param_str'] = setup_experiment(args, overwrite=True)
    print('running experiment: ' + args['param_str'])

    # prep data and set results metric names
    genes, y_labels, ids = get_gbm_data(data_dir=args['data'],
                                        data_type=args['dataType'],
                                        label_name=args['label'])

    if args['pretrain'] is not None:  # use pretrained autoencoder?
        print('pretraining')
        args['ae_model'] = os.path.join(args['pretrain'], 'model_retrained.h5')  # get ae
        args['ae_scaler'] = os.path.join(args['pretrain'], 'geneTrainScalers.pkl')  # get ae's preprocessing

    genes = genes.values
    y_names = y_labels.columns
    y_labels = y_labels.values

    if args['predType'] == 'regression':
        y = y_labels
        y_labels = None
        classes = None
        args['nonneg'] = get_nonneg(y)
        res_colnames, res_rename = get_regression_names()
    else:
        y, le, classes = encode_y(y_labels)
        y_labels = y_labels.reshape(y_labels.shape[0], )  # reshape for compute_sample_weight

        pickle.dump(le, open(args['exp_dir']+'/labelencoder.pkl', 'wb'))

        if args['predType'] == 'multiClass':
            print('multiclass')
            y = np_utils.to_categorical(y)  # ints to one-hot encoding

        for i in classes:
            print(list(y_labels).count(i), '\t', i)
        res_colnames = None

    # TRAIN
    if args['retrain']:
        print('retraining')
        args['save'] = True
        folds = None
        args['model_name'] = 'model_retrained.h5'
        args['logger'] = 'logger_retrained.csv'

        model, t_metrics, _, train_preds, _ = fit_model(args,
                                                        x_train=genes,
                                                        y_train=y,
                                                        x_val=None,
                                                        y_val=None,
                                                        sample_weights=None)
        # save retrain results
        if args['predType'] == 'regression':
            # save individual feature metrics
            t_metrics = np.array(t_metrics)
            t_avg = np.mean(t_metrics, axis=0)
            y_names.append('avg')
            results = pd.DataFrame(data=np.hstack((t_metrics, t_avg)),
                                   columns=res_colnames[0: (len(t_avg)/2)],
                                   index=y_names)
            results.to_csv(args['exp_dir'] + '/retrain_scores.txt', sep='\t')
            cnames = y_names
        else:
            # get roc plots
            t_plot = os.path.join(args['exp_dir'], 'retrain')
            if args['predType'] == 'multiClass':
                t_roc, _ = multiclass_metrics(classes,
                                              y_truth=y,
                                              y_preds=train_preds,
                                              fn=t_plot,
                                              t=args['label'])
                cnames = classes
            else:
                t_roc, _ = binaryclass_metrics(y_truth=y,
                                               y_preds=train_preds,
                                               fn=t_plot,
                                               t=args['label'])
                cnames = [args['label']]

        # save predictions
        fn = os.path.join(args['exp_dir'], 'retrain_preds.csv')
        pd.DataFrame(data=train_preds, columns=cnames, index=ids).to_csv(fn)

    else:
        print('cross_validation')
        args['save'] = False  # don't save individual fold predictions and metrics, just get averaged performances

        if args['bs_iter'] is None:
            # cv
            cv = get_split(x=genes,
                           y=y,
                           y_labels=y_labels,
                           pred_type=args['predType'],
                           seed=args['seed'],
                           folds=args['folds'])
        else:
            # bootstrap cv
            seeds = get_seeds(seed=args['seed'], num_seeds=args['num_boots'])
            cv = get_split(x=genes,
                           y=y,
                           y_labels=y_labels,
                           folds=args['folds'],
                           seed=seeds[args['bs_iter']],  # use another split seed
                           pred_type=args['predType'])

            # call a bootstrap generator
            if args['bs_method'] == 1:
                cv = bootstap_gen_cv(cv_split=cv,
                                     seed=seeds[args['bs_iter']],
                                     y=y,
                                     classes=classes)
            if args['bs_method'] == 2:
                cv = bootstap_gen_cv_class(cv_split=cv,
                                           seed=seeds[args['bs_iter']],
                                           y=y,
                                           folds=args['folds'])

        cv_logger, cv_metrics = cross_validate(args,
                                               x=genes,
                                               y=y,
                                               cv=cv,
                                               model_type=model_type,
                                               y_labels=y_labels,
                                               y_names=y_names,
                                               ids=ids,
                                               classes=classes)

        cv_logger.to_csv(args['exp_dir'] + '/cv_logger.txt', sep='\t')
        if cv_metrics is not None:
            cv_metrics.to_csv(args['exp_dir'] + '/cv_metrics.txt', sep='\t')

    # logging
    elapsed_time = time.time() - start_time

    with open(args['log_file'], 'w') as f:
        print(str(datetime.now()), file=f)
        print('\n', file=f)
        print('x shape:    \t', genes.shape, file=f)
        print('y label:    \t', args['label'], file=f)
        print('\n', file=f)
        print('param str \t', args['param_str'], file=f)
        print('patience \t', args['patience'], file=f)
        print('folds   \t', args['folds'], file=f)
        print('retrain \t', args['retrain'], file=f)
        print('pretrain \t', args['pretrain'], file=f)
        print('seed    \t', args['seed'], file=f)
        print('\n', file=f)
        print('tot secs \t', elapsed_time, file=f)
        if args['predType'] != 'regression':
            for i in classes:
                print(list(y_labels).count(i), '\t', i, file=f)


if __name__ == "__main__":
    main()
