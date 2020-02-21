'''
train other model(s)
'''

import sys
import shutil
import logging
from sklearn.externals import joblib

from other_models_utils import *
from setup import get_gbm_data
from neuralnet import get_split, encode_y
from bootstrap import *


def main():
    s = time.time()

    args = get_args()
    args['seed'] = 425  # constants
    args['num_boots'] = 500  # bootstrap num samples

    # boostrapping?
    if args['bs_iter'] is not None:
        args['exp'] = args['exp'] + '/boot_' + str(args['bs_iter'])

    # setup experiment output
    print(args['exp'])
    print(args['model'])

    args['exp'] = os.path.join(args['dir'], args['exp'], args['model'], args['label'])
    if os.path.exists(args['exp']):
        shutil.rmtree(args['exp'])
    os.makedirs(args['exp'])
    print(args['exp'])

    # get data
    genes, y, ids = get_gbm_data(data_dir=args['data'],
                                 data_type=args['dataType'],
                                 label_name=args['label'])

    # format data
    binary = False
    pos_weight = 1

    if args['predType'] == 'regression':
        # format radiomics
        if args['dataType'] == 'radiomic':
            img_scaler = get_img_scaler(y)
            fn = os.path.join(args['exp'], 'imgTrainScalers.csv')
            if args['bs_iter'] is None:
                img_scales = pd.DataFrame(img_scaler.max_abs_, index=y.columns.values)  # save for later
                img_scales.to_csv(fn)
            y = y.values
            y = img_scaler.transform(y)  # apply transformer
    else:
        y = y.values
        if args['predType'] == 'binaryClass':
            binary = True

        y_labels, le, classes = encode_y(y)  # convert to 0's and 1's for xgboost to calc sample weighting
        pos_weight = float(np.sum(y_labels == 0)) / np.sum(y_labels == 1)
        # no pos_weight definition for multiclass in XGBClassifier

    y = y.ravel()  # only one col

    # modeling

    # param search
    n, m, p = get_model_params(args['model'], binary=binary, scale_pos_weight=pos_weight)

    if args['bs_iter'] is None:
        # cv
        cv = get_split(x=genes,
                       y=y,
                       y_labels=None,
                       folds=args['folds'],
                       seed=args['seed'],
                       pred_type=args['predType'])
    else:
        # bootstrap
        seeds = get_seeds(seed=args['seed'], num_seeds=args['num_boots'])
        cv = get_split(x=genes,
                       y=y,
                       y_labels=None,
                       folds=args['folds'],
                       seed=seeds[args['bs_iter']],  # use another split seed
                       pred_type=args['predType'])

        # reset cv generator to one that bootstraps based classes for each cv split

        # BS METHOD 1
        if args['bs_method'] == 1:
            cv = bootstap_gen_cv(cv_split=cv,
                                 seed=seeds[args['bs_iter']],
                                 y=y,
                                 classes=classes)
        # BS METHOD 2
        if args['bs_method'] == 2:
            cv = bootstap_gen_cv_class(cv_split=cv,
                                       seed=seeds[args['bs_iter']],
                                       y=y,
                                       folds=args['folds'])

        p = get_selected_params(file=args['params'],  # reset to a single param
                                model=args['model'],
                                label=args['label'])

    cv_grid, t = cross_val(x=genes,
                           y=y,
                           model_name=n,
                           model=m,
                           param_grid=p,
                           cv=cv,
                           pred_type=args['predType'],
                           n_jobs=args['cpus'])
    end = (time.time() - s) / 60
    print(end)

    # save experiment
    res_fn = os.path.join(args['exp'], 'cv_results.csv')
    pd.DataFrame(cv_grid.cv_results_).to_csv(res_fn, sep=',')  # parameter tuning

    if args['bs_iter'] is None:
        mod_fn = os.path.join(args['exp'], 'cv_best_estimator.pkl')
        joblib.dump(cv_grid.best_estimator_, mod_fn, compress=1)  # model

    # log
    logging.basicConfig(filename=os.path.join(args['exp'], 'exp.log'), level=logging.INFO)
    logging.info('model: %s', args['model'])
    logging.info('label: %s', args['label'])
    logging.info('gene shape: (%d, %d)', genes.shape[0], genes.shape[1])
    logging.info('label shape: %d', y.shape[0])
    logging.info('exp min: %.3f', end)
    # logging.info('boot iteration: %f', args['bs_iter'])


if __name__ == "__main__":
    main()
