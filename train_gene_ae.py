'''
Create a deep autoencoder using gene expression as inputs.
'''

import sys
sys.path.insert(0, '../')

from neuralnet import *
from setup import *
import time


def main():
    start_time = time.time()
    model_type = 'ae'

    seed = 504  # constants for cross-validation

    args = get_args()
    args['exp_dir'], args['tb_dir'], args['log_file'], \
        args['scale_in'], args['scale_out'], args['param_str'] = setup_experiment(args, overwrite=True)

    # load ae data
    genes = pd.read_table(os.path.join(args['data'], 'gene_expression.txt'))

    # remove patients with vasari data
    v = os.path.join(args['data'], 'vasari_annotations.csv')
    v = pd.read_table(v, index_col=0, sep=',')
    v = v.drop(columns=['research.id', 'scan.type', 'comments', 'study.date'])

    genes = genes[~genes.index.isin(v.index)]
    genes = genes.values
    print(genes.shape)
    del v

    args['nonneg'] = False  # genes are mean centered and std normalized

    # TRAIN
    if args['retrain']:
        print('retraining')
        args['save'] = True
        args['model_name'] = 'model_retrained.h5'
        args['logger'] = 'logger_retrained.csv'
        ae, train_metrics, _, train_preds, _ = fit_autoencoder(args,
                                                               x_train=genes,
                                                               x_val=None)
        df = pd.DataFrame(train_metrics)
        df.columns = ['mae', 'mape', 'r2']
        df.to_csv(args['exp_dir'] + '/retrain_metrics.txt', sep='\t')

        del ae
        K.clear_session()
    else:
        print('cross_validation')
        args['save'] = False
        cv = get_split(x=genes,
                       y=genes,
                       y_labels=None,
                       folds=args['folds'],
                       pred_type=args['predType'],
                       seed=seed)
        cv_log, cv_metrics = cross_validate(args,
                                            x=genes,
                                            y=genes,
                                            cv=cv,
                                            model_type=model_type)

        cv_log.to_csv(args['exp_dir'] + '/cv_logger.txt', sep='\t')
        if cv_metrics is not None:
            cv_metrics.to_csv(args['exp_dir'] + '/cv_metrics.txt', sep='\t')

    # logging
    elapsed_time = time.time() - start_time
    with open(args['log_file'], 'w') as f:
        print(str(datetime.now()), file=f)
        print('\n', file=f)
        print('x shape:    \t', genes.shape, file=f)
        print('\n', file=f)
        print('param str \t', args['param_str'], file=f)
        print('patience \t', args['patience'], file=f)
        print('folds   \t', args['folds'], file=f)
        print('retrain \t', args['retrain'], file=f)
        print('pretrain \t', args['pretrain'], file=f)
        print('seed    \t', seed, file=f)
        print('\n', file=f)
        print('tot secs \t', elapsed_time, file=f)


if __name__ == "__main__":
    main()
