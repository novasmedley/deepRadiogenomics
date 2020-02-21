"""
Methods to get data for modeling and saving experiments
"""
import os
import sys
import pandas as pd


def float_formatter():
    return lambda x: "%.3f" % x


def get_param_str(args):
    param_str = ( str(args['h1']) + '_' + str(args['h2']) + '_' + str(args['h3']) + '_' + str(args['h4']) + '_'
                  + str(args['h5']) + '_' + args['act']
                  + '_decay_' + str(args['decay'])
                  + '_drop_' + str(args['drop'])
                  + '_opt_' + str(args['opt'])
                  + '_loss_' + args['loss']
                  + '_bat_' + str(args['batch'])
                  + '_eph_' + str(args['epoch'])
                  )
    if args['pretrain'] is not None:
        param_str += '_ae_nl_' + str(args['num_ae_layers'])
        param_str += '_freeze_' + str(args['freeze'])
    return param_str


'''
GBM - TCGA-TCIA dataset
'''


def get_gbm_data(data_dir, data_type, label_name):
    """

    Args:
        data_dir: path to model data
        data_type:
        label_name:

    Returns:

    """
    gene = os.path.join(data_dir, 'gene_expression.txt')
    mr_vasari = os.path.join(data_dir, 'vasari_annotations.csv')
    subtype = os.path.join(data_dir, 'TCGA_unified_CORE_ClaNC840.txt')  # verhaak subtypes
    ignore = None

    # get outputs
    if data_type == 'vasari':
        y = pd.read_table(mr_vasari, index_col=0, sep=',')
        y = y.drop(columns=['research.id', 'scan.type', 'comments', 'study.date'])  # TODO include tumor location information
        ignore = ['n/a', 'indeterminate']  # types of categories to ignore

    elif data_type == 'verhaak':
        y = pd.read_table(subtype)
        y = y.iloc[0:1, 2:].transpose()
        y.index = [i[0:12] for i in y.index]
        y.columns = ['subtype']

    else:
        print('incorrect datatype')
        return None

    # select outcomes
    if label_name != 'all':
        these = [label_name in i for i in y.columns.values]
        if these is None:
            sys.exit("label, " + label_name + " not found in dataset")

        output_names = y.columns[these]
        y = y[output_names]

    y = y.dropna()  # remove patients with missing data

    # ignore and merge certain outcomes
    if ignore is not None:
        for i in ignore:
            y = y[(y != i).all(axis=1)]

    # binarize
    bin_these = ['f5.proportion.enhancing', 'f6.proportion.ncet', 'f7.proportion.necrosis', 'f14.proportion.edema']
    for i in bin_these:
        if i in y.columns.values:
            y.loc[y[i] == 'Minimal'] = 'Less than 1/3'
            y.loc[y[i] == 'More than 2/3'] = 'More than 1/3'
            y.loc[y[i] == 'Between 1/3 and 2/3'] = 'More than 1/3'

    if 'f9.multifocal.or.multicentric' in y.columns.values:
        y.loc[y['f9.multifocal.or.multicentric'] != 'Focal'] = 'non-focal'

    if 'f10.t1.flair.ratio' in y.columns.values:
        y.loc[y['f10.t1.flair.ratio'] != 'expansive (T1~FLAIR)'] = 'infiltrative'

    # get inputs
    genes = pd.read_table(gene)
    
    # get cohort
    cohort = genes.index.intersection(y.index)
    genes = genes.loc[cohort]
    y = y.loc[cohort]
    print(genes.shape)
    print(y.shape)
    print(y.columns)
    print('indices still match?: ', genes.index.equals(y.index))  # check row names match
    ids = genes.index.values

    return genes, y, ids