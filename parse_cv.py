"""
Parse cv results
"""

import os
import re
import time
import pandas as pd
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', help='full path to experiment directory, where grid search results were saved')
    parser.add_argument('--model', help='neural networks or other models?', choices=['nn', 'other'], default='nn')

    args = parser.parse_args()
    args = vars(args)  # convert to dictionary
    return args


def get_val_loss(log, fn):
    mi = log['val_loss'].idxmin(axis=1)
    log = log.loc[mi:mi]
    log.index = [fn]
    return log


def read_cv_scores(data_dir, label_name, exp_name, fn):
    p = os.path.join(data_dir, label_name, 'neuralnets', exp_name, fn)
    s = pd.read_csv(p, sep='\t', header=0, index_col=0)
    s = s.loc[['avg']]
    s.index = [exp_name]
    s['label'] = label_name
    return s


def read_cv_scores_other(data_dir, model_name, label_name, fn):
    p = os.path.join(data_dir, model_name, label_name, fn)
    s = pd.read_csv(p, sep=',', header=0, index_col=0)
    cols = [ i for i in s.columns if i in ['params', 'split0_test_score',
       'split1_test_score', 'split2_test_score', 'split3_test_score',
       'split4_test_score', 'split5_test_score', 'split6_test_score',
       'split7_test_score', 'split8_test_score', 'split9_test_score',
       'mean_test_score', 'std_test_score', 'rank_test_score', 'mean_fit_time']]
    s = s[cols]
    s = s.sort_values(by=['rank_test_score'], ascending=True)
    s = s
    s['label'] = label_name
    s['model'] = model_name
    s = s.loc[s['rank_test_score'] == 1]
    s = s.loc[[s['mean_fit_time'].idxmin()]]

    return s


def parse_param_str(param_str):
    archit = re.findall('^\d*_\d*_\d*_\d*_\d*', param_str)[0]
    act = re.findall('sigmoid|relu|tanh', param_str)[0]
    drop = re.findall('drop_(.*?)_', param_str)[0]
    opt = re.findall('opt_(.*?)_', param_str)[0]
    loss = re.findall('loss_(.*?)_', param_str)[0]
    nl = re.findall('nl_(.*?)$', param_str)

    return [archit, act, drop, opt, loss, nl]


def get_cv_scores(data_dir, fn='cv_logger.txt', these_labs=None):
    # fn = 'cv_logger.txt' # file to parse
    labels = [i for i in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, i))]  # sub directories
    print(labels)
    if these_labs is not None:
        labels = these_labs
        print(labels)
    cv = [read_cv_scores(data_dir, i, j, fn) for i in labels for j in
          os.listdir(os.path.join(data_dir, i, 'neuralnets'))]
    cv = pd.concat(cv)

    a = pd.DataFrame([parse_param_str(i) for i in cv.index])
    a.columns = ['archit', 'act', 'drop', 'opt', 'loss', 'nl']
    a.index = cv.index
    cv = pd.concat([cv, a], axis=1)
    cv = cv.reset_index(drop=True)

    return cv


def main():
    print('reading labels...')
    start_time = time.time()
    args = get_args()

    if args['model'] == 'nn':
        cv_scores = get_cv_scores(data_dir=args['dir'])
    else:
        om = [i for i in os.listdir(args['dir']) if os.path.isdir(os.path.join(args['dir'], i))]
        labels = os.listdir(os.path.join(args['dir'], om[0]))
        cv_scores = [read_cv_scores_other(args['dir'], i, j, 'cv_results.csv') for i in om for j in labels]
        cv_scores = pd.concat(cv_scores)

    cv_scores.to_csv(os.path.join(args['dir'], 'parsed_cv_results_test.csv'))
    print(time.time() - start_time)


if __name__ == "__main__":
    main()
