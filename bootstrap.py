"""
Methods for boostrapping model performance
"""

import random
import numpy as np
import pandas as pd
from ast import literal_eval
from sklearn.utils import resample


def get_seeds(seed, num_seeds):
    random.seed(seed)
    return random.sample(range(0, 2 ** 32), num_seeds)


def get_resample_index(seed, indices):
    return resample(indices, n_samples=len(indices), random_state=seed)


def get_resample(genes, y_labels, ids, seed, num_seeds, bs_iter):
    seeds = get_seeds(seed=seed,  num_seeds=num_seeds)
    indices = get_resample_index(seed=seeds[bs_iter],
                                 indices=ids)
    return genes.loc[indices], y_labels.loc[indices], indices


def get_selected_params(file, model, label):
    cv = pd.read_csv(file, sep=',', header=0, index_col=0)
    params = cv['params'].loc[(cv['label'] == label) & (cv['model'] == model)].values
    params = literal_eval(params[0])  # to dict
    for k, v in params.items():  # put each param value in a list
        params[k] = [v]

    return params


def bootstap_gen_cv(cv_split, seed, y, classes):
    """f or each cv split, resample from train and validation partitions

    Args:
        cv_split: fold indices generator
        seed: int seed
        y: int ndarray y class labels
        classes: str ndarray y class names

    Returns: a bootstrap generator, for input to gridsearchcv

    """

    for train_index, val_index in cv_split:
        train = resample(train_index, n_samples=len(train_index), random_state=seed)
        if classes is None:
            val = resample(val_index, n_samples=len(val_index), random_state=seed)
        else:
            # if classes are highly imbalanced, make sure they appear in validation still
            n = 0
            s = seed
            while n < len(np.unique(y)):
                val = resample(val_index, n_samples=len(val_index), random_state=s)
                n = len(np.unique(y[val]))
                s += 1
        yield train, val


def bootstap_gen_cv_class(cv_split, seed, y, folds):
    """for each cv split, for each partition, separate classes and sample from each class

    Args:
        cv_split: fold indices generator
        seed: int seed
        y: int ndarray y class labels
        folds: int folds

    Returns: a bootstrap generator, for input to gridsearchcv

    """
    seeds = get_seeds(seed=seed, num_seeds=folds)
    i = 0
    for train_index, val_index in cv_split:
        s = seeds[i]
        train = get_boot_idx(train_index, y, s)
        val = get_boot_idx(val_index, y, s)
        i += 1
        yield train, val


def get_boot_idx(index, y, seed):
    """for each class in subset, sample and return indices

    Args:
        index: subset of y_labels via indices
        y: class labels
        seed: int seed

    Returns:

    """

    sub_y = y[index]  # subset
    classes, counts = np.unique(sub_y, return_counts=True)  # subset class and counts
    bs = []
    for c, n, in zip(classes, counts):
        # print(c, n)
        a = resample(index[sub_y == c], n_samples=n, random_state=seed)
        bs.extend(a)
    return np.array(bs)
