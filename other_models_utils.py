"""
Helper functions for training other types of classiers
"""
import os
import time
import argparse
import numpy as np
import pandas as pd

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, MaxAbsScaler
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_absolute_error, f1_score, make_scorer

import xgboost as xgb
from sklearn.svm import SVR, SVC
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.linear_model import ElasticNet, LogisticRegression
from keras.utils import np_utils

from neuralnet import multiclass_metrics


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataType', help='output data type', choices=['radiomic',
                                                                        'histology', 'stage',
                                                                        'vasari', 'verhaak'], default='radiomic')
    parser.add_argument('--predType', help='prediction type', choices=['regression', 'binaryClass', 'multiClass'],
                        default='regression')
    parser.add_argument('--exp', help='experiment name')
    parser.add_argument('--dir', help='full path to experiment directory, where all data will be saved')
    parser.add_argument('--data', help='full path to data directory')
    parser.add_argument('--label', help='name of output to find, must be at least a partial match in data files')
    parser.add_argument('--model', help='model',
                        choices=['rfr', 'lr', 'gbr', 'svmr',
                                 'rfc', 'logit1', 'logit2', 'gbc', 'svmc'])
    parser.add_argument('--folds', help='number of cross validation folds', type=int, default=10)
    parser.add_argument('--bs_iter', help='boostrap ith iteration', type=int, default=None)
    parser.add_argument('--params', help='parsed cv results with params for model and label', type=str, default=None)
    parser.add_argument('--bs_method', help='boostrap method', type=int, choices=[1, 2], default=1)
    parser.add_argument('--cpus', help='number of cpus', type=int, default=1)

    args = parser.parse_args()
    args = vars(args)  # convert to dictionary
    return args


def get_model_params(model_name, binary=True, nthread=1, scale_pos_weight=None):
    return {
        'rfr': get_rf_r(),
        'lr': get_lr(),
        'gbr': get_gbr(),
        'svmr': get_svm_r(),

        'rfc': get_rf_c(),
        'logit1': get_logit(l1=True, binary=binary),
        'logit2': get_logit(l1=False, binary=binary),
        'gbc': get_gbc(binary=binary, nthread=nthread, scale_pos_weight=scale_pos_weight),
        'svmc': get_svm_c()
    }[model_name]


def get_img_scaler(img_features):
    scaler = MaxAbsScaler()
    scaler.fit(img_features)
    return scaler


def get_pipe(name, model):
    return Pipeline(steps=[('preprocess', StandardScaler()), (name, model)])


def get_multiclass_scorer():
    return make_scorer(f1_score, greater_is_better=True, average='micro')


# CV with grid search and refitting
def cross_val(x, y, model_name, model, param_grid, cv, pred_type, n_jobs=1):
    if pred_type == 'regression':
        scorer = 'r2'
    else:
        scorer = 'roc_auc'  # binary
        if pred_type == 'multiClass':
            scorer = get_multiclass_scorer()  # gridsearchCV doesn't do multiclass, to use make_scorer()

    pipe = get_pipe(model_name, model)
    s = time.time()
    clf = GridSearchCV(estimator=pipe,
                       param_grid=param_grid,
                       cv=cv,
                       scoring=scorer,
                       return_train_score=True,
                       n_jobs=n_jobs,
                       verbose=1)
    clf.fit(X=x, y=y)
    t = (time.time() - s)

    return clf, t


def get_results(args, model, x, y, y_one_hot=None, fit=False, classes=None, labelencoder=None):
    # if multiclass, y is expected to be one-hot label
    if fit:
        model.fit(X=x, y=y)

    pred = model.predict(x)  # classes as ints
    score = model.score(X=x, y=y)  # R^2 in regression, mean acc in classification

    if classes is None:
        mae = mean_absolute_error(y_true=y, y_pred=pred)
        preds = pd.DataFrame(data={'y_pred': pred,
                                   'y_true': y})
    else:
        mae = None
        if len(classes) == 2:
            preds = pd.DataFrame(data={'y_pred': pred,
                                       'y_true': y})
            get_roc_results(y_preds=preds,
                            y_truth=y,
                            plotname='train_',
                            classes=classes,
                            label=args['label'],
                            exp_folder=args['exp'])
        else:
            # multiclass
            preds = pd.DataFrame(data={'y_pred': pred,
                                       'y_true': y,
                                       'y_pred_class': labelencoder.inverse_transform(pred),  # reverse to class name
                                       'y_true_class': labelencoder.inverse_transform(y_one_hot)})
            get_roc_results(y_preds=np_utils.to_categorical(preds['y_pred'], num_classes=len(classes)),
                            y_truth=y_one_hot,
                            plotname='train_',
                            classes=classes,
                            label=args['label'],
                            exp_folder=args['exp'])

    return model, preds, mae, score


def get_roc_results(y_preds, y_truth, plotname, classes, label, exp_folder):
    roc_scores, pr_scores = multiclass_metrics(classes=classes,
                                               y_truth=y_truth,
                                               y_preds=y_preds,
                                               fn=os.path.join(exp_folder, plotname),
                                               t=label)

    roc_scores = pd.DataFrame(roc_scores, index=[0])
    roc_scores.columns.values[0:3] = classes
    roc_scores.to_csv(os.path.join(exp_folder, 'roc_scores_' + plotname + '.csv'))


def get_rf_r():
    n = 'RandomForest'
    m = RandomForestRegressor()
    p = {'RandomForest__max_features': ['sqrt', 'log2', None],  # features
         'RandomForest__criterion': ['mae', 'mse'],
         'RandomForest__max_depth': [None],
         'RandomForest__n_estimators': list(np.arange(50, 2050, 50))  # 40 trees
         }
    return n, m, p


def get_rf_c():
    n = 'RandomForest'
    m = RandomForestClassifier()
    p = {'RandomForest__max_features': ['sqrt', 'log2', None],  # features
         'RandomForest__criterion': ['gini', 'entropy'],
         'RandomForest__max_depth': [None],
         'RandomForest__n_estimators': list(np.arange(50, 2050, 50)),  # 40 trees
         'RandomForest__class_weight': ['balanced']
         }
    return n, m, p


def get_lr():
    n = 'LR'
    m = ElasticNet()
    lmax = 2
    p = {'LR__alpha': list(np.exp(np.linspace(np.log(0.001), np.log(lmax), 100))),  # lambda in glmnet
         'LR__l1_ratio': list(np.arange(0.1, 1.1, 0.1)),  # alpha in glmnet
         'LR__max_iter': [2000]
         }
    return n, m, p


def get_logit(binary=True, l1=True):
    n = 'Logit'
    m = LogisticRegression()
    p = {'Logit__class_weight': ['balanced'],
         'Logit__C': list(np.logspace(-3, 3, 1000)),
         'Logit__max_iter': [2000]
         }

    if binary:
        if l1:
            p['Logit__solver'] = ['liblinear']
        else:
            p['Logit__solver'] = ['lbfgs', 'sag', 'newton-cg']

    else:
        p['Logit__multi_class'] = ['multinomial']
        if l1:
            p['Logit__solver'] = ['saga']
        else:
            p['Logit__solver'] = ['lbfgs', 'sag', 'newton-cg']

    return n, m, p


def get_gbr(nthread=1):
    n = 'GBR'
    m = xgb.XGBRegressor()
    lmax = 2
    p = {'GBR__max_depth': list(np.arange(1, 5, 1)),  #4 , more = overfit
         'GBR__n_estimators': list(np.arange(50, 1050, 50)),  #20
         'GBR__booster': ['gbtree'],
         'GBR__learning_rate': list(np.arange(0.01, 0.51, 0.05)) # 10
         }
    return n, m, p


def get_gbc(binary=True, nthread=1, scale_pos_weight=1):
    n = 'GBC'
    m = xgb.XGBClassifier()
    lmax = 2
    p = {'GBC__max_depth': list(np.arange(1, 5, 1)),   #4 , more = overfit
         'GBC__n_estimators': list(np.arange(50, 1000, 50)),  #20, trees
         'GBC__booster': ['gbtree'],
         'GBC__gamma': [0.01],
         'GBC__learning_rate': list(np.arange(0.01, 0.51, 0.05)),
         'GBC__nthread': [nthread],
         'GBC__scale_pos_weight': [scale_pos_weight]
         }
    if binary:
        p['GBC__objective']: ['binary:logistic']
    else:
        p['GBC__objective']: ['multi:softmax']
    return n, m, p


def get_svm_r():
    n = 'SVM'
    m = SVR()
    p = {'SVM__C': list(np.logspace(-6, 1, 1000)),
         'SVM__kernel': ['linear', 'poly', 'rbf', 'sigmoid'],
         'SVM__cache_size': [1000]
         }
    return n, m, p


def get_svm_c():
    n = 'SVM'
    m = SVC()
    p = {'SVM__C': list(np.logspace(-6, 6, 2000)),
         'SVM__kernel': ['linear', 'poly', 'rbf', 'sigmoid'],
         'SVM__cache_size': [1000],
         'SVM__decision_function_shape': ['ovo'],
         'SVM__class_weight': ['balanced']
        }
    return n, m, p
