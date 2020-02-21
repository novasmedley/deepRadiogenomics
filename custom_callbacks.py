"""
Get roc auc or r2 as callback metric during training
Save roc auc and average precision in logger
"""

from sklearn.metrics import roc_auc_score, average_precision_score
from keras.callbacks import Callback


class roc_callback(Callback):
    # see, https://github.com/keras-team/keras/issues/3230#issuecomment-319208366
    def __init__(self, training_data, validation_data, binary):
        self.x = training_data[0]
        self.y = training_data[1]
        self.x_val = validation_data[0]
        self.y_val = validation_data[1]
        self.binary = binary

    def on_train_begin(self, logs={}):
        return

    def on_train_end(self, logs={}):
        return

    def on_epoch_begin(self, epoch, logs={}):
        return

    def on_epoch_end(self, epoch, logs={}):
        y_pred = self.model.predict(self.x)

        # add metrics to logger
        if self.binary:
            avg = 'macro'
        else:
            avg = 'micro'

        logs['roc'] = roc_auc_score(self.y, y_pred, average=avg)
        logs['pr'] = average_precision_score(y_true=self.y, y_score=y_pred, average=avg)

        if self.x_val is not None:
            y_pred_val = self.model.predict(self.x_val)
            logs['val_roc'] = roc_auc_score(self.y_val, y_pred_val, average=avg)
            logs['val_pr'] = average_precision_score(y_true=self.y_val, y_score=y_pred_val, average=avg)

        if self.x_val is not None:
            print('\rroc: %s - val_roc: %s' % (str(round(logs['roc'], 4)), str(round(logs['val_roc'], 4))), end=100 * ' ' + '\n')
        else:
            print('\rroc: %s ' % (str(round(logs['roc'], 4))), end=100 * ' ' + '\n')

        return

    def on_batch_begin(self, batch, logs={}):
        return

    def on_batch_end(self, batch, logs={}):
        return


def get_monitor_info(prediction_type, val, patience=200):
    if prediction_type == 'regression':
        if val is None:
            vm = 'r2'  # just monitor training
        else:
            vm = 'val_r2'
        vmode = 'max'
        min_d = 0.005
        p = patience
    else:
        if val is None:
            vm = 'roc'  # just monitor training
        else:
            vm = 'val_roc'
        vmode = 'max'
        min_d = 0.005
        p = patience
    return vm, vmode, min_d, p
