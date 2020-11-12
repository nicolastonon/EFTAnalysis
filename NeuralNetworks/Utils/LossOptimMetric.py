#Define the loss function, optimizer and metric to be used to train the network

import tensorflow as tf
from tensorflow.keras.optimizers import SGD, Adam,Nadam, RMSprop, Adadelta
from tensorflow.keras import backend as K
from Utils.ColoredPrintout import colors
from Utils.adabound_tf import AdaBound
# from focal_loss import BinaryFocalLoss #Requires installation

# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
 #######  ########  ######## #### ##     ## #### ######## ######## ########
##     ## ##     ##    ##     ##  ###   ###  ##       ##  ##       ##     ##
##     ## ##     ##    ##     ##  #### ####  ##      ##   ##       ##     ##
##     ## ########     ##     ##  ## ### ##  ##     ##    ######   ########
##     ## ##           ##     ##  ##     ##  ##    ##     ##       ##   ##
##     ## ##           ##     ##  ##     ##  ##   ##      ##       ##    ##
 #######  ##           ##    #### ##     ## #### ######## ######## ##     ##
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Focal loss definition, taken from: https://github.com/mkocabas/focal-loss-keras/blob/master/focal_loss.py, based on paper 'Focal Loss for Dense Object Detection' (Lin et. al.)
'''
def focal_loss(gamma=2., alpha=.25):
    def focal_loss_fixed(y_true, y_pred):
        pt_1 = tf.where(tf.equal(y_true, 1), y_pred, tf.ones_like(y_pred))
        pt_0 = tf.where(tf.equal(y_true, 0), y_pred, tf.zeros_like(y_pred))
        return -K.mean(alpha * K.pow(1. - pt_1, gamma) * K.log(pt_1)) - K.mean((1 - alpha) * K.pow(pt_0, gamma) * K.log(1. - pt_0))

    return focal_loss_fixed
'''


# //--------------------------------------------
# //--------------------------------------------


#See: https://www.tensorflow.org/guide/keras/train_and_evaluate
# mean squared error, but with an added term that will de-incentivize prediction values far from 0.5 --> Prevents network from being 'over-confident'
# Example: model.compile(optimizer=keras.optimizers.Adam(), loss=CustomMSE())
class CustomMSE(tf.keras.losses.Loss):
    def __init__(self, regularization_factor=0.1, name="custom_mse"):
        super().__init__(name=name)
        self.regularization_factor = regularization_factor

    def call(self, y_true, y_pred):
        mse = tf.math.reduce_mean(tf.square(y_true - y_pred))
        reg = tf.math.reduce_mean(tf.square(0.5 - y_pred))
        return mse + reg * self.regularization_factor


# //--------------------------------------------
# //--------------------------------------------


#See: https://www.tensorflow.org/guide/keras/train_and_evaluate
# Custom metric that counts how many samples were correctly classified as belonging to a given class
# Example: model.compile(..., metrics=[CategoricalTruePositives()])
class CategoricalTruePositives(tf.keras.metrics.Metric):
    def __init__(self, name="categorical_true_positives", **kwargs):
        super(CategoricalTruePositives, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name="ctp", initializer="zeros")

    def update_state(self, y_true, y_pred, sample_weight=None):
        y_pred = tf.reshape(tf.argmax(y_pred, axis=1), shape=(-1, 1))
        values = tf.cast(y_true, "int32") == tf.cast(y_pred, "int32")
        values = tf.cast(values, "float32")
        if sample_weight is not None:
            sample_weight = tf.cast(sample_weight, "float32")
            values = tf.multiply(values, sample_weight)
        self.true_positives.assign_add(tf.reduce_sum(values))

    def result(self):
        return self.true_positives

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)


# //--------------------------------------------
# //--------------------------------------------

#Define here the loss function / optimizer / metrics to be used to train the model
def Get_Loss_Optim_Metrics(opts):

    _lr = opts["learnRate"]

    #-- These values are hard-coded here
    _momentum = 0.5 #helps preventing oscillations. Usually 0.5 - 0.9
    _decay = 0.0 #Decreases the _lr by specified amount after each epoch. Used in similar way as LearningRateScheduler
    _nesterov = False #improved momentum

    #-- Some possible choices of optimizers
    if opts["optimizer"] == "Adam": optim = Adam(lr=_lr, decay=_decay) #default lr=0.001
    elif opts["optimizer"] == "Nadam": optim = Nadam(learning_rate=_lr, beta_1=0.9, beta_2=0.999, decay=_decay) #default lr=0.001
    elif opts["optimizer"] == "Adadelta": optim = Adadelta(learning_rate=_lr, rho=0.95, epsilon=1e-07) #default lr=0.001
    elif opts["optimizer"] == "AdaBound": optim = AdaBound(learning_rate=_lr, final_learning_rate=0.1, gamma=1e-03, weight_decay=0., amsbound=False)
    elif opts["optimizer"] == "RMSprop": optim = RMSprop(lr=_lr)
    elif opts["optimizer"] == "SGD": optim = SGD(lr=_lr, decay=_decay, momentum=_momentum, nesterov=_nesterov)
    else: print(colors.fg.red, "ERROR: unknown optimizer algorithm ", opts["optimizer"], colors.reset); exit(1)

    #LOSS -- function used to optimize the model (minimized). Must be differentiable (for gradient method)
    # Examples : binary_crossentropy, categorical_crossentropy, mean_squared_error, , ...

    #METRICS -- used to judge performance of model. Not related to training. Can be used (e.g. with Keras' callbacks) to assess model's performance at given stages, or to stop training at some point (EarlyStopping, etc.)
    # metrics = 'accuracy'
    # metrics = 'binary_accuracy' #Calculates the mean accuracy rate across all predictions for binary classification problems.
    # metrics = 'categorical_accuracy'#Calculates the mean accuracy rate across all predictions for multiclass classification problems.
    # metrics = 'mean_squared_error' #Calculates the mean squared error (mse) rate between predicted and target values.
    # metrics = 'mean_squared_logarithmic_error'
    # metrics = 'mean_absolute_error' #Calculates the mean absolute error (mae) rate between predicted and target values.
    # metrics = 'hinge' #Calculates the hinge loss, which is defined as max(1 - y_true * y_pred, 0).
    # metrics = 'binary_crossentropy' #Calculates the cross-entropy value for binary classification problems.
    # metrics = 'AUC' #tf.metrics.AUC #Computes the approximate AUC (Area under the curve) via a Riemann sum.

    if opts["regress"]==False: #Classification
        if opts["nofOutputNodes"] == 1: #BINARY
            loss = 'binary_crossentropy'
            metrics = 'binary_accuracy'
            # metrics = 'AUC'

        elif opts["nofOutputNodes"] > 1: #MULTI
            loss = 'categorical_crossentropy'
            metrics = 'categorical_accuracy'
            # metrics = 'AUC'

    else: #Regression
        loss = 'mean_squared_error'
        # loss = 'mean_squared_logarithmic_error'
        # loss = 'mean_absolute_error' #More robust to outliers

        # if opts["strategy"] is "ROLR": loss = clipped_mse #use custom (clipped MSE) loss to avoid huge loss values dominating the training

        metrics = loss

    lossWeights = None
    if opts["strategy"] == "RASCAL": #2 outputs : [r,t]
        loss = ['mean_squared_error', 'mean_squared_error']
        # loss = ['mean_squared_logarithmic_error', 'mean_squared_logarithmic_error']
        lossWeights = [1, opts["score_lossWeight"]] # Apply scale factor to score loss weight

    # loss = BinaryFocalLoss(gamma=2) #Larger gamma <-> more focus on hard-to-classify events

    return loss, optim, metrics, lossWeights

# //--------------------------------------------
# //--------------------------------------------

# def clipped_mse(y_true, y_pred, min=-10., max=10.):
def clipped_mse(y_true, y_pred):

    return K.mean(K.square(K.clip(y_pred, -10., 10.) - K.clip(y_true, -10., 10.)), axis=-1)

# //--------------------------------------------
# //--------------------------------------------
