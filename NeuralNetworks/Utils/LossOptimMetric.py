#Define the loss function, optimizer and metric to be used to train the network

from tensorflow.keras.optimizers import SGD, Adam, RMSprop, Adadelta
from tensorflow.keras import backend as K
from Utils.ColoredPrintout import colors
# from keras_adabound import AdaBound
from adabound import AdaBound

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

#Define here the loss function / optimizer / metrics to be used to train the model
def Get_Loss_Optim_Metrics(opts):

    _lr = opts["learnRate"]

    #-- These values are hard-coded here
    _momentum = 0.5 #helps preventing oscillations. Usually 0.5 - 0.9
    _decay = 0.0 #Decreases the _lr by specified amount after each epoch. Used in similar way as LearningRateScheduler
    _nesterov = False #improved momentum

    #-- Some possible choices of optimizers
    if opts["optimizer"] == "Adam": optim = Adam(lr=_lr, decay=_decay) #default lr=0.001
    elif opts["optimizer"] == "Adadelta": optim = Adadelta(learning_rate=_lr, rho=0.95, epsilon=1e-07) #default lr=0.001
    elif opts["optimizer"] == "AdaBound": optim = AdaBound(lr=1e-03, final_lr=0.1, gamma=1e-03, weight_decay=0., amsbound=False)
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
        if opts["nofOutputNodes"] > 1:
            loss = 'categorical_crossentropy'
            metrics = 'categorical_accuracy'
            # metrics = 'AUC'
        elif opts["nofOutputNodes"] == 1:
            loss = 'binary_crossentropy'
            metrics = 'binary_accuracy'
            # metrics = 'AUC'

    else: #Regression
        loss = 'mean_squared_error'
        # loss = 'mean_squared_logarithmic_error'
        # loss = 'mean_absolute_error' #More robust to outliers
        # loss = 'mean_squared_logarithmic_error'

        # if opts["strategy"] is "ROLR": loss = clipped_mse #use custom (clipped MSE) loss to avoid huge loss values dominating the training

        # metrics = 'mean_absolute_error' #More robust to outliers
        metrics = loss

    lossWeights = None
    if opts["strategy"] == "RASCAL": #2 outputs : [r,t]
        loss = ['mean_squared_error', 'mean_squared_error']
        # loss = ['mean_squared_logarithmic_error', 'mean_squared_logarithmic_error']
        lossWeights = [1, opts["score_lossWeight"]] # Apply scale factor to score loss weight

    return loss, optim, metrics, lossWeights

# //--------------------------------------------
# //--------------------------------------------

# def clipped_mse(y_true, y_pred, min=-10., max=10.):
def clipped_mse(y_true, y_pred):

    return K.mean(K.square(K.clip(y_pred, -10., 10.) - K.clip(y_true, -10., 10.)), axis=-1)

# //--------------------------------------------
# //--------------------------------------------
