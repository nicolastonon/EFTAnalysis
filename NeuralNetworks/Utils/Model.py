# Define the NN architecture and training settings

import numpy as np
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Lambda, Input, Dense, Dropout, AlphaDropout, Activation, BatchNormalization, LeakyReLU, PReLU
from tensorflow.keras.activations import relu
from tensorflow.keras import regularizers
from tensorflow.keras.regularizers import l2
from keras.initializers import Constant
from tensorflow.keras.callbacks import TensorBoard, EarlyStopping, LambdaCallback, LearningRateScheduler, ReduceLROnPlateau
from tensorflow.keras import backend as K
from Utils.Helper import normalize

# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
 ######  ########  ########    ###    ######## ########    ##     ##  #######  ########  ######## ##
##    ## ##     ## ##         ## ##      ##    ##          ###   ### ##     ## ##     ## ##       ##
##       ##     ## ##        ##   ##     ##    ##          #### #### ##     ## ##     ## ##       ##
##       ########  ######   ##     ##    ##    ######      ## ### ## ##     ## ##     ## ######   ##
##       ##   ##   ##       #########    ##    ##          ##     ## ##     ## ##     ## ##       ##
##    ## ##    ##  ##       ##     ##    ##    ##          ##     ## ##     ## ##     ## ##       ##
 ######  ##     ## ######## ##     ##    ##    ########    ##     ##  #######  ########  ######## ########
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Define here the NN model
def Create_Model(opts, outdir, list_features, shifts, scales, NN_name="NN"):

# //--------------------------------------------
# //--------------------------------------------
# Architecture options (set by user in main code)

    nHiddenLayers = opts["nHiddenLayers"]
    nNeuronsAllHiddenLayers = opts["nNeuronsAllHiddenLayers"]
    activInputLayer = opts["activInputLayer"]
    activHiddenLayers = opts["activHiddenLayers"]
    use_normInputLayer = opts["use_normInputLayer"]
    use_batchNorm = opts["use_batchNorm"]
    dropoutRate = opts["dropoutRate"]
    regress_onLogr = opts["regress_onLogr"]

    use_dropout = False
    if dropoutRate > 0: use_dropout = True

    nNeuronsPerHiddenLayer = []
    if "nNeuronsPerHiddenLayer" in opts: nNeuronsPerHiddenLayer = opts["nNeuronsPerHiddenLayer"]

    nof_outputs = opts["nofOutputNodes"]

# //--------------------------------------------
# //--------------------------------------------
# Automatically set some variables

    if activHiddenLayers in ['lrelu','prelu']: activHiddenLayers = None #These must be included as a layer, not an activation

    num_input_variables = len(list_features) #Nof input variables to be read by the NN
    # print(num_input_variables)

    # if activHiddenLayers is "lrelu": #Not recognized by some later code when freezing/saving model... ?
        # activHiddenLayers = LeakyReLU(alpha=0.01)
        # activHiddenLayers = lambda x: relu(x, alpha=0.1)

    kernInit = 'he_normal' #Hard-coded here
    # kernInit = 'glorot_normal'
    # kernInit = 'glorot_uniform'
    # kernInit = 'lecun_normal'

    #Regularizers
    #NB : "In practice, if you are not concerned with explicit feature selection, L2 regularization can be expected to give superior performance over L1."
    # reg = regularizers.l1(0.001) #Default 0.001
    # reg = regularizers.l2(0.0001) #Default 0.001
    # reg = regularizers.l1_l2(l1=0.01, l2=0.01)

    reg = None
    if opts["regularizer"][0] is 'l1': reg = regularizers.l1(opts["regularizer"][1])
    elif opts["regularizer"][1] is 'l2': reg = regularizers.l2(opts["regularizer"][1])
    elif opts["regularizer"][0] is 'l1l2': reg = regularizers.l1_l2(opts["regularizer"][1])

    #Examples of advanced activations (should be added as layers, after dense layers)
    # model.add(LeakyReLU(alpha=0.1))
    # model.add(PReLU(alpha_initializer=kernInit))
    # model.add(Activation('selu'))

# //--------------------------------------------
# //--------------------------------------------

# //--------------------------------------------
# INPUT LAYER

    inp = Input(shape=num_input_variables, name="MYINPUT") #(Inactive) input layer
    X = inp #First building block of model

    if use_normInputLayer == True :
        X = Lambda(normalize, arguments={'shift': shifts, 'scale': scales}, name="Feature_normalization",dtype='float32')(X)

# //--------------------------------------------
# HIDDEN LAYERS

    for iLayer in range(nHiddenLayers):

        # Can choose to use different activation for first layer reading input features
        activ_tmp = activHiddenLayers
        if iLayer == 0 and activInputLayer != '': activ_tmp = activInputLayer

        # Number of neurons for current hidden layer
        nNeurons = nNeuronsAllHiddenLayers
        if len(nNeuronsPerHiddenLayer) > 0: nNeurons = nNeuronsPerHiddenLayer[iLayer]

        # Dense hidden layers
        X = Dense(nNeurons, activation=activ_tmp, activity_regularizer=reg, kernel_initializer=kernInit)(X)
        if activ_tmp == 'lrelu': X = LeakyReLU(alpha=0.1)(X) #Arbitrary alpha parameter
        elif activ_tmp == 'prelu': X = PReLU(alpha_initializer=Constant(value=0.25))(X) #Arbitrary alpha init. as proposed by He et al. (2015)

        # Dropout
        if use_dropout==True and iLayer < nHiddenLayers-1: #Don't apply dropout after last hidden layer #CHANGED -- apply dropout before BN (else, still get some info from blind nodes passed through BN)
            X = Dropout(dropoutRate)(X)

        # Batch normalization
        if use_batchNorm==True and iLayer < nHiddenLayers-1: #CHANGED -- never add batchNorm before output layer ? From this question (to be verified!): https://stats.stackexchange.com/questions/361700/lack-of-batch-normalization-before-last-fully-connected-layer
            X = BatchNormalization()(X)
# //--------------------------------------------
# OUTPUT LAYER

    activOutput = "sigmoid" #Default activation

    if opts["regress"] == False: #Classification (sum(P)=1)

        if nof_outputs == 1: activOutput = "sigmoid" #Binary
        else: activOutput = "softmax" #Multiclass

        out = Dense(nof_outputs, kernel_initializer=kernInit, activation=activOutput, name="MYOUTPUT")(X)
        model = Model(inputs=[inp], outputs=[out])

    else: #Regression

        if opts["strategy"] in ["ROLR", "RASCAL"]:

            if regress_onLogr: #Apply linear activation and exponentiate the result to match the target <-> train on log(r)
                logr = Dense(1, activation="linear", name="linear")(X)
                r = Lambda(lambda v: K.exp(v), name="likelihood_ratio")(logr)
            else: #Apply linear activation to match the target <-> train on r. Also define log(r), may be needed to compute score
                r = Dense(1, activation="linear", name="likelihood_ratio")(X)
                logr = Lambda(lambda v: K.log(v), name="log")(r)

            if opts["strategy"] == "RASCAL": #Ratio+Score. Differentiate 'logr' layer to get the score
                t = Lambda(lambda_layer_score, arguments={"theta_dim": len(opts["listOperatorsParam"])}, name="score")([logr, inp])
                model = Model(inputs=[inp], outputs=[r, t]) #1 output for r, n_operator outputs for t
            else: #Ratio only
                # rclip = K.clip(r, min_value=-10, max_value=10) #NB: could clip gradient as such to avoid huge values (but also truncates prediction values... !)
                model = Model(inputs=[inp], outputs=[r])

        elif opts["strategy"] is "regressor":
            out = Dense(nof_outputs, activation="linear", name="linear")(X)
            model = Model(inputs=[inp], outputs=[out])

        else: print("ERROR: no regressor model defined for this strategy"); exit(1)

# //--------------------------------------------
# //--------------------------------------------

    #Model visualization
    print(model.summary())

    return model

# //--------------------------------------------
# //--------------------------------------------

# Returns the gradients of loss w.r.t. variables
# Here: v[0] is output of logr layer, v[1] are inputs
def lambda_layer_score(v, theta_dim):

    # grad = K.gradients(loss=v[0], variables=v[1][:-theta_dim])[0]
    grad = K.gradients(loss=v[0], variables=v[1])[0]

    if grad is None:
        grad = K.zeros_like(theta_dim)

    return grad[:, -theta_dim:] #Only return gradients w.r.t. theory parameters (**must be last input features !!**) --> Nof score outputs = nof operators

# //--------------------------------------------
# //--------------------------------------------
