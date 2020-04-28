# Nicolas TONON (DESY)
# Train fully-connected neural networks with Keras (tensorflow back-end) for classification and regression
# //--------------------------------------------

# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
 #######  ########  ######## ####  #######  ##    ##  ######
##     ## ##     ##    ##     ##  ##     ## ###   ## ##    ##
##     ## ##     ##    ##     ##  ##     ## ####  ## ##
##     ## ########     ##     ##  ##     ## ## ## ##  ######
##     ## ##           ##     ##  ##     ## ##  ####       ##
##     ## ##           ##     ##  ##     ## ##   ### ##    ##
 #######  ##           ##    ####  #######  ##    ##  ######
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

# Analysis options
# //--------------------------------------------
# -- Choose here what data you want to consider (separate ntuples per year) ; same convention as for main analysis code
# Naming convention enforced : 2016+2017 <-> "201617" ; etc.; 2016+2017+2018 <-> "Run2" # NB : years must be placed in the right order !
_list_lumiYears = []
# _list_lumiYears.append("2016")
_list_lumiYears.append("2017")
# _list_lumiYears.append("2018")

#Signal process must be first
_list_processClasses = []
# _list_processClasses.append(["tZq"])
_list_processClasses.append(["PrivMC_tZq"])
# _list_processClasses.append(["PrivMC_tZq_ctz"])
# _list_processClasses.append(["PrivMC_tZq_ctw"])
# _list_processClasses.append(["ttZ"])
# _list_processClasses.append(["ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep"])
# _list_processClasses.append(["ttZ", "ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep",])

_list_labels = []
# _list_labels.append("tZq")
_list_labels.append("PrivMC_tZq")
# _list_labels.append("PrivMC_tZq_ctz")
# _list_labels.append("PrivMC_tZq_ctw")
# _list_labels.append("ttZ")
# _list_labels.append("Backgrounds")
# //--------------------------------------------

#--- Training options
# //--------------------------------------------
optsTrain = {

#Classification / regression
"regress": False, #True <-> DNN used for regression ; False <-> classification
"target": "class", #'None', 'class'

#EFT
"parameterizedDNN": True, #True <-> DNN is parameterized on the Wilson coefficients of the EFT operators
"listOperatorsParam": ['ctZ','ctW'], #None <-> parameterize on all possible operators
# "listOperatorsParam": ['ctZ','ctW', 'cpQM', 'cpQ3', 'cpt'], #None <-> parameterize on all possible operators
"nPointsPerOperator": 20, "minWC": -5, "maxWC": 5, #Interval [min,max,step] in which EFT points get sampled uniformly to train the DNN on
"nEventsPerPoint": 10000, #max nof events to be used for each EFT point (for parameterized DNN only) ; -1 <-> use all available events

#Hyperparameters
"nepochs": 20, #Number of training epochs (<-> nof times the full training dataset is shown to the NN)
"batchSize": 5000, #Batch size (<-> nof events fed to the network before its parameter get updated)

# "maxEvents": 100, #max nof events to be used for each process class, *or for each EFT point* (in parameterized DNN) ; -1 <-> use all available events
"maxEventsPerClass": 100000, #max nof events to be used for each process class (non-parameterized DNN only) ; -1 <-> use all available events
"nEventsTot_train": -1, "nEventsTot_test": -1, #nof events to be used for training & testing ; -1 <-> use _maxEvents & _splitTrainEventFrac params instead
"splitTrainEventFrac": 0.8, #Fraction of events to be used for training (1 <-> use all requested events for training)

#Event preselection
"cuts": "1", #Event selection, both for train/test ; "1" <-> no cut
}

# _startFromExistingModel = False #True <-> Skip training, load latest checkpoint model and create perf plots #not used yet
# //--------------------------------------------

# Define list of input variables
# //--------------------------------------------
_list_features = []
_list_features.append("maxDijetDelR")
_list_features.append("dEtaFwdJetBJet")
_list_features.append("dEtaFwdJetClosestLep")
_list_features.append("mHT")
_list_features.append("mTW")
_list_features.append("Mass_3l")
_list_features.append("forwardJetAbsEta")
_list_features.append("jPrimeAbsEta")
_list_features.append("maxDeepCSV")
_list_features.append("delRljPrime")
_list_features.append("lAsymmetry")
_list_features.append("maxDijetMass")
_list_features.append("maxDelPhiLL")











# //--------------------------------------------
#Filtering out manually some unimportant warnings
# import warnings
# warnings.filterwarnings("ignore", message="tensorflow:sample_weight modes were coerced")

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import argparse
# //--------------------------------------------
import tensorflow
import keras
import numpy as np
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score
# from tensorflow.keras.models import load_model

from Utils.FreezeSession import FreezeSession_and_SaveModel #freeze_session
from Utils.Helper import *
from Utils.Model import Create_Model
from Utils.Callbacks import Get_Callbacks
from Utils.GetData import Get_Data
from Utils.Optimizer import Get_Loss_Optim_Metrics
from Utils.ColoredPrintout import colors
from Utils.Validation_Control import *
from Utils.Predictions import *
from Utils.DataGenerator import *
# //--------------------------------------------
# //--------------------------------------------





# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
######## ########     ###    #### ##    ##
   ##    ##     ##   ## ##    ##  ###   ##
   ##    ##     ##  ##   ##   ##  ####  ##
   ##    ########  ##     ##  ##  ## ## ##
   ##    ##   ##   #########  ##  ##  ####
   ##    ##    ##  ##     ##  ##  ##   ###
   ##    ##     ## ##     ## #### ##    ##

######## ########  ######  ########
   ##    ##       ##    ##    ##
   ##    ##       ##          ##
   ##    ######    ######     ##
   ##    ##             ##    ##
   ##    ##       ##    ##    ##
   ##    ########  ######     ##

######## ##     ##    ###    ##
##       ##     ##   ## ##   ##
##       ##     ##  ##   ##  ##
######   ##     ## ##     ## ##
##        ##   ##  ######### ##
##         ## ##   ##     ## ##
########    ###    ##     ## ########
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Main function, calling sub-functions to perform all necessary actions
def Train_Test_Eval_DNN(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features):

 # #    # # #####
 # ##   # #   #
 # # #  # #   #
 # #  # # #   #
 # #   ## #   #
 # #    # #   #

    #-- Initialization, sanity checks
    optsTrain, _lumiName, _weightDir, _ntuplesDir, _h5modelName = Initialization_And_SanityChecks(optsTrain, _list_lumiYears, _list_processClasses, _list_labels)


                                       #
 ##### #####    ##   # #    #         #     ##### ######  ####  #####
   #   #    #  #  #  # ##   #        #        #   #      #        #
   #   #    # #    # # # #  #       #         #   #####   ####    #
   #   #####  ###### # #  # #      #          #   #           #   #
   #   #   #  #    # # #   ##     #           #   #      #    #   #
   #   #    # #    # # #    #    #            #   ######  ####    #

    print('\n\n')
    print(colors.bg.orange, colors.bold, "=====================================", colors.reset)
    print('\t', colors.fg.orange, colors.bold, "DNN Training", colors.reset)
    print(colors.bg.orange, colors.bold, "=====================================", colors.reset, '\n\n')

    #-- Get data
    print(colors.fg.lightblue, "--- Read and shape the data...", colors.reset); print('\n')
    x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, shifts, scales, x_control_firstNEvents, xTrainRescaled = Get_Data(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features, _weightDir, _ntuplesDir, _lumiName)

    # From there, for parameterized DNN, the different 'classes' correspond to different operators, and must include their WCs as inputs
    if optsTrain["parameterizedDNN"]==True :
        _list_labels = optsTrain["listOperatorsParam"][:]; _list_labels.insert(0, "SM") #Specify '[:]' to create a copy, not a reference
        _list_features = np.append(_list_features, optsTrain["listOperatorsParam"])

    #-- Plot input features distributions, after applying to train data same rescaling as will be done by first DNN layer (-> check rescaling)
    Plot_Input_Features(optsTrain, xTrainRescaled, y_process_train, PhysicalWeights_train, _list_features, _weightDir, True)

    print('\n'); print(colors.fg.lightblue, "--- Define the loss function & metrics...", colors.reset); print('\n')
    _loss, _optim, _metrics = Get_Loss_Optim_Metrics(optsTrain['regress'], optsTrain['nofOutputNodes'])

    #-- Get model and compile it
    print('\n'); print(colors.fg.lightblue, "--- Create the Keras model...", colors.reset); print('\n')
    model = Create_Model(optsTrain, _weightDir, _list_features, shifts, scales)

    # Can printout the output of the ith layer here for N events, e.g. to verify that the normalization layer works properly
    # Printout_Outputs_FirstLayer(model, ilayer=0, xx=x[0:5])

    print('\n'); print(colors.fg.lightblue, "--- Compile the Keras model...", colors.reset); print('\n')
    model.compile(loss=_loss, optimizer=_optim, metrics=[_metrics]) #For multiclass classification

    #-- Define list of callbacks
    callbacks_list = Get_Callbacks(_weightDir)
    # ckpt_dir = os.path.dirname(ckpt_path); history = 0

    #-- Fit model (TRAIN)
    print('\n'); print(colors.fg.lightblue, "--- Train (fit) DNN on training sample...", colors.reset, " (may take a while)"); print('\n')
    if optsTrain["parameterizedDNN"]==False:
        history = model.fit(x_train, y_train, sample_weight=LearningWeights_train, validation_data=(x_test, y_test, PhysicalWeights_test), epochs=optsTrain["nepochs"], batch_size=optsTrain["batchSize"], callbacks=callbacks_list, shuffle=True, verbose=1)

        # Evaluate the neural network's performance (evaluate metrics on validation or test dataset)
        print('\n'); print(colors.fg.lightblue, "--- Evaluate DNN performance on test sample...", colors.reset); print('\n')
        score = model.evaluate(x_test, y_test, batch_size=_batchSize, sample_weight=PhysicalWeights_test, verbose=1)

    else:
        my_training_batch_generator = DataGenerator(x_train, y_train, LearningWeights_train, optsTrain["batchSize"])
        my_validation_batch_generator = DataGenerator(x_test, y_test, PhysicalWeights_test, optsTrain["batchSize"])
        _steps_per_epoch = np.ceil(len(x_train) / optsTrain["batchSize"]); _steps_per_epoch_val = np.ceil(len(x_test)/ optsTrain["batchSize"])
        history = model.fit(my_training_batch_generator, steps_per_epoch=_steps_per_epoch, validation_data=my_validation_batch_generator, validation_steps=_steps_per_epoch_val, epochs=optsTrain["nepochs"], callbacks=callbacks_list, verbose=1)

        print('\n'); print(colors.fg.lightblue, "--- Evaluate DNN performance on test sample...", colors.reset); print('\n')
        score = model.evaluate(my_validation_batch_generator, steps=_steps_per_epoch_val, verbose=1)

# //--------------------------------------------
    #-- Can access weights and biases of any layer
    # weights_layer, biases_layer = model.layers[0].get_weights(); print(weights_layer.shape); print(biases_layer.shape); print(weights_layer); print(biases_layer[0:2])
    #-- Loads the latest checkpoint weights
    # latest = tensorflow.train.latest_checkpoint(ckpt_dir)
    # tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    # model = load_model(_h5modelName) # model has to be re-loaded
    # model.load_weights(latest)
# //--------------------------------------------


  ####    ##   #    # ######    #    #  ####  #####  ###### #
 #       #  #  #    # #         ##  ## #    # #    # #      #
  ####  #    # #    # #####     # ## # #    # #    # #####  #
      # ###### #    # #         #    # #    # #    # #      #
 #    # #    #  #  #  #         #    # #    # #    # #      #
  ####  #    #   ##   ######    #    #  ####  #####  ###### ######

    print('\n'); print(colors.fg.lightblue, "--- Save model...", colors.reset);

    # Serialize model to HDF5
    model.save(_h5modelName)

    # Save the model architecture
    with open(_weightDir + 'arch_DNN.json', 'w') as json_file:
        json_file.write(model.to_json())

    # Convert model to estimator and save model as frozen graph for c++
    with tensorflow.compat.v1.Session() as sess: FreezeSession_and_SaveModel(optsTrain, sess, _weightDir, _h5modelName) #Must first open a new session #Can't manage to run code below without this... (why?)


 #    #   ##   #      # #####    ##   ##### #  ####  #    #
 #    #  #  #  #      # #    #  #  #    #   # #    # ##   #
 #    # #    # #      # #    # #    #   #   # #    # # #  #
 #    # ###### #      # #    # ######   #   # #    # #  # #
  #  #  #    # #      # #    # #    #   #   # #    # #   ##
   ##   #    # ###### # #####  #    #   #   #  ####  #    #

    print('\n\n')
    print(colors.bg.orange, colors.bold, "##############################################", colors.reset)
    print(colors.fg.orange, '\t Results & Control Plots', colors.reset)
    print(colors.bg.orange, colors.bold, "##############################################", colors.reset, '\n')

    #-- Get control results (printouts, plots, histos)
    print('\n', colors.fg.lightblue, "--- Apply model to train & test data...", colors.reset, " (may take a while)\n")
    list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses = Apply_Model_toTrainTestData(optsTrain, _list_labels, x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, _h5modelName, x_control_firstNEvents)

    Control_Printouts(optsTrain["nofOutputNodes"], score, _list_labels, y_test, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses)

    Make_TrainTestPrediction_Histograms(optsTrain["nofOutputNodes"], _lumiName, _list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, _metrics)

    Create_Control_Plots(optsTrain, _list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, x_train, x_test, y_train, y_test, y_process_train, y_process_test, _h5modelName, _metrics, _weightDir, history)

    Create_Correlation_Plot(x, _list_features, _weightDir)

    Plot_Input_Features(optsTrain, x, y_process, PhysicalWeights_allClasses, _list_features, _weightDir, False)

# //--------------------------------------------
# //--------------------------------------------
























# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##
##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ##
##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ##
######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##
##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####
##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ###
##        #######  ##    ##  ######     ##    ####  #######  ##    ##


 ######     ###    ##       ##        ######
##    ##   ## ##   ##       ##       ##    ##
##        ##   ##  ##       ##       ##
##       ##     ## ##       ##        ######
##       ######### ##       ##             ##
##    ## ##     ## ##       ##       ##    ##
 ######  ##     ## ######## ########  ######
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#----------  Manual call to DNN training function
Train_Test_Eval_DNN(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features)

# //--------------------------------------------
#-- Set up the command line arguments
# if __name__ == "__main__":
# parser = argparse.ArgumentParser()
# parser.add_argument("xxx", metavar="xxx", help="help")
