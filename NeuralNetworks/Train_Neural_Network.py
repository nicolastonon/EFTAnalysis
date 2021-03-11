# Nicolas TONON (DESY)
# Train fully-connected neural networks with Keras (tensorflow back-end) for classification and regression, with advanced functionnalities for EFT inference
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

# Training options
# //--------------------------------------------
optsTrain = {

#=== NTuple settings ===#
"TTree": 'result', # Name of the root TTree containing input features
"eventWeightName": '', #'' <-> hardcoded var name for my own NTuples; otherwise, use the specified var for per-event weights

#=== NN strategy ===#
"parameterizedNN": False, #True <-> include WCs of selected EFT operators as additional input features (--> parameterize the NN on the WC values) #Only valid for NN strategies requiring to train over mixture of SMEFT hypotheses (e.g. CARL, ...) #NB: much simpler, but sensible for multi-dim parameter space... ?

# "strategy": "classifier", # <-> Regular classifier: separates events from different samples [central or pure-EFT samples only]
# "strategy": "regressor", # <-> Regular regressor: regress some quantity for different samples #CHOOSE MODE IN Get_Data.py !
# "strategy": "CARL_singlePoint", # <-> Calibrated Classifier: separates SM from single EFT point [EFT samples only]
"strategy": "CARL", # <-> Calibrated Classifier: separates points in EFT phase space via classification, single output node [EFT samples only, parameterized]
# "strategy": "CARL_multiclass", # <-> Calibrated Classifier: separates points in EFT phase space via classification, 1 output node per EFT operator [EFT samples only, parameterized] #NB: useless for multi-dim fits, since can't 'classify' EFT operators at mixed point; only advantage would be for 1D limits (replace N training by 1, dealing with all N operators)... #NB: may nto work if 'parameterizedNN==False'
# "strategy": "ROLR", # <-> Ratio Regression: regresses likelihood ratio between ref point and any EFT point [EFT samples only, parameterized]
# "strategy": "RASCAL", # <-> Ratio+Score Regression: same as ROLR, but also include score info in training [EFT samples only, parameterized]
# "strategy": "CASCAL", # <-> Calibrated Classifier (CARL) + score regression; still in test phase

#=== General training/architecture settings ===#
"nEpochs": 50, #Number of training epochs (<-> nof times the full training dataset is shown to the NN)
"splitTrainValTestData": [0.80, 0.0, 0.20], #Fractions of events to be used for the training / validation (evaluation after each epoch) / test (final evaluation) datasets respectively #If frac_val=0, only split between train/test data (not ideal but may be necessary if stat. is too low) #Superseeded by 'nEventsTot_train/val/test' (when trainAtManyEFTpoints==False)
# "splitTrainEventFrac": 0.80, #Fraction of events to be used for training (1 <-> use all requested events for training)
"nHiddenLayers": 3, #Number of hidden layers
"nNeuronsAllHiddenLayers": 50, #Number of neurons per same-size hidden layer
# "nNeuronsPerHiddenLayer": [128,64,32,16], #Number of neurons per same-size hidden layer
"activInputLayer": 'relu', #Activation function for 1st hidden layer (connected to input layer) # '' <-> use same as for activHiddenLayers #NB: don't use lrelu/prelu/... for first layer (neglect info. ?) !
"activHiddenLayers": 'relu', #Activation function for hidden layers #sigmoid,tanh,relu,lrelu,prelu,selu,...
"use_normInputLayer": True, #True <-> add a transformation layer to rescale input features
"use_batchNorm": True, #True <-> apply batch normalization after each hidden layer
"dropoutRate": 0.4, #Dropout rate (0 <-> disabled) #Use to avoid overtraining for complex architectures only, and with sufficient nof epochs
"regularizer": ['L2', 0.0001], #Weight regularization: '' (<-> None), 'L1','L2','L1L2' <-> apply value given in 2nd arg.
"optimizer": "Adam", #Optimization algorithm: 'SGD', 'RMSprop', 'Adam', 'Nadam','Adadelta','AdaBound',... #See basic explanations here: https://medium.com/@sdoshi579/optimizers-for-training-neural-network-59450d71caf6
"learnRate": 0.001, #Learning rate (initial value) of optimizer. Too low -> weights don't update. Too large -> Unstable, no convergence #Default (Adam): 0.001
"balancedClasses": True, #True <-> apply event SFs in training to balance the total weights of all classes
"earlyStopping": True, #True <-> use Keras' early stopping

#=== Settings for NN *not* trained at many different SMEFT points (parameterized or not) ===#
"maxEventsPerClass": -1, #Max. nof events to be used for each process class (non-parameterized NN only) ; -1 <-> use all available events #NB: setting a max nof events may bias the training ? Since it may correspond to ~full stat for some samples and to very partial stat for others <-> normalization yields get biased !
"nEventsTot_train": -1, "nEventsTot_val": -1, "nEventsTot_test": -1, #total nof events to be used for train / val / test; -1 <-> use _maxEvents & splitTrainValTestData params instead
"batchSizeClass": 1000, #Batch size (<-> nof events fed to the network before its parameter get updated)

#=== Settings for CARL/ROLR/RASCAL strategies ===#
"refPoint": "SM", #Reference point used e.g. to compute likelihood ratios. Must be "SM" for CARL_multiclass strategy (<-> separate SM from EFT). Must be != "SM" for CARL_singlePoint strategy (<-> will correspond to the single hypothesis to separate from SM). Follow naming convention from MG, e.g.: 'ctZ_-3.5_ctp_2.6'
# "refPoint": "rwgt_cpqm_15",
# "refPoint": "rwgt_ctw_2",
# "refPoint": "rwgt_ctw_5",
# "refPoint": "rwgt_ctz_5",
# "listOperatorsParam": ['ctz', 'ctw', 'cpqm', 'cpq3', 'cpt'], #None <-> parameterize on all possible operators
# "listOperatorsParam": ['ctz','ctw', 'cpq3'], #None <-> parameterize on all possible operators
# "listOperatorsParam": ['ctz', 'ctw'], #None <-> parameterize on all possible operators
"listOperatorsParam": ['ctw'], #None <-> parameterize on all possible operators
"nPointsPerOperator": 50, #Number of EFT points from which to sample events (times the number of selected operators)
"minWC": -5, "maxWC": 5, #Interval [min,max] in which EFT points get sampled uniformly to train the NN on
# "listMinMaxWC": [-5,5, -5,5, -15,15, -10,10, -15,15], #If activated, and len(listMinMaxWC)=2*len(listOperatorsParam), will be interpreted as a list of min/max values for each operator selected above for NN parameterization (superseeds minWC/maxWC values) #NB: must keep ranges for all operators (found in samples), even if considering only a subset of operators !
"nEventsPerPoint": 10000, #max nof events to be used for each EFT point (for parameterized NN only) ; -1 <-> use all available events
"batchSizeEFT": 1000, #Batch size (<-> nof events fed to the network before its parameter get updated)
"score_lossWeight": 0.5, #Apply scale factor to score term in loss function
"regress_onLogr": False, #True <-> NN will regress on log(r) instead of r

#=== Settings for regressor strategy ===#
"targetVarIdx": -1, #List of indices *in the list of input features* (NB: only for convenience) of variable(s) to use as target(s) for regression; the var(s) get removed from training and from the list later. If multiple indices provided, multiple are regressed. If set to < 0, the target will be defined in the Get_Targets() function
"comparVarIdx": -1, #Index *in the list of input features* of a var to compare to predictions in some validation plots (e.g.: Truth vs Pred vs kinReco). If < 0, not used

#=== Event preselection ===#
# "cuts": "1", #Event selection, both for train/test ; "1" <-> no cut
"cuts": "is_signal_SR",

#=== Input features ===#
"useHardCodedListInputFeatures": True, #True <-> use list of input features hard-coded in 'InputFeatures.py' (can define several for specific cases); otherwise, use the list of features defined here below
"useLowLevelFeatures": False, #True <-> include P4 vectors corresponding to 3 selected leptons, and up to 4 hardest jets #(+ btagging score) not for now

#=== OTHERS ===#
"makeValPlotsOnly": False, #True <-> load pre-existing model, skip train/test phase, create validation plots directly. Get data first (needed for plots)
"testToy1D": False, #True <-> Testing (expert) mode: try to replicate 1D toy example from arXiv:1601.07913, to debug/understand basic paramNN
"storeInTestDirectory": True, #True <-> all results (weights, plots, etc.) overwrite existing files in a common dir.; False <-> store results in specific sub-dir., depending on user-options, following path conventions of main analysis code
"storePerOperatorSeparately": True, #True <-> when training on SM vs EFT and considering a single operator, will store all outputs in an operatgor-specific dir. (allows to then consider different trainigs for different operators)
"useFakeableNPL": True, #True <-> if running on MC NPL samples (hardcoded names), will consider events from both SR+sideband (i.e. events with ==3FO,2 or 3 tights instead of ==3 tight leptons); for normalization, use that of ==3 tight lepton events only (otherwise need FR weights)
"displayImages": False, #True <-> will display (pop-up) some output figures (metrics, loss, ROC, ...)
"shapPlots": False, #True <-> create SHAP validation plots
}

# Analysis options
# //-------------------------------------------

# -- Choose the data to consider #NB: same convention as for main analysis code. Naming convention enforced : 2016+2017 <-> "201617" ; etc.; 2016+2017+2018 <-> "Run2" #NB: years must be placed in the right order !
_list_lumiYears = []
_list_lumiYears.append("2016")
_list_lumiYears.append("2017")
_list_lumiYears.append("2018")

#-- Choose the classes of processes to consider #NB: can group several physics processes in same process class #NB: place main signal in first position
_list_processClasses = []
# _list_processClasses.append(["tZq"])
# _list_processClasses.append(["ttZ"])
# _list_processClasses.append(["tZq", "ttZ"])
_list_processClasses.append(["PrivMC_tZq_noPS"])
# _list_processClasses.append(["PrivMC_ttZ_noPS"])
# _list_processClasses.append(["PrivMC_tZq"])
# _list_processClasses.append(["PrivMC_ttZ"])
# _list_processClasses.append(["PrivMC_tZq_TOP19001"])
# _list_processClasses.append(["PrivMC_ttZ_TOP19001"])
# _list_processClasses.append(["PrivMC_tZq_ctz"])
# _list_processClasses.append(["ttW", "ttH", "WZ", "ZZ4l"]) #Bkg
# _list_processClasses.append(["tX", "WZ", "VVV"]) #Bkg -- grouped
# _list_processClasses.append(["TTbar_DiLep", "DY"])
# _list_processClasses.append(["ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep"])
# _list_processClasses.append(["ttZ", "ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep"])
# _list_processClasses.append(["ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep", "DY"]) #WARNING: too low statistics for DY/ZGToLLG_01J, degrades performance

#-- Define labels associated with each process class #NB: keyword 'PrivMC' is used to denote private EFT samples
_list_labels = []
# _list_labels.append("tZq")
# _list_labels.append("ttZ")
_list_labels.append("PrivMC_tZq")
# _list_labels.append("PrivMC_ttZ")
# _list_labels.append("PrivMC_tZq_TOP19001")
# _list_labels.append("PrivMC_ttZ_TOP19001")
# _list_labels.append("PrivMC_tZq_ctz")
# _list_labels.append("PrivMC_ttZ_ctz")
# _list_labels.append("SM")
# _list_labels.append("Backgrounds")

# //--------------------------------------------

#-- Choose input features x
_list_features = []

_list_features.append("recoZ_Pt")
_list_features.append("recoZ_Eta")
# _list_features.append("recoZ_dPhill")
# _list_features.append("recoLepTop_Pt")

_list_features.append("dEta_tjprime")
_list_features.append("dR_bjprime")
_list_features.append("dR_jprimeClosestLep")
_list_features.append("dR_lWjprime")
_list_features.append("dR_tZ")
_list_features.append("lAsymmetry")

# '''
# _list_features.append("lAsymmetry")
# _list_features.append("jPrimeAbsEta")
# _list_features.append("mHT")
# _list_features.append("mTW")
# _list_features.append("Mass_3l")
# _list_features.append("maxEtaJet")
# _list_features.append("maxDeepJet")
# _list_features.append("njets")
# _list_features.append("nbjets")
# '''

'''
_list_features.append("cosThetaStarPolTop")
_list_features.append("cosThetaStarPolZ")
'''

'''
_list_features.append("recoLepTop_Pt") #!
_list_features.append("recoLepTop_Eta") #!
_list_features.append("TopZsystem_M") #!
_list_features.append("mbjMax")
_list_features.append("maxDiJet_Pt")
_list_features.append("maxDelRbL")
_list_features.append("minDelRbL")
_list_features.append("dR_ZlW")
_list_features.append("dR_blW")
_list_features.append("dR_tClosestJet")
_list_features.append("dR_bW") #!
_list_features.append("dEta_jprimeClosestLep")
'''

'''
_list_features.append("recoLepTopLep_Pt")
_list_features.append("dR_tjprime")
_list_features.append("dEta_bjprime") #!
_list_features.append("dR_bjprime") #!
_list_features.append("dR_lWjprime") #!
_list_features.append("dR_Zjprime") #!
_list_features.append("dEta_lWjprime") #!
'''




# //--------------------------------------------
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

from Utils.FreezeSave import FreezeSession_and_SaveModel #freeze_session
from Utils.Helper import *
from Utils.Model import Create_Model
from Utils.Callbacks import Get_Callbacks
from Utils.GetData import Get_Data
from Utils.LossOptimMetric import Get_Loss_Optim_Metrics
from Utils.ColoredPrintout import colors
from Utils.Validation_Control import *
from Utils.Predictions import *
from Utils.DataGenerator import *
import Utils.InputFeatures

# from trains import Task #Uncomment to use Trains to log/display some training parameters
# task = Task.init(project_name="my project", task_name="my task")

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
# 0 = all messages are logged (default behavior)
# 1 = INFO messages are not printed
# 2 = INFO and WARNING messages are not printed
# 3 = INFO, WARNING, and ERROR messages are not printed
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
def Train_Test_Eval_NN(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features):

    print('\n\n')
    print(colors.bg.orange, colors.bold, "=====================================", colors.reset)
    print('\t', colors.fg.orange, colors.bold, "NN Training", colors.reset)
    print(colors.bg.orange, colors.bold, "=====================================", colors.reset, '\n\n')


 # #    # # #####
 # ##   # #   #
 # # #  # #   #
 # #  # # #   #
 # #   ## #   #
 # #    # #   #

    #-- Initialization, sanity checks
    _lumiName, _weightDir, _h5modelName, _batchSize, _list_features = Initialization_And_SanityChecks(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features)
    print(colors.fg.lightgrey, '\n===> Saving NN settings to: ', _weightDir + "NN_settings.txt", colors.reset)
    print(colors.fg.lightgrey, '\n===> Saving NN features list and node names to: ', _weightDir + "NN_info.txt", colors.reset)

                                       #
 ##### #####    ##   # #    #         #     ##### ######  ####  #####
   #   #    #  #  #  # ##   #        #        #   #      #        #
   #   #    # #    # # # #  #       #         #   #####   ####    #
   #   #####  ###### # #  # #      #          #   #           #   #
   #   #   #  #    # # #   ##     #           #   #      #    #   #
   #   #    # #    # # #    #    #            #   ######  ####    #

    #-- Get data
    print(colors.fg.lightblue, "\n\n--- Get the data...\n", colors.reset)
    x_train, x_val, x_test, y_train, y_val, y_test, y_process_train, y_process_val, y_process_test, PhysicalWeights_train, PhysicalWeights_val, PhysicalWeights_test, LearningWeights_train, LearningWeights_val, LearningWeights_test, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, shifts, scales, xTrainRescaled, _list_labels, _list_features, jointLR_allClasses, scores_allClasses_eachOperator = Get_Data(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features, _weightDir, ntuplesDir, _lumiName)

    #-- Plot input features distributions, after applying to train data same rescaling as will be done by first NN layer (-> check rescaling)
    Plot_Input_Features(optsTrain, xTrainRescaled, y_process_train, PhysicalWeights_train, _list_features, _weightDir, True)

    print('\n'); print(colors.fg.lightblue, "--- Define the loss function & metrics...", colors.reset); print('\n')
    _loss, _optim, _metrics, _lossWeights = Get_Loss_Optim_Metrics(optsTrain)

    if optsTrain["makeValPlotsOnly"] == False: #Train, test, save, model
        #-- Get model and compile it
        print('\n'); print(colors.fg.lightblue, "--- Create the Keras model...", colors.reset); print('\n')
        model = Create_Model(optsTrain, _weightDir, _list_features, shifts, scales)

        print('\n'); print(colors.fg.lightblue, "--- Compile the Keras model...", colors.reset); print('\n')
        model.compile(loss=_loss, loss_weights=_lossWeights, optimizer=_optim, metrics=[_metrics]) #For multiclass classification

        #-- Define list of callbacks
        callbacks_list = Get_Callbacks(optsTrain, _weightDir)
        # ckpt_dir = os.path.dirname(ckpt_path); history = 0

        #-- Fit model (TRAIN)
        #NB: need to use weighted abs weights for val. data too... (as for training) !
        print('\n'); print(colors.fg.lightblue, "--- Train (fit) NN on training sample...", colors.reset); print('\n')
        if optsTrain["trainAtManyEFTpoints"]==False and len(x_train)<200000: #Use the default data generator (load entire samples at once)
            history = model.fit(x_train, y_train, sample_weight=LearningWeights_train, validation_data=(x_val, y_val, LearningWeights_val), epochs=optsTrain["nEpochs"], batch_size=_batchSize, callbacks=callbacks_list, shuffle=True, verbose=1) #Use learning weights (abs, rescaled) for training/validation

            #-- Evaluate the neural network's performance (evaluate metrics on validation or test dataset) #Use physical weights
            print('\n'); print(colors.fg.lightblue, "--- Evaluate NN performance on test sample...", colors.reset); print('\n')
            score = model.evaluate(x_test, y_test, batch_size=_batchSize, sample_weight=PhysicalWeights_test, verbose=1) #Use physical weights for testing

        else: #Use custom data generator, fit data in small batches, necessary for large samples
            my_training_batch_generator = DataGenerator(x_train, y_train, LearningWeights_train, optsTrain["strategy"], _batchSize, returnWeights=True) #Use learning weights (abs, rescaled) for training
            my_validation_batch_generator = DataGenerator(x_val, y_val, LearningWeights_val, optsTrain["strategy"], _batchSize, returnWeights=True) #Use learning weights (abs, rescaled) for validation
            my_testing_batch_generator = DataGenerator(x_test, y_test, PhysicalWeights_test, optsTrain["strategy"], _batchSize, returnWeights=True) #Use physical weights for testing

            #-- Cross check: compare mean x/y/w for train/val/test sets
            # for ifeat in range(len(_list_features)): print('Mean x feature', _list_features[ifeat], '-->', np.mean(x_train[:,ifeat]), '/', np.mean(x_val[:,ifeat]), '/', np.mean(x_test[:,ifeat]), '(train/val/test)')
            # print('Mean y -->', np.mean(y_train), '/', np.mean(y_val), '/', np.mean(y_test), '(train/val/test)')
            # print('Mean learning w -->', np.mean(LearningWeights_train), '/', np.mean(LearningWeights_val), '/', np.mean(LearningWeights_test), '(train/val/test)')
            # batch_x, batch_y, batch_weights = my_training_batch_generator.__getitem__(0); print(batch_x, batch_y)
            # batch_x, batch_y, batch_weights = my_validation_batch_generator.__getitem__(0); print(batch_x, batch_y)
            # exit(1)

            _steps_per_epoch_train = np.ceil(len(x_train) / _batchSize); _steps_per_epoch_val = np.ceil(len(x_val)/ _batchSize) ; _steps_per_epoch_test = np.ceil(len(x_test)/ _batchSize)
            # batch_x, batch_y = my_training_batch_generator.__getitem__(0); print(batch_x); print(batch_y)

            history = model.fit(my_training_batch_generator, steps_per_epoch=_steps_per_epoch_train, validation_data=my_validation_batch_generator, validation_steps=_steps_per_epoch_val, epochs=optsTrain["nEpochs"], callbacks=callbacks_list, verbose=1)

            print('\n'); print(colors.fg.lightblue, "--- Evaluate NN performance on test sample...", colors.reset); print('\n')
            score = model.evaluate(my_testing_batch_generator, steps=_steps_per_epoch_test, verbose=1)

        #-- Can printout the output of the i-th layer here for N events, e.g. to verify that the normalization layer works properly
        # for ilayer in range(len(model.layers)): Printout_Outputs_Layer(model, ilayer, xx=x[0:2])
        # Printout_Weights_Layer(model)


                                      #
  ####    ##   #    # ######         #     #       ####    ##   #####
 #       #  #  #    # #             #      #      #    #  #  #  #    #
  ####  #    # #    # #####        #       #      #    # #    # #    #
      # ###### #    # #           #        #      #    # ###### #    #
 #    # #    #  #  #  #          #         #      #    # #    # #    #
  ####  #    #   ##   ######    #          ######  ####  #    # #####

        print('\n'); print(colors.fg.lightblue, "--- Save model...", colors.reset);

        # Serialize full model (arch+weights+config+state of optimizer) to HDF5
        #Can then get this compiled model with 'model = load_model('xxx.h5')' (see: https://keras.io/getting-started/faq/#how-can-i-save-a-keras-model)
        model.save(_h5modelName)

        # Save the model architecture only to json format
        with open(_weightDir + 'arch_NN.json', 'w') as json_file:
            json_file.write(model.to_json())

        # Convert model to estimator and save model as frozen graph for c++
        with tensorflow.compat.v1.Session() as sess: FreezeSession_and_SaveModel(optsTrain, sess, _weightDir, _h5modelName, _list_labels) #Must first open a new session

    elif optsTrain["makeValPlotsOnly"] == True: #Load pre-existing model

        model = Load_PreExisting_Model(_h5modelName)
        score = None; history = None


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

    print(colors.fg.lightgrey, "\nTo open Tensorboard dir:",colors.reset,colors.ital,'tensorboard --logdir='+_weightDir+'logs'+' --port 0\n', colors.reset)

    #-- Get control results (printouts, plots, histos)
    list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, list_yTrain_allClasses, list_yTest_allClasses, list_xTrain_allClasses, list_xTest_allClasses = Apply_Model_toTrainTestData(optsTrain, _weightDir, _list_processClasses, _list_labels, x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, _h5modelName)

    #-- Store NN predictions for train/test datasets, for later use
    Store_TrainTestPrediction_Histograms(optsTrain, _lumiName, _list_features, _list_labels, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTest_allClasses, list_xTest_allClasses, list_predictions_train_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_xTrain_allClasses)

    #-- Create several validation plots automatically
    Make_Default_Validation_Plots(optsTrain, _list_features, _list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, PhysicalWeights_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, x, y_train, y_test, y_process, y_process_train, y_process_test, list_yTrain_allClasses, list_yTest_allClasses, list_xTrain_allClasses, list_xTest_allClasses, _metrics, _weightDir, scores_allClasses_eachOperator, model, score, history)

    NN_settings_filepath = "{}NN_settings.txt".format(_weightDir)
    Write_Timestamp_toLogfile(NN_settings_filepath, 1) #Write final timestamp before exit

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

#----------  Manual call to NN training function
if __name__ == "__main__":

    Train_Test_Eval_NN(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features)

# //--------------------------------------------
#-- Set up the command line arguments
# parser = argparse.ArgumentParser()
# parser.add_argument("xxx", metavar="xxx", help="help")
