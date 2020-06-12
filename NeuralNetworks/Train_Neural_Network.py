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

#=== NN strategy ===#
"strategy": "classifier", # <-> Regular classifier: separates events from different samples [central or pure-EFT samples only]
# "strategy": "regressor", # <-> Regular regressor: regress some quantity for different samples. Only label regression supported yet [central or pure-EFT samples only]
# "strategy": "CARL_singlePoint", # <-> Calibrated Classifier: separates SM from single EFT point [EFT samples only]
# "strategy": "CARL", # <-> Calibrated Classifier: separates points in EFT phase space via classification, single output node [EFT samples only, parameterized]
# "strategy": "CARL_multiclass", #BUGGED # <-> Calibrated Classifier: separates points in EFT phase space via classification, 1 output node per EFT operator [EFT samples only, parameterized]
# "strategy": "ROLR", # <-> Ratio Regression: regresses likelihood ratio between ref point and any EFT point [EFT samples only, parameterized]
# "strategy": "RASCAL", # <-> Ratio+Score Regression: same as ROLR, but also include score info in training [EFT samples only, parameterized]

#=== General training settings ===#
"nEpochs": 50, #Number of training epochs (<-> nof times the full training dataset is shown to the NN)
"splitTrainEventFrac": 0.8, #Fraction of events to be used for training (1 <-> use all requested events for training)

"nHiddenLayers": 4, #Number of hidden layers
"nNeuronsPerLayer": 100, #Number of neurons per hidden layer
"activInputLayer": 'tanh', #Activation function of input layer
"activHiddenLayers": 'relu', #Activation function of hidden layers #tanh,relu, ...
"use_normInputLayer": True, #True <-> add a transformation layer to rescale input features
"use_batchNorm": True, #True <-> apply batch normalization after each hidden layer
"dropoutRate": 0.4, #Dropout rate (0 <-> disabled) #Use to avoid overtraining for complex architectures only, and with sufficient nof epochs
"regularizer": ['L2', 0.0001], #Weight regularization: '' (<-> None), 'L1','L2','L1L2' <-> apply value given in 2nd arg.

#=== Settings for non-parameterized NN ===# (separate processes, or SM/pure-EFT)
"maxEventsPerClass": -1, #max nof events to be used for each process class (non-parameterized NN only) ; -1 <-> use all available events
"nEventsTot_train": -1, "nEventsTot_test": -1, #total nof events to be used for training & testing ; -1 <-> use _maxEvents & _splitTrainEventFrac params instead
"batchSizeClass": 512, #Batch size (<-> nof events fed to the network before its parameter get updated)

#=== Settings for CARL/ROLR/RASCAL strategies ===#
# "listOperatorsParam": ['ctz','ctw', 'cpqm', 'cpq3', 'cpt'], #None <-> parameterize on all possible operators
"listOperatorsParam": ['ctz'], #None <-> parameterize on all possible operators
# "listOperatorsParam": ['cpqm', 'cpt'], #None <-> parameterize on all possible operators
# "listOperatorsParam": ['ctZ','ctW', 'cpQM', 'cpQ3', 'cpt'], #None <-> parameterize on all possible operators
"nPointsPerOperator": 100, "minWC": -5, "maxWC": 5, #Interval [min,max,step] in which EFT points get sampled uniformly to train the NN on
# "listMinMaxWC": [-2,2,-2,2,-15,15,-15,15,-15,15], #If activated, and len(listMinMaxWC)=2*len(listOperatorsParam), will be interpreted as a list of min/max values for each operator selected above for NN parameterization (superseeds minWC/maxWC values)
"nEventsPerPoint": 5000, #max nof events to be used for each EFT point (for parameterized NN only) ; -1 <-> use all available events
"batchSizeEFT": 5000, #Batch size (<-> nof events fed to the network before its parameter get updated)
"refPoint": "SM", #Reference point used e.g. to compute likelihood ratios. Must be "SM" for CARL_multiclass strategy (<-> separate SM from EFT). Must be != "SM" for CARL_singlePoint strategy (<-> will correspond to the single hypothesis to separate from SM). Follow naming convention from MG, e.g.: 'ctZ_-3.5_ctp_2.6'
# "refPoint": "rwgt_ctZ_5",
# "refPoint": "rwgt_ctZ_3_ctW_0_cpQM_0_cpQ3_0_cpt_0",
"score_lossWeight": 1, #Apply scale factor to score term in loss function
"regress_onLogr": True, #True <-> NN will regress on log(r) instead of r

#=== Event preselection ===#
"cuts": "1", #Event selection, both for train/test ; "1" <-> no cut

#=== OTHERS ===#
"makeValPlotsOnly": False, #True <-> load pre-existing model, skip train/test phase, create validation plots directly. Get data first (needed for plots)
}


# Analysis options
# //--------------------------------------------

# -- Choose the data to consider #NB: same convention as for main analysis code. Naming convention enforced : 2016+2017 <-> "201617" ; etc.; 2016+2017+2018 <-> "Run2" #NB: years must be placed in the right order !
_list_lumiYears = []
# _list_lumiYears.append("2016")
_list_lumiYears.append("2017")
# _list_lumiYears.append("2018")

#-- Choose the classes of processes to consider #NB: can group several physics processes in same process class #NB: place main signal in first position
_list_processClasses = []
_list_processClasses.append(["tZq"])
_list_processClasses.append(["ttZ"])
# _list_processClasses.append(["PrivMC_tZq_fullsim"])
# _list_processClasses.append(["PrivMC_tZq_training"])
# _list_processClasses.append(["PrivMC_ttZ_training"])
# _list_processClasses.append(["PrivMC_ttZ_v3"])
# _list_processClasses.append(["PrivMC_tZq_ctz"])
# _list_processClasses.append(["PrivMC_tZq_ctw"])
# _list_processClasses.append(["PrivMC_ttZ_ctz"])
# _list_processClasses.append(["PrivMC_ttZ_ctw"])
# _list_processClasses.append(["PrivMC_tZq_ctz", "PrivMC_ttZ_ctz"])
# _list_processClasses.append(["PrivMC_tZq_ctw", "PrivMC_ttZ_ctw"])
# _list_processClasses.append(["ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep"])
# _list_processClasses.append(["ttZ", "ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep",])

#-- Define labels associated with each process class #NB: keyword 'PrivMC' is used to denote private EFT samples
_list_labels = []
_list_labels.append("tZq")
_list_labels.append("ttZ")
# _list_labels.append("PrivMC_tZq")
# _list_labels.append("PrivMC_ttZ")
# _list_labels.append("PrivMC_ttZ_top19001")
# _list_labels.append("PrivMC_tZq_ctz")
# _list_labels.append("PrivMC_tZq_ctw")
# _list_labels.append("PrivMC_ttZ_ctz")
# _list_labels.append("PrivMC_ttZ_ctw")
# _list_labels.append("Backgrounds")
# //--------------------------------------------

#-- Choose input features x
_list_features = []
_list_features.append("mTW")
_list_features.append("mHT")
_list_features.append("Mass_3l")
_list_features.append("jPrimeAbsEta")
_list_features.append("lAsymmetry")
_list_features.append("maxDelPhiLL")
_list_features.append("maxDeepCSV")
_list_features.append("deepCSV_2nd")
# '''
_list_features.append("njets")
_list_features.append("nbjets")
_list_features.append("cosThetaStarPolTop")
_list_features.append("cosThetaStarPolZ")

_list_features.append("recoZ_Pt")
_list_features.append("recoZ_Eta")
_list_features.append("recoZ_M")
_list_features.append("recoLepTop_Pt")
_list_features.append("recoLepTop_Eta")
_list_features.append("recoLepTop_M")
_list_features.append("TopZsystem_Pt")
_list_features.append("TopZsystem_M")
_list_features.append("jprime_Pt")
_list_features.append("recoLepTopLep_Pt")
_list_features.append("recoLepTopLep_Eta")

_list_features.append("maxEtaJet")
_list_features.append("mbjMax")
_list_features.append("maxDiJet_pt")
_list_features.append("maxDelRbL")
_list_features.append("dR_tZ")
_list_features.append("dR_ZlW")
_list_features.append("dR_blW")
_list_features.append("dR_bW")
_list_features.append("dR_tClosestLep")
_list_features.append("dR_jprimeClosestLep")
_list_features.append("dEta_jprimeClosestLep")
_list_features.append("dR_tClosestJet")
_list_features.append("dR_tjprime")
_list_features.append("dR_bjprime")
_list_features.append("dR_lWjprime")
_list_features.append("dR_Zjprime")
_list_features.append("maxDiJet_m")
_list_features.append("dEta_tjprime")
_list_features.append("dEta_bjprime")
_list_features.append("dEta_lWjprime")
_list_features.append("dEta_Zjprime")
# '''

# _list_features.append("maxDiJet_dPhi")
# _list_features.append("maxDiJet_dEta")
# _list_features.append("maxDiJet_dR")

# _list_features.append("recoZLepMinus_Pt")
# _list_features.append("recoZLepMinus_Eta")
# _list_features.append("recoZLepMinus_Phi")
# _list_features.append("recoZLepPlus_Pt")
# _list_features.append("recoZLepPlus_Eta")
# _list_features.append("recoZLepPlus_Phi")
# _list_features.append("recoLepTopLep_Phi")

# _list_features.append("jprime_Eta")
# _list_features.append("jprime_Phi")
# _list_features.append("TopZsystem_Eta")
# _list_features.append("TopZsystem_Phi")
# _list_features.append("recoLepTopB_Pt")
# _list_features.append("recoLepTopB_Eta")
# _list_features.append("recoLepTopB_Phi")
# _list_features.append("recoZ_Phi")
# _list_features.append("recoLepTop_Phi")

#== OR USE 'recoZLepPlus', 'recoZLepMinus', 'recoTopLep' for leptons ? etc.
# _list_features.append("lep1_pt")
# _list_features.append("lep2_pt")
# _list_features.append("lep3_pt");
# _list_features.append("lep1_eta")
# _list_features.append("lep2_eta")
# _list_features.append("lep3_eta")
# _list_features.append("lep1_phi")
# _list_features.append("lep2_phi")
# _list_features.append("lep3_phi")
# _list_features.append("jet1_pt")
# _list_features.append("jet2_pt")
# _list_features.append("jet3_pt")
# _list_features.append("jet1_eta")
# _list_features.append("jet2_eta")
# _list_features.append("jet3_eta")
# _list_features.append("jet1_phi")
# _list_features.append("jet2_phi")
# _list_features.append("jet3_phi")
# _list_features.append("jet1_DeepCSV")
# _list_features.append("jet2_DeepCSV")
# _list_features.append("jet3_DeepCSV")
# _list_features.append("jet4_pt")
# _list_features.append("jet4_eta")
# _list_features.append("jet4_phi")
# _list_features.append("jet4_DeepCSV")



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

from Utils.FreezeSession import FreezeSession_and_SaveModel #freeze_session
from Utils.Helper import *
from Utils.Model import Create_Model
from Utils.Callbacks import Get_Callbacks
from Utils.GetData import Get_Data
from Utils.LossOptimMetric import Get_Loss_Optim_Metrics
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
    _lumiName, _weightDir, _h5modelName, _ntuplesDir, _batchSize = Initialization_And_SanityChecks(optsTrain, _list_lumiYears, _list_processClasses, _list_labels)
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
    x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, shifts, scales, xTrainRescaled, _list_labels, _list_features = Get_Data(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features, _weightDir, _ntuplesDir, _lumiName)

    #-- Plot input features distributions, after applying to train data same rescaling as will be done by first NN layer (-> check rescaling)
    Plot_Input_Features(optsTrain, xTrainRescaled, y_process_train, PhysicalWeights_train, _list_features, _weightDir, True)

    print('\n'); print(colors.fg.lightblue, "--- Define the loss function & metrics...", colors.reset); print('\n')
    _loss, _optim, _metrics, _lossWeights = Get_Loss_Optim_Metrics(optsTrain)

    if optsTrain["makeValPlotsOnly"] == False: #Train, test, save, model
        #-- Get model and compile it
        print('\n'); print(colors.fg.lightblue, "--- Create the Keras model...", colors.reset); print('\n')
        model = Create_Model(optsTrain, _weightDir, _list_features, shifts, scales)

        # Can printout the output of the ith layer here for N events, e.g. to verify that the normalization layer works properly
        # Printout_Outputs_FirstLayer(model, ilayer=0, xx=x[0:5])

        print('\n'); print(colors.fg.lightblue, "--- Compile the Keras model...", colors.reset); print('\n')
        model.compile(loss=_loss, loss_weights=_lossWeights, optimizer=_optim, metrics=[_metrics]) #For multiclass classification

        #-- Define list of callbacks
        callbacks_list = Get_Callbacks(_weightDir)
        # ckpt_dir = os.path.dirname(ckpt_path); history = 0

        #-- Fit model (TRAIN)
        print('\n'); print(colors.fg.lightblue, "--- Train (fit) NN on training sample...", colors.reset); print('\n')
        if optsTrain["parameterizedNN"]==False:
            history = model.fit(x_train, y_train, sample_weight=LearningWeights_train, validation_data=(x_test, y_test, LearningWeights_test), epochs=optsTrain["nEpochs"], batch_size=_batchSize, callbacks=callbacks_list, shuffle=True, verbose=1) #CHANGED -- need to use training weights for valdata too !
            # history = model.fit(x_train, y_train, sample_weight=LearningWeights_train, validation_data=(x_test, y_test, PhysicalWeights_test), epochs=optsTrain["nEpochs"], batch_size=_batchSize, callbacks=callbacks_list, shuffle=True, verbose=1)

            # Evaluate the neural network's performance (evaluate metrics on validation or test dataset)
            print('\n'); print(colors.fg.lightblue, "--- Evaluate NN performance on test sample...", colors.reset); print('\n')
            score = model.evaluate(x_test, y_test, batch_size=_batchSize, sample_weight=PhysicalWeights_test, verbose=1)

        else:
            my_training_batch_generator = DataGenerator(x_train, y_train, LearningWeights_train, optsTrain["strategy"], _batchSize, returnWeights=False)
            my_validation_batch_generator = DataGenerator(x_test, y_test, PhysicalWeights_test, optsTrain["strategy"], _batchSize, returnWeights=False)
            _steps_per_epoch = np.ceil(len(x_train) / _batchSize); _steps_per_epoch_val = np.ceil(len(x_test)/ _batchSize)
            # batch_x, batch_y = my_training_batch_generator.__getitem__(0); print(batch_x); print(batch_y)

            history = model.fit(my_training_batch_generator, steps_per_epoch=_steps_per_epoch, validation_data=my_validation_batch_generator, validation_steps=_steps_per_epoch_val, epochs=optsTrain["nEpochs"], callbacks=callbacks_list, verbose=1)

            print('\n'); print(colors.fg.lightblue, "--- Evaluate NN performance on test sample...", colors.reset); print('\n')
            score = model.evaluate(my_validation_batch_generator, steps=_steps_per_epoch_val, verbose=1)


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
        with tensorflow.compat.v1.Session() as sess: FreezeSession_and_SaveModel(optsTrain, sess, _weightDir, _h5modelName) #Must first open a new session #Can't manage to run code below without this... (why?)

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

    #-- Get control results (printouts, plots, histos)
    list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, list_yTrain_allClasses, list_yTest_allClasses, list_xTrain_allClasses, list_xTest_allClasses = Apply_Model_toTrainTestData(optsTrain, _list_labels, x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, _h5modelName)

    #-- Store NN predictions for train/test datasets, for later use
    Store_TrainTestPrediction_Histograms(optsTrain, _lumiName, _list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses)

    #-- Create several validation plots automatically
    Make_Default_Validation_Plots(optsTrain, _list_features, _list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, PhysicalWeights_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, x, y_train, y_test, y_process, y_process_train, y_process_test, list_yTrain_allClasses, list_yTest_allClasses, list_xTrain_allClasses, list_xTest_allClasses, _metrics, _weightDir, score, history)

    Write_Timestamp_toLogfile(_weightDir, 1) #Write final timestamp before exit

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
# if __name__ == "__main__":
# parser = argparse.ArgumentParser()
# parser.add_argument("xxx", metavar="xxx", help="help")
