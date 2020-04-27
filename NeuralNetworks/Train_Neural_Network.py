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
_lumi_years = []
# _lumi_years.append("2016")
_lumi_years.append("2017")
# _lumi_years.append("2018")

#Signal process must be first
_processClasses_list = []
# _processClasses_list.append(["tZq"])
_processClasses_list.append(["PrivMC_tZq"])
# _processClasses_list.append(["PrivMC_tZq_ctz"])
# _processClasses_list.append(["PrivMC_tZq_ctw"])
# _processClasses_list.append(["ttZ"])
# _processClasses_list.append(["ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep"])
# _processClasses_list.append(["ttZ", "ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep",])

_labels_list = []
# _labels_list.append("tZq")
_labels_list.append("PrivMC_tZq")
# _labels_list.append("PrivMC_tZq_ctz")
# _labels_list.append("PrivMC_tZq_ctw")
# _labels_list.append("ttZ")
# _labels_list.append("Backgrounds")

cuts = "passedBJets==1" #Event selection, both for train/test ; "1" <-> no cut
# //--------------------------------------------

#--- Training options
# //--------------------------------------------
_regress = False #True <-> DNN used for regression ; False <-> classification
_target = "class"

_parameterizedDNN = True #True <-> DNN is parameterized on the Wilson coefficients of the EFT operators
_listOperatorsParam = ['ctZ','ctW'] #None <-> parameterize on all possible operators
# _listOperatorsParam = ['ctZ','ctW', 'cpQM', 'cpQ3', 'cpt'] #None <-> parameterize on all possible operators

_nepochs = 1 #Number of training epochs (<-> nof times the full training dataset is shown to the NN)
_batchSize = 50000 #Batch size (<-> nof events fed to the network before its parameter get updated)

_maxEvents_perClass = -1 #max nof events to be used for each process class ; -1 <-> all events
_nEventsTot_train = -1; _nEventsTot_test = -1  #nof events to be used for training & testing ; -1 <-> use _maxEvents_perClass & _splitTrainEventFrac params instead
_splitTrainEventFrac = 0.8 #Fraction of events to be used for training (1 <-> use all requested events for training)

# _startFromExistingModel = False #True <-> Skip training, load latest checkpoint model and create perf plots #not used yet
# //--------------------------------------------

# Define list of input variables
# //--------------------------------------------
_var_list = []
_var_list.append("maxDijetDelR")
_var_list.append("dEtaFwdJetBJet")
_var_list.append("dEtaFwdJetClosestLep")
_var_list.append("mHT")
_var_list.append("mTW")
_var_list.append("Mass_3l")
_var_list.append("forwardJetAbsEta")
_var_list.append("jPrimeAbsEta")
_var_list.append("maxDeepCSV")
_var_list.append("delRljPrime")
_var_list.append("lAsymmetry")
_var_list.append("maxDijetMass")
_var_list.append("maxDelPhiLL")









# //--------------------------------------------
#Filtering out manually some unimportant warnings
# import warnings
# warnings.filterwarnings("ignore", message="tensorflow:sample_weight modes were coerced")

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os
import argparse
# //--------------------------------------------
import tensorflow
import keras
import numpy as np
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score
from tensorflow.keras.models import load_model

from Utils.FreezeSession import freeze_session
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
def Train_Test_Eval_DNN(_regress, _target, _parameterizedDNN, _listOperatorsParam, _lumi_years, _processClasses_list, _labels_list, _var_list, cuts, _nepochs, _batchSize, _maxEvents_perClass, _splitTrainEventFrac, _nEventsTot_train, _nEventsTot_test):

 # #    # # #####
 # ##   # #   #
 # # #  # #   #
 # #  # # #   #
 # #   ## #   #
 # #    # #   #

    #Sanity chek of input args
    _parameterizedDNN = SanityChecks_Parameters(_processClasses_list, _labels_list, _regress, _target, _parameterizedDNN, _listOperatorsParam)

    #Read luminosity choice
    lumiName = Get_LumiName(_lumi_years)

    # Set main output paths
    weight_dir = "../weights/DNN/" + lumiName + '/'
    os.makedirs(weight_dir, exist_ok=True)

    #Top directory containing all input ntuples
    _ntuples_dir = "../input_ntuples/"

    #Model output name
    h5model_outname = weight_dir + 'model.h5'

    #Determine/store number of process classes
    _nof_output_nodes = len(_processClasses_list) #1 output node per class
    if _regress or _nof_output_nodes == 2: #Special case : 2 classes -> binary classification -> 1 output node only
        _nof_output_nodes = 1
    if _parameterizedDNN:
        _nof_output_nodes = len(_listOperatorsParam)+1 #1 output node for SM and each EFT operator

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
    _transfType = 'quantile' #Feature norm. method -- 'range', 'gauss', 'quantile'
    x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, shifts, scales, x_control_firstNEvents, xTrainRescaled = Get_Data(_regress, weight_dir, _lumi_years, _ntuples_dir, _processClasses_list, _labels_list, _var_list, cuts, _nof_output_nodes, _maxEvents_perClass, _splitTrainEventFrac, _nEventsTot_train, _nEventsTot_test, lumiName, _parameterizedDNN, _listOperatorsParam, _transfType)

    # From there, for parameterized DNN, the different 'classes' correspond to different operators, and must include their WCs as inputs
    if _parameterizedDNN==True :
        _labels_list = _listOperatorsParam[:]; _labels_list.insert(0, "SM") #Specify '[:]' to create a copy, not a reference
        _var_list = np.append(_var_list, _listOperatorsParam)

    #-- Plot input features distributions, after applying to train data same rescaling as will be done by first DNN layer (-> check rescaling)
    Plot_Input_Features(xTrainRescaled, y_process_train, PhysicalWeights_train, _var_list, weight_dir, _nof_output_nodes, _parameterizedDNN, True)

    print('\n'); print(colors.fg.lightblue, "--- Define the loss function & metrics...", colors.reset); print('\n')
    _loss, _optim, _metrics = Get_Loss_Optim_Metrics(_regress, _nof_output_nodes)

    #-- Get model and compile it
    print('\n'); print(colors.fg.lightblue, "--- Create the Keras model...", colors.reset); print('\n')
    model = Create_Model(_regress, weight_dir, _nof_output_nodes, _var_list, shifts, scales, _parameterizedDNN, _listOperatorsParam) #-- add default args

    # Can printout the output of the ith layer here for N events, e.g. to verify that the normalization layer works properly
    # Printout_Outputs_FirstLayer(model, ilayer=0, xx=x[0:5])

    print('\n'); print(colors.fg.lightblue, "--- Compile the Keras model...", colors.reset); print('\n')
    model.compile(loss=_loss, optimizer=_optim, metrics=[_metrics]) #For multiclass classification

    #-- Define list of callbacks
    callbacks_list = Get_Callbacks(weight_dir)
    # ckpt_dir = os.path.dirname(ckpt_path); history = 0

    #-- Fit model (TRAIN)
    print('\n'); print(colors.fg.lightblue, "--- Train (fit) DNN on training sample...", colors.reset, " (may take a while)"); print('\n')
    if _parameterizedDNN==False:
        history = model.fit(x_train, y_train, sample_weight=LearningWeights_train, validation_data=(x_test, y_test, PhysicalWeights_test), epochs=_nepochs, batch_size=_batchSize, callbacks=callbacks_list, shuffle=True, verbose=1)

        # Evaluate the neural network's performance (evaluate metrics on validation or test dataset)
        print('\n'); print(colors.fg.lightblue, "--- Evaluate DNN performance on test sample...", colors.reset); print('\n')
        score = model.evaluate(x_test, y_test, batch_size=_batchSize, sample_weight=PhysicalWeights_test, verbose=1)

    else:
        my_training_batch_generator = DataGenerator(x_train, y_train, LearningWeights_train, _batchSize)
        my_validation_batch_generator = DataGenerator(x_test, y_test, PhysicalWeights_test, _batchSize)
        _steps_per_epoch = np.ceil(len(x_train) / _batchSize); _steps_per_epoch_val = np.ceil(len(x_test)/ _batchSize)
        history = model.fit(my_training_batch_generator, steps_per_epoch=_steps_per_epoch, validation_data=my_validation_batch_generator, validation_steps=_steps_per_epoch_val, epochs=_nepochs, callbacks=callbacks_list, verbose=1)

        print('\n'); print(colors.fg.lightblue, "--- Evaluate DNN performance on test sample...", colors.reset); print('\n')
        score = model.evaluate(my_validation_batch_generator, steps=_steps_per_epoch_val, verbose=1)

# //--------------------------------------------
    #-- Can access weights and biases of any layer
    # weights_layer, biases_layer = model.layers[0].get_weights(); print(weights_layer.shape); print(biases_layer.shape); print(weights_layer); print(biases_layer[0:2])
    #-- Loads the latest checkpoint weights
    # latest = tensorflow.train.latest_checkpoint(ckpt_dir)
    # tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    # model = load_model(h5model_outname) # model has to be re-loaded
    # model.load_weights(latest)
# //--------------------------------------------


  ####    ##   #    # ######    #    #  ####  #####  ###### #
 #       #  #  #    # #         ##  ## #    # #    # #      #
  ####  #    # #    # #####     # ## # #    # #    # #####  #
      # ###### #    # #         #    # #    # #    # #      #
 #    # #    #  #  #  #         #    # #    # #    # #      #
  ####  #    #   ##   ######    #    #  ####  #####  ###### ######

    print('\n'); print(colors.fg.lightblue, "--- Save model...", colors.reset);

    #Serialize model to HDF5
    model.save(h5model_outname)

    # Save the model architecture
    with open(weight_dir + 'arch_DNN.json', 'w') as json_file:
        json_file.write(model.to_json())


 ###### #####  ###### ###### ###### ######     ####  #####    ##   #####  #    #
 #      #    # #      #          #  #         #    # #    #  #  #  #    # #    #
 #####  #    # #####  #####     #   #####     #      #    # #    # #    # ######
 #      #####  #      #        #    #         #  ### #####  ###### #####  #    #
 #      #   #  #      #       #     #         #    # #   #  #    # #      #    #
 #      #    # ###### ###### ###### ######     ####  #    # #    # #      #    #

# --- Convert model to estimator and save model as frozen graph for c++

    #FIXME -- can move to other class ?
    with tensorflow.compat.v1.Session() as sess: #Must first open a new session #Can't manage to run code below without this... (why?)

        print('\n'); print(colors.fg.lightblue, "--- Freeze graph...", colors.reset); print('\n')

        tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
        model = load_model(h5model_outname) # model has to be re-loaded

        # tensorflow.compat.v1.keras.backend.clear_session() #Closing the last session avoids that node names get a suffix appened when opening a new session #Does not work?
        # sess = tensorflow.compat.v1.keras.backend.get_session()
        # graph = sess.graph

        inputs_names = [input.op.name for input in model.inputs]
        outputs_names = [output.op.name for output in model.outputs]
        # print('\ninputs: ', model.inputs)
        print('\n')
        print(colors.fg.lightgrey, '--> inputs_names: ', inputs_names[0], colors.reset, '\n')
        # print('\noutputs: ', model.outputs)
        print(colors.fg.lightgrey, '--> outputs_names: ', outputs_names[0], colors.reset, '\n')
        # tf_node_list = [n.name for n in  tensorflow.compat.v1.get_default_graph().as_graph_def().node]; print('nodes list : ', tf_node_list)
        frozen_graph = freeze_session(sess, output_names=[output.op.name for output in model.outputs])
        tensorflow.io.write_graph(frozen_graph, weight_dir, 'model.pbtxt', as_text=True)
        tensorflow.io.write_graph(frozen_graph, weight_dir, 'model.pb', as_text=False)
        print('\n'); print(colors.fg.lightgrey, '===> Successfully froze graph :', colors.reset, weight_dir+'model.pb', '\n')

        #Also append the names of the input/output nodes in the file "DNN_info.txt" containing input features names, etc. (for later use in C++ code)
        text_file = open(weight_dir + "DNN_infos.txt", "a") #Append mode
        text_file.write(inputs_names[0]); text_file.write(' -1 -1 \n'); #use end values as flags to signal these lines
        text_file.write(outputs_names[0]); text_file.write(' -2 -2 \n');
        text_file.write(str(_nof_output_nodes)); text_file.write(' -3 -3 \n');
        text_file.close()


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

        if _nof_output_nodes == 1:
            loss = score[0]
            accuracy = score[1]
            print(colors.fg.lightgrey, '** Loss :', float('%.4g' % loss), colors.reset)
            print(colors.fg.lightgrey, '** Metrics :', float('%.4g' % accuracy), colors.reset)

        # if len(np.unique(y_train)) > 1: # prevent bug in roc_auc_score, need >=2 unique values (at least sig+bkg classes)
        #     auc_score = roc_auc_score(y_test, model.predict(x_test))
        #     auc_score_train = roc_auc_score(y_train, model.predict(x_train))
        #     print('\n'); print(colors.fg.lightgrey, '**** AUC scores ****', colors.reset)
        #     print(colors.fg.lightgrey, "-- TEST SAMPLE  \t==> " + str(auc_score), colors.reset)
        #     print(colors.fg.lightgrey, "-- TRAIN `SAMPLE \t==> " + str(auc_score_train), colors.reset); print('\n')

# //--------------------------------------------
        #-- Get control results (printouts, plots, histos)
        print('\n', colors.fg.lightblue, "--- Apply model to train & test data...", colors.reset, " (may take a while)\n")
        list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses = Apply_Model_toTrainTestData(_nof_output_nodes, _labels_list, x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, h5model_outname, x_control_firstNEvents, _parameterizedDNN, _listOperatorsParam, _maxEvents_perClass)

        Control_Printouts(_nof_output_nodes, _labels_list, y_test, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses)

        Create_TrainTest_ROC_Histos(lumiName, _nof_output_nodes, _labels_list, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, _metrics)

        Create_Control_Plots(_regress, _nof_output_nodes, _labels_list, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, x_train, x_test, y_train, y_test, y_process_train, y_process_test, model, _metrics, _nof_output_nodes, weight_dir, history)

        Create_Correlation_Plot(x, _var_list, weight_dir)

        Plot_Input_Features(x, y_process, PhysicalWeights_allClasses, _var_list, weight_dir, _nof_output_nodes, _parameterizedDNN, False)
# //--------------------------------------------

    #End [with ... as sess]
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
Train_Test_Eval_DNN(_regress, _target, _parameterizedDNN, _listOperatorsParam, _lumi_years, _processClasses_list, _labels_list, _var_list, cuts, _nepochs, _batchSize, _maxEvents_perClass, _splitTrainEventFrac, _nEventsTot_train, _nEventsTot_test)

# //--------------------------------------------
#-- Set up the command line arguments
# if __name__ == "__main__":
# parser = argparse.ArgumentParser()
# parser.add_argument("xxx", metavar="xxx", help="help")
