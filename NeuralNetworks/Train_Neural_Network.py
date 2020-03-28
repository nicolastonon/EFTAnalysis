# Nicolas TONON (DESY)
# Train fully-connected neural networks with Keras (tensorflow back-end)
# //--------------------------------------------

'''
#TODO#
- argparse (model, ...)
- Very long to load full Run2 dataset into RAM... load by batches ?
'''


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
_processClasses_list = [
                # ["PrivMC_tZq"],
                ["tZq"],
                # ["ttZ"]]
                ["ttZ"], ["ttW", "ttH", "WZ", "ZZ4l", "DY", "TTbar_DiLep"]]
                # ["ttZ", "ttW", "ttH", "WZ", "ZZ4l", "DY", "TTbar_DiLep",]]

_labels_list =  ["tZq",
                "ttZ", "Backgrounds"]
                # "Backgrounds"]

cuts = "passedBJets==1" #Event selection, both for train/test ; "1" <-> no cut
# //--------------------------------------------

#--- Training options
# //--------------------------------------------
_nepochs = 5 #Number of training epochs (<-> nof times the full training dataset is shown to the NN)
_batchSize = 512 #Batch size (<-> nof events fed to the network before its parameter get updated)
# _nof_output_nodes = 3 #1 (binary) or N (multiclass)

_maxEvents_perClass = -1 #max nof events to be used for each process ; -1 <-> all events
_nEventsTot_train = -1; _nEventsTot_test = -1  #nof events to be used for training & testing ; -1 <-> use _maxEvents_perClass & _splitTrainEventFrac params instead
_splitTrainEventFrac = 0.8 #Fraction of events to be used for training (1 <-> use all requested events for training)

# _startFromExistingModel = False #True <-> Skip training, load latest checkpoint model and create perf plots #not used yet
# //--------------------------------------------

# Define list of input variables
# //--------------------------------------------
var_list = []
var_list.append("maxDijetDelR")
var_list.append("dEtaFwdJetBJet")
var_list.append("dEtaFwdJetClosestLep")
var_list.append("mHT")
var_list.append("mTW")
var_list.append("Mass_3l")
var_list.append("forwardJetAbsEta")
var_list.append("jPrimeAbsEta")
var_list.append("maxDeepCSV")
var_list.append("delRljPrime")
var_list.append("lAsymmetry")
var_list.append("maxDijetMass")
var_list.append("maxDelPhiLL")









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
# //--------------------------------------------
import tensorflow
import keras
import numpy as np
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score
from tensorflow.keras.models import load_model

from Utils.FreezeSession import freeze_session
from Utils.Helper import batchOutput, Write_Variables_To_TextFile, TimeHistory, Get_LumiName, SanityChecks_Parameters
from Utils.Model import Create_Model
from Utils.Callbacks import Get_Callbacks
from Utils.GetData import Get_Data_For_DNN_Training
from Utils.Optimizer import Get_Loss_Optim_Metrics
from Utils.ColoredPrintout import colors
from Utils.Output_Plots_Histos import Create_TrainTest_ROC_Histos, Create_Control_Plots
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
def Train_Test_Eval_PureKeras(_lumi_years, _processClasses_list, _labels_list, var_list, cuts, _nepochs, _batchSize, _maxEvents_perClass, _splitTrainEventFrac, _nEventsTot_train, _nEventsTot_test):

 # #    # # #####
 # ##   # #   #
 # # #  # #   #
 # #  # # #   #
 # #   ## #   #
 # #    # #   #

    #Sanity chek of input args
    SanityChecks_Parameters(_processClasses_list, _labels_list)

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
    if _nof_output_nodes == 2: #Special case : 2 classes -> binary classification -> 1 output node only
        _nof_output_nodes = 1

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

    #Get data
    print(colors.fg.lightblue, "--- Read and shape the data...", colors.reset); print('\n')
    x_train, y_train, x_test, y_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, x, y, PhysicalWeights_allClasses, LearningWeights_allClasses, means, stddev, x_control_firstNEvents = Get_Data_For_DNN_Training(weight_dir, _lumi_years, _ntuples_dir, _processClasses_list, _labels_list, var_list, cuts, _nof_output_nodes, _maxEvents_perClass, _splitTrainEventFrac, _nEventsTot_train, _nEventsTot_test, lumiName)

    print('\n'); print(colors.fg.lightblue, "--- Define the loss function & metrics...", colors.reset); print('\n')
    _loss, _optim, _metrics = Get_Loss_Optim_Metrics(_nof_output_nodes)

    #Get model and compile it
    print('\n'); print(colors.fg.lightblue, "--- Create the Keras model...", colors.reset); print('\n')
    model = Create_Model(weight_dir, "DNN", _nof_output_nodes, var_list, means, stddev) #-- add default args
    print('\n'); print(colors.fg.lightblue, "--- Compile the Keras model...", colors.reset); print('\n')
    model.compile(loss=_loss, optimizer=_optim, metrics=[_metrics]) #For multiclass classification

    #Define list of callbacks
    callbacks_list = Get_Callbacks(weight_dir)
    # ckpt_dir = os.path.dirname(ckpt_path)
    history = 0

    #Fit model (TRAIN)
    print('\n'); print(colors.fg.lightblue, "--- Train (fit) DNN on training sample...", colors.reset, " (may take a while)"); print('\n')
    history = model.fit(x_train, y_train, validation_data=(x_test, y_test), epochs=_nepochs, batch_size=_batchSize, sample_weight=LearningWeights_train, callbacks=callbacks_list, shuffle=True, verbose=1)

    #-- Can access weights and biases of any layer
    # weights_layer, biases_layer = model.layers[0].get_weights(); print(weights_layer.shape); print(biases_layer.shape); print(weights_layer); print(biases_layer[0:2])

    # else:
    #     # Loads the latest checkpoint weights
    #     latest = tensorflow.train.latest_checkpoint(ckpt_dir)
    #     tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    #     model = load_model(h5model_outname) # model has to be re-loaded
    #     model.load_weights(latest)

    # Evaluate the neural network's performance (evaluate metrics on validation or test dataset)
    print('\n'); print(colors.fg.lightblue, "--- Evaluate DNN performance on test sample...", colors.reset); print('\n')
    score = model.evaluate(x_test, y_test, batch_size=_batchSize, sample_weight=PhysicalWeights_test)
    # print(score)


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

    #Save list of variables #Done in data transformation function now
    # Write_Variables_To_TextFile(weight_dir, var_list)


 ###### #####  ###### ###### ###### ######     ####  #####    ##   #####  #    #
 #      #    # #      #          #  #         #    # #    #  #  #  #    # #    #
 #####  #    # #####  #####     #   #####     #      #    # #    # #    # ######
 #      #####  #      #        #    #         #  ### #####  ###### #####  #    #
 #      #   #  #      #       #     #         #    # #   #  #    # #      #    #
 #      #    # ###### ###### ###### ######     ####  #    # #    # #      #    #

# --- Convert model to estimator and save model as frozen graph for c++

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


 #####  ######  ####  #    # #      #####  ####
 #    # #      #      #    # #        #   #
 #    # #####   ####  #    # #        #    ####
 #####  #           # #    # #        #        #
 #   #  #      #    # #    # #        #   #    #
 #    # ######  ####   ####  ######   #    ####

        print('\n\n')
        print(colors.bg.orange, colors.bold, "##############################################", colors.reset)
        print(colors.fg.orange, '\t Results & Control Plots', colors.reset)
        print(colors.bg.orange, colors.bold, "##############################################", colors.reset, '\n')

        if _nof_output_nodes == 1:
            loss = score[0]
            accuracy = score[1]
            print(colors.fg.lightgrey, '** Accuracy :', str(accuracy), colors.reset)
            print(colors.fg.lightgrey, '** Loss', str(loss), colors.reset)

        # if len(np.unique(y_train)) > 1: # prevent bug in roc_auc_score, need >=2 unique values (at least sig+bkg classes)
        #     auc_score = roc_auc_score(y_test, model.predict(x_test))
        #     auc_score_train = roc_auc_score(y_train, model.predict(x_train))
        #     print('\n'); print(colors.fg.lightgrey, '**** AUC scores ****', colors.reset)
        #     print(colors.fg.lightgrey, "-- TEST SAMPLE  \t==> " + str(auc_score), colors.reset)
        #     print(colors.fg.lightgrey, "-- TRAIN `SAMPLE \t==> " + str(auc_score_train), colors.reset); print('\n')

        #Get control results
        list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses = Apply_Model_toTrainTestData(_nof_output_nodes, _processClasses_list, _labels_list, x_train, y_train, x_test, y_test, PhysicalWeights_train, PhysicalWeights_test, h5model_outname, x_control_firstNEvents)

        Create_TrainTest_ROC_Histos(lumiName, _nof_output_nodes, _labels_list, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, _metrics)

        Create_Control_Plots(_nof_output_nodes, _labels_list, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, x_train, y_train, y_test, x_test, model, _metrics, _nof_output_nodes, weight_dir, history)

    #End [with ... as sess]
# //--------------------------------------------
# //--------------------------------------------















# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
   ###    ########  ########  ##       ##    ##
  ## ##   ##     ## ##     ## ##        ##  ##
 ##   ##  ##     ## ##     ## ##         ####
##     ## ########  ########  ##          ##
######### ##        ##        ##          ##
##     ## ##        ##        ##          ##
##     ## ##        ##        ########    ##

##     ##  #######  ########  ######## ##
###   ### ##     ## ##     ## ##       ##
#### #### ##     ## ##     ## ##       ##
## ### ## ##     ## ##     ## ######   ##
##     ## ##     ## ##     ## ##       ##
##     ## ##     ## ##     ## ##       ##
##     ##  #######  ########  ######## ########
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

def Apply_Model_toTrainTestData(nof_output_nodes, processClasses_list, labels_list, x_train, y_train, x_test, y_test, PhysicalWeights_train, PhysicalWeights_test, savedModelName, x_control_firstNEvents):

    print('\n', colors.fg.lightblue, "--- Apply model to train & test data...", colors.reset, " (may take a while)\n")

    # print('x_test:\n', x_test[:10]); print('y_test:\n', y_test[:10]); print('x_train:\n', x_train[:10]); print('y_train:\n', y_train[:10])

    list_xTrain_allClasses = []
    list_xTest_allClasses = []
    list_yTrain_allClasses = []
    list_yTest_allClasses = []
    list_PhysicalWeightsTrain_allClasses = []
    list_PhysicalWeightsTest_allClasses = []

    if nof_output_nodes == 1: #Binary
        list_xTrain_allClasses.append(x_train[y_train==1]); list_yTrain_allClasses.append(y_train[y_train==1]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_train==1])
        list_xTrain_allClasses.append(x_train[y_train==0]); list_yTrain_allClasses.append(y_train[y_train==0]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_train==0])
        list_xTest_allClasses.append(x_test[y_test==1]); list_yTest_allClasses.append(y_test[y_test==1]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_test==1])
        list_xTest_allClasses.append(x_test[y_test==0]); list_yTest_allClasses.append(y_test[y_test==0]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_test==0])
    else: #Multiclass
        for i in range(len(_processClasses_list)):
            list_xTrain_allClasses.append(x_train[y_train[:,i]==1]); list_yTrain_allClasses.append(y_train[:,i]==1); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_train[:,i]==1])
            list_xTest_allClasses.append(x_test[y_test[:,i]==1]); list_yTest_allClasses.append(y_test[:,i]==1); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_test[:,i]==1])


 #####  #####  ###### #####  #  ####  ##### #  ####  #    #  ####
 #    # #    # #      #    # # #    #   #   # #    # ##   # #
 #    # #    # #####  #    # # #        #   # #    # # #  #  ####
 #####  #####  #      #    # # #        #   # #    # #  # #      #
 #      #   #  #      #    # # #    #   #   # #    # #   ## #    #
 #      #    # ###### #####  #  ####    #   #  ####  #    #  ####

    #--- Load model
    tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    model = load_model(savedModelName)

    #Application (can also use : predict_classes, predict_proba)
    list_predictions_train_allNodes_allClasses = []
    list_predictions_test_allNodes_allClasses = []
    for inode in range(nof_output_nodes):

        list_predictions_train_allClasses = []
        list_predictions_test_allClasses = []
        for iclass in range(len(processClasses_list)):
            list_predictions_train_allClasses.append(model.predict(list_xTrain_allClasses[iclass])[:,inode])
            list_predictions_test_allClasses.append(model.predict(list_xTest_allClasses[iclass])[:,inode])

        list_predictions_train_allNodes_allClasses.append(list_predictions_train_allClasses)
        list_predictions_test_allNodes_allClasses.append(list_predictions_test_allClasses)

    # -- Printout of some predictions
    # np.set_printoptions(threshold=5) #If activated, will print full numpy arrays
    print("-------------- FEW EXAMPLES... --------------")
    for i in range(10):
        if nof_output_nodes == 1:
            if y_test[i]==1:
                true_label = "signal"
            else:
                true_label = "background"
            print("===> Prediction for %s event : %s" % (true_label, (list_predictions_test_allClasses[0])[i]))

        else:
            for j in range(len(processClasses_list)):
                if y_test[i][j]==1:
                    true_label = labels_list[j]
        print("===> Outputs nodes predictions for %s event : %s" % (true_label, (list_predictions_test_allClasses[j])[i]) )
    print("--------------\n")

    #-- Print predictions for first few events of first process => can compare with predictions obtained for same DNN/events using another code
    # for j in range(x_control_firstNEvents.shape[0]):
    #     print(x_control_firstNEvents[j])
    #     print("===> Prediction for event", j," :", model.predict(x_control_firstNEvents)[j][0], '\n')

    # return list_predictions_train_allClasses, list_predictions_test_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses
    return list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses

















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

# if len(sys.argv) == 2:
    # nLep = "3l" if sys.argv[1] == True else False

#----------  Manual call to DNN training function
Train_Test_Eval_PureKeras(_lumi_years, _processClasses_list, _labels_list, var_list, cuts, _nepochs, _batchSize, _maxEvents_perClass, _splitTrainEventFrac, _nEventsTot_train, _nEventsTot_test)
