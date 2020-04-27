#Apply the trained NN model to train and test data, to obtain its predictions (for control plots, etc.)

import numpy as np
import keras
import tensorflow
import math
from Utils.Helper import Printout_Outputs_FirstLayer
from Utils.ColoredPrintout import colors
from tensorflow.keras.models import load_model

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

def Apply_Model_toTrainTestData(nof_output_nodes, labels_list, x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, savedModelName, x_control_firstNEvents, parameterizedDNN, listOperatorsParam, maxEvents_perClass):

    # print('x_test:\n', x_test[:10]); print('y_test:\n', y_test[:10]); print('x_train:\n', x_train[:10]); print('y_train:\n', y_train[:10])

    print('...Order data by target class...')
    list_xTrain_allClasses = []; list_xTest_allClasses = []
    list_yTrain_allClasses = []; list_yTest_allClasses = []
    list_PhysicalWeightsTrain_allClasses = []; list_PhysicalWeightsTest_allClasses = []
    if nof_output_nodes == 1: #Binary
        list_xTrain_allClasses.append(x_train[y_process_train==1][:maxEvents_perClass]); list_yTrain_allClasses.append(y_train[y_process_train==1][:maxEvents_perClass]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train==1][:maxEvents_perClass])
        list_xTrain_allClasses.append(x_train[y_process_train==0][:maxEvents_perClass]); list_yTrain_allClasses.append(y_train[y_process_train==0][:maxEvents_perClass]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train==0][:maxEvents_perClass])
        list_xTest_allClasses.append(x_test[y_process_test==1][:maxEvents_perClass]); list_yTest_allClasses.append(y_test[y_process_test==1][:maxEvents_perClass]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test==1][:maxEvents_perClass])
        list_xTest_allClasses.append(x_test[y_process_test==0][:maxEvents_perClass]); list_yTest_allClasses.append(y_test[y_process_test==0][:maxEvents_perClass]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test==0][:maxEvents_perClass])
    else: #Multiclass
            for i in range(len(labels_list)):
                list_xTrain_allClasses.append(x_train[y_process_train[:,i]==1][:maxEvents_perClass]); list_yTrain_allClasses.append(y_train[y_process_train[:,i]==1][:maxEvents_perClass]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train[:,i]==1][:maxEvents_perClass])
                list_xTest_allClasses.append(x_test[y_process_test[:,i]==1][:maxEvents_perClass]); list_yTest_allClasses.append(y_test[y_process_test[:,i]==1][:maxEvents_perClass]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test[:,i]==1][:maxEvents_perClass])

    #-- Sanity checks
    assert all(len(l) for l in list_xTrain_allClasses); assert all(len(l) for l in list_xTest_allClasses)

 #####  #####  ###### #####  #  ####  ##### #  ####  #    #  ####
 #    # #    # #      #    # # #    #   #   # #    # ##   # #
 #    # #    # #####  #    # # #        #   # #    # # #  #  ####
 #####  #####  #      #    # # #        #   # #    # #  # #      #
 #      #   #  #      #    # # #    #   #   # #    # #   ## #    #
 #      #    # ###### #####  #  ####    #   #  ####  #    #  ####

    #--- Load model
    tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    model = load_model(savedModelName)

    #--- Printout : check the outputs of each layer for 2 events
    check_outputs_eachLayer = False
    if check_outputs_eachLayer == True:
        print(x_train[0:1])
        for ilayer in range(0, len(model.layers)):
            Printout_Outputs_FirstLayer(model, ilayer, x_train[0:1])

    #Application (can also use : predict_classes, predict_proba)
    print('...Compute model predictions...')
    list_predictions_train_allNodes_allClasses = []
    list_predictions_test_allNodes_allClasses = []
    for inode in range(nof_output_nodes):

        list_predictions_train_allClasses = []
        list_predictions_test_allClasses = []
        for iclass in range(len(labels_list)):
            # print('inode', inode, 'iclass', iclass)
            list_predictions_train_allClasses.append(model.predict(list_xTrain_allClasses[iclass])[:,inode])
            list_predictions_test_allClasses.append(model.predict(list_xTest_allClasses[iclass])[:,inode])

        list_predictions_train_allNodes_allClasses.append(list_predictions_train_allClasses)
        list_predictions_test_allNodes_allClasses.append(list_predictions_test_allClasses)

    # -- Printout of some predictions
    # np.set_printoptions(threshold=5) #If activated, will print full numpy arrays
    # print("\n-------------- FEW EXAMPLES... --------------")
    # for i in range(10):
    #     if nof_output_nodes == 1:
    #         if y_test[i]==1:
    #             true_label = "signal"
    #         else:
    #             true_label = "background"
    #         print("===> Prediction for %s event : %s" % (true_label, (list_predictions_test_allClasses[0])[i]))
    #
    #     else:
    #         true_label = ''
    #         for j in range(len(labels_list)):
    #             if y_test[i][j]==1:
    #                 true_label = labels_list[j]
    #         print("===> Outputs nodes predictions for ", true_label, " event :", (list_predictions_test_allClasses[0])[i] )
    # print("--------------\n")

    #-- Print predictions for first few events of first process => can compare with predictions obtained for same DNN/events using another code
    # for j in range(x_control_firstNEvents.shape[0]):
    #     print(x_control_firstNEvents[j])
    #     print("===> Prediction for event", j," :", model.predict(x_control_firstNEvents)[j][0], '\n')

    return list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses
