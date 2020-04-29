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

def Apply_Model_toTrainTestData(opts, list_labels, x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, savedModelName, x_control_firstNEvents):

    # print('x_test:\n', x_test[:10]); print('y_test:\n', y_test[:10]); print('x_train:\n', x_train[:10]); print('y_train:\n', y_train[:10])

    maxEvents = 500000 #Upper limit on nof events per class, else validation too slow (problematic for parameterized DNN with huge training stat.)

    # print('...Order data by target class...')
    list_xTrain_allClasses = []; list_xTest_allClasses = []
    list_yTrain_allClasses = []; list_yTest_allClasses = []
    list_truth_Train_allClasses = []; list_truth_Test_allClasses = []
    list_PhysicalWeightsTrain_allClasses = []; list_PhysicalWeightsTest_allClasses = []
    if opts["nofOutputNodes"] == 1: #Binary
        list_xTrain_allClasses.append(x_train[y_process_train==1][:maxEvents]); list_yTrain_allClasses.append(y_train[y_process_train==1][:maxEvents]); list_truth_Train_allClasses.append(y_process_train[y_process_train==1][:maxEvents]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train==1][:maxEvents])
        list_xTrain_allClasses.append(x_train[y_process_train==0][:maxEvents]); list_yTrain_allClasses.append(y_train[y_process_train==0][:maxEvents]); list_truth_Train_allClasses.append(y_process_train[y_process_train==0][:maxEvents]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train==0][:maxEvents])
        list_xTest_allClasses.append(x_test[y_process_test==1][:maxEvents]); list_yTest_allClasses.append(y_test[y_process_test==1][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[y_process_test==1][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test==1][:maxEvents])
        list_xTest_allClasses.append(x_test[y_process_test==0][:maxEvents]); list_yTest_allClasses.append(y_test[y_process_test==0][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[y_process_test==0][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test==0][:maxEvents])
    else: #Multiclass
        for i in range(len(list_labels)):
            if opts["parameterizedDNN"] == False or i==0: #Use all events from process class for control plots (for non-parameterized DNN, and for SM point)
                list_xTrain_allClasses.append(x_train[y_process_train[:,i]==1][:maxEvents]); list_yTrain_allClasses.append(y_train[y_process_train[:,i]==1][:maxEvents]); list_truth_Train_allClasses.append(y_process_train[y_process_train[:,i]==1][:maxEvents]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train[:,i]==1][:maxEvents])
                list_xTest_allClasses.append(x_test[y_process_test[:,i]==1][:maxEvents]); list_yTest_allClasses.append(y_test[y_process_test[:,i]==1][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[y_process_test[:,i]==1][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test[:,i]==1][:maxEvents])
            else: #For each EFT operator, only use events generated at most extreme EFT point
                maxWC = opts["maxWC"]
                j = -len(opts["listOperatorsParam"]) -1 + i #Index j <-> column in array of features corresponding to current operator ==> Only select events for which this operator has its maximum WC value (extremum)
                # print(j)
                list_xTrain_allClasses.append(x_train[np.logical_and(y_process_train[:,i]==1, x_train[:,j]==maxWC)][:maxEvents]); list_yTrain_allClasses.append(y_train[np.logical_and(y_process_train[:,i]==1, x_train[:,j]==maxWC)]); list_truth_Train_allClasses.append(y_process_train[np.logical_and(y_process_train[:,i]==1, x_train[:,j]==maxWC)]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[np.logical_and(y_process_train[:,i]==1, x_train[:,j]==maxWC)][:maxEvents])
                list_xTest_allClasses.append(x_test[np.logical_and(y_process_test[:,i]==1, x_test[:,j]==maxWC)][:maxEvents]); list_yTest_allClasses.append(y_test[np.logical_and(y_process_test[:,i]==1, x_test[:,j]==maxWC)][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[np.logical_and(y_process_test[:,i]==1, x_test[:,j]==maxWC)][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[np.logical_and(y_process_test[:,i]==1, x_test[:,j]==maxWC)][:maxEvents])

    #-- Sanity checks
    # assert all(len(l) for l in list_xTrain_allClasses); assert all(len(l) for l in list_xTest_allClasses)

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
    # print('...Compute model predictions...')
    list_predictions_train_allNodes_allClasses = []
    list_predictions_test_allNodes_allClasses = []
    for inode in range(opts["nofOutputNodes"]):

        list_predictions_train_class = []
        list_predictions_test_class = []
        for iclass in range(len(list_labels)):
            # print('inode', inode, 'iclass', iclass)
            list_predictions_train_class.append(model.predict(list_xTrain_allClasses[iclass])[:,inode])
            list_predictions_test_class.append(model.predict(list_xTest_allClasses[iclass])[:,inode])

        list_predictions_train_allNodes_allClasses.append(list_predictions_train_class)
        list_predictions_test_allNodes_allClasses.append(list_predictions_test_class)

    # -- Printout of some predictions
    # np.set_printoptions(threshold=5) #If activated, will print full numpy arrays
    # print("\n-------------- FEW EXAMPLES... --------------")
    # for i in range(10):
    #     if opts["nofOutputNodes"] == 1:
    #         if y_test[i]==1:
    #             true_label = "signal"
    #         else:
    #             true_label = "background"
    #         print("===> Prediction for %s event : %s" % (true_label, (list_predictions_test_allClasses[0])[i]))
    #
    #     else:
    #         true_label = ''
    #         for j in range(len(list_labels)):
    #             if y_test[i][j]==1:
    #                 true_label = list_labels[j]
    #         print("===> Outputs nodes predictions for ", true_label, " event :", (list_predictions_test_allClasses[0])[i] )
    # print("--------------\n")

    #-- Print predictions for first few events of first process => can compare with predictions obtained for same DNN/events using another code
    # for j in range(x_control_firstNEvents.shape[0]):
    #     print(x_control_firstNEvents[j])
    #     print("===> Prediction for event", j," :", model.predict(x_control_firstNEvents)[j][0], '\n')

    return list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses
