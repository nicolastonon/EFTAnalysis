#Apply the trained NN model to train and test data, to obtain its predictions (for control plots, etc.)

import numpy as np
# import keras
import tensorflow as tf
import math
from Utils.Helper import Printout_Outputs_Layer
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

def Apply_Model_toTrainTestData(opts, weightDir, list_processClasses, list_labels, x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, savedModelName):

    # print('x_test:\n', x_test[:10]); print('y_test:\n', y_test[:10]); print('x_train:\n', x_train[:10]); print('y_train:\n', y_train[:10])

    maxEvents = 100000 #Upper limit on nof events per class, else validation too slow (problematic for parameterized NN with huge training stat.)

    useMostExtremeWCvaluesOnly = True #True <-> for 'EFT' class, will only consider points generated at the most extreme WC values included during training (not all the intermediate points) #Use this to create more "representative" val plots, in which only a few specific WC values are included instead of all points

# //--------------------------------------------

    minWC = opts["minWC"]
    maxWC = opts["maxWC"]

    if opts["parameterizedNN"] == False or "listMinMaxWC" in opts: useMostExtremeWCvaluesOnly = False #Can only choose to consider 'extreme' (min/max) WC values for parameterized NN (since we cut on the corresponding input features)

    #-- Store the training and testing events into lists, depending on their 'true class' values
    #NB: For 'EFT' class, don't just consider all non-SM events (too much). For now, only consider events generated at boundary values (for all operators)
    #NB: the order of the 'append' commands matters ! (determines which class will be first, impacting plotting and such)
    list_xTrain_allClasses = []; list_xTest_allClasses = []
    list_yTrain_allClasses = []; list_yTest_allClasses = []
    list_truth_Train_allClasses = []; list_truth_Test_allClasses = []
    list_PhysicalWeightsTrain_allClasses = []; list_PhysicalWeightsTest_allClasses = []
    if opts["nofOutputNodes"] == 1 or opts["regress"] == True: #Binary class label

        if useMostExtremeWCvaluesOnly is False: #Get predictions for all events
            list_xTrain_allClasses.append(x_train[y_process_train==1][:maxEvents]); list_yTrain_allClasses.append(y_train[y_process_train==1][:maxEvents]); list_truth_Train_allClasses.append(y_process_train[y_process_train==1][:maxEvents]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train==1][:maxEvents])
            list_xTrain_allClasses.append(x_train[y_process_train==0][:maxEvents]); list_yTrain_allClasses.append(y_train[y_process_train==0][:maxEvents]); list_truth_Train_allClasses.append(y_process_train[y_process_train==0][:maxEvents]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train==0][:maxEvents])

            list_xTest_allClasses.append(x_test[y_process_test==1][:maxEvents]); list_yTest_allClasses.append(y_test[y_process_test==1][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[y_process_test==1][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test==1][:maxEvents])
            list_xTest_allClasses.append(x_test[y_process_test==0][:maxEvents]); list_yTest_allClasses.append(y_test[y_process_test==0][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[y_process_test==0][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test==0][:maxEvents])
        else:
            #-----
            #Get relevant indices or slices: only retain EFT events whose operators are at boundaries
            sl = np.s_[x_train.shape[1]-len(opts["listOperatorsParam"]):] #Specify slice corresponding to indices of theory parameters features (placed last)

            # indices_train_class0 = np.where(y_process_train==0) #All events drawn at SM point #CHANGED -- must evaluate SM events at same points as EFT !
            indices_train_class0 = np.where(np.logical_and(y_process_train==0, np.logical_or(np.all(x_train[:,sl]==minWC,axis=1),np.all(x_train[:,sl]==maxWC,axis=1))) ) #All events drawn at SM point corresponding to min or max boundaries (arbitrary choice)
            indices_train_class1 = np.where(np.logical_and(y_process_train==1, np.logical_or(np.all(x_train[:,sl]==minWC,axis=1),np.all(x_train[:,sl]==maxWC,axis=1))) ) #All events drawn at EFT point corresponding to min or max boundaries (arbitrary choice)

            # indices_test_class0 = np.where(y_process_test==0) #All events drawn at SM point #CHANGED -- must evaluate SM events at same points as EFT !
            indices_test_class0 = np.where(np.logical_and(y_process_test==0, np.logical_or(np.all(x_test[:,sl]==minWC,axis=1),np.all(x_test[:,sl]==maxWC,axis=1))) ) #All events drawn at EFT point corresponding to min or max boundaries (arbitrary choice)
            indices_test_class1 = np.where(np.logical_and(y_process_test==1, np.logical_or(np.all(x_test[:,sl]==minWC,axis=1),np.all(x_test[:,sl]==maxWC,axis=1))) ) #All events drawn at EFT point corresponding to min or max boundaries (arbitrary choice)
            # print(indices_train_class0); print(indices_train_class1); print(indices_test_class0); print(indices_test_class1)
            #-----

            list_xTrain_allClasses.append(x_train[indices_train_class1][:maxEvents]); list_yTrain_allClasses.append(y_train[indices_train_class1][:maxEvents]); list_truth_Train_allClasses.append(y_process_train[indices_train_class1][:maxEvents]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[indices_train_class1][:maxEvents])
            list_xTrain_allClasses.append(x_train[indices_train_class0][:maxEvents]); list_yTrain_allClasses.append(y_train[indices_train_class0][:maxEvents]); list_truth_Train_allClasses.append(y_process_train[indices_train_class0][:maxEvents]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[indices_train_class0][:maxEvents])

            list_xTest_allClasses.append(x_test[indices_test_class1][:maxEvents]); list_yTest_allClasses.append(y_test[indices_test_class1][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[indices_test_class1][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[indices_test_class1][:maxEvents])
            list_xTest_allClasses.append(x_test[indices_test_class0][:maxEvents]); list_yTest_allClasses.append(y_test[indices_test_class0][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[indices_test_class0][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[indices_test_class0][:maxEvents])

    else: #Multi-output

        for inode in range(len(list_labels)): #For each node

            if opts["parameterizedNN"] == False or inode==0 or useMostExtremeWCvaluesOnly == False: #Use all events for control plots (for non-parameterized NN, and for SM point). ((Specify 'or inode==0' so that 'SM' events fall in this condition and not the next for 'CARL_multiclass' strategy))
                list_xTrain_allClasses.append(x_train[y_process_train[:,inode]==1][:maxEvents]); list_yTrain_allClasses.append(y_train[y_process_train[:,inode]==1][:maxEvents]); list_truth_Train_allClasses.append(y_process_train[y_process_train[:,inode]==1][:maxEvents]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[y_process_train[:,inode]==1][:maxEvents])
                list_xTest_allClasses.append(x_test[y_process_test[:,inode]==1][:maxEvents]); list_yTest_allClasses.append(y_test[y_process_test[:,inode]==1][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[y_process_test[:,inode]==1][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[y_process_test[:,inode]==1][:maxEvents])
            elif opts["strategy"] is "CARL_multiclass": #For each EFT operator (not SM = first label!), only use events generated at most extreme EFT point
                j = -len(opts["listOperatorsParam"]) -1 + inode #Index j <-> column in array of features corresponding to current operator ==> Only select events for which this operator has its maximum WC value (extremum)
                list_xTrain_allClasses.append(x_train[np.logical_and(y_process_train[:,inode]==1, np.logical_or(x_train[:,j]==minWC,x_train[:,j]==maxWC))][:maxEvents]); list_yTrain_allClasses.append(y_train[np.logical_and(y_process_train[:,inode]==1, np.logical_or(x_train[:,j]==minWC,x_train[:,j]==maxWC))]); list_truth_Train_allClasses.append(y_process_train[np.logical_and(y_process_train[:,inode]==1, np.logical_or(x_train[:,j]==minWC,x_train[:,j]==maxWC))]); list_PhysicalWeightsTrain_allClasses.append(PhysicalWeights_train[np.logical_and(y_process_train[:,inode]==1, np.logical_or(x_train[:,j]==minWC,x_train[:,j]==maxWC))][:maxEvents])
                list_xTest_allClasses.append(x_test[np.logical_and(y_process_test[:,inode]==1, np.logical_or(x_test[:,j]==minWC,x_test[:,j]==maxWC))][:maxEvents]); list_yTest_allClasses.append(y_test[np.logical_and(y_process_test[:,inode]==1, np.logical_or(x_test[:,j]==minWC,x_test[:,j]==maxWC))][:maxEvents]); list_truth_Test_allClasses.append(y_process_test[np.logical_and(y_process_test[:,inode]==1, np.logical_or(x_test[:,j]==minWC,x_test[:,j]==maxWC))][:maxEvents]); list_PhysicalWeightsTest_allClasses.append(PhysicalWeights_test[np.logical_and(y_process_test[:,inode]==1, np.logical_or(x_test[:,j]==minWC,x_test[:,j]==maxWC))][:maxEvents])

    #-- Sanity checks: make sure no class is empty
    #NB: if 1 of the 'assert' fails, it may be because the train/test datasets do not contain 1 process class, e.g. because the data were not shuffled before splitting
    if opts["trainAtManyEFTpoints"] is True or len(list_processClasses)>1: #E.g. for a regressor trained on a single sample, the 'bkg' class will be empty <-> don't check
        assert all(len(l) for l in list_xTrain_allClasses)
        assert all(len(l) for l in list_xTest_allClasses)

    # print(list_xTest_allClasses[0][:15])


 #####  #####  ###### #####  #  ####  ##### #  ####  #    #  ####
 #    # #    # #      #    # # #    #   #   # #    # ##   # #
 #    # #    # #####  #    # # #        #   # #    # # #  #  ####
 #####  #####  #      #    # # #        #   # #    # #  # #      #
 #      #   #  #      #    # # #    #   #   # #    # #   ## #    #
 #      #    # ###### #####  #  ####    #   #  ####  #    #  ####

    #--- Load model
    tf.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    model = load_model(savedModelName, compile=False) #compile=False <-> does not need to define any custom loss, since not needed for testing

    # for i in range(len(list_xTest_allClasses[1][-100:] )):
    #     if list_yTest_allClasses[1][i] > 15:
    #         print(list_xTest_allClasses[1][i])
    #         print(list_yTest_allClasses[1][i])
    #         print(list_xTest_allClasses[1].shape)
    #         print(model.predict(list_xTest_allClasses[1])[i,:])

    #--- Printout : check the outputs of each layer for 2 events
    # for ilayer in range(0, len(model.layers)):
    #     Printout_Outputs_Layer(model, ilayer, x_train[:1])

    #-- Debug: check model predictions for different WC values
    # list_xTest_allClasses[0][:,-len(opts["listOperatorsParam"]):] = 0
    # print(model.predict(list_xTest_allClasses[0][:10,:]))
    # list_xTest_allClasses[0][:,-len(opts["listOperatorsParam"]):] = 5
    # print(model.predict(list_xTest_allClasses[0][:10,:]))
    # print(list_yTest_allClasses[0][:15])
    # print(model.predict(list_xTest_allClasses[0][:15]))

    #-- Get model predictions for all events, as found in lists created above #Other available functions are: predict_classes, predict_proba, ...
    #-- Store predictions in 2D lists; 1st element = node ; 2nd element = class (e.g. 'tZq'/'ttZ', or 'SM'/'EFT')
    #-- NB: from there, all automated validation plots should use these lists ! (<-> consistent event ordering in all lists, by construction)
    list_predictions_train_allNodes_allClasses = []; list_predictions_test_allNodes_allClasses = []
    for inode in range(opts["nofOutputNodes"]):
        list_predictions_train_class = []; list_predictions_test_class = []
        for iclass in range(len(list_labels)):
            # print('inode', inode, 'iclass', iclass)

            if opts["nofOutputNodes"] == 1: #Single output node
            # if opts["strategy"] is "regressor" and iclass>0: continue #hard-coded fix, because in that case the different 'labels' already correspond to different nodes... #need to change labels/nodes ?

                list_predictions_train_class.append(np.squeeze(model.predict(list_xTrain_allClasses[iclass])) ) #Squeeze: 2D (single column) -> 1D
                list_predictions_test_class.append(np.squeeze(model.predict(list_xTest_allClasses[iclass])) )

            else: #Multinode
                if opts["strategy"] is "RASCAL": #For RASCAL there are 2 asymmetric outputs : 1 single output for r, and 1 output with N sub-outputs for t
                    if inode == 0:
                        list_predictions_train_class.append(np.squeeze(model.predict(list_xTrain_allClasses[iclass])[0]))
                        list_predictions_test_class.append(np.squeeze(model.predict(list_xTest_allClasses[iclass])[0]))
                    else:
                        list_predictions_train_class.append(model.predict(list_xTrain_allClasses[iclass])[1][:,inode-1])
                        list_predictions_test_class.append(model.predict(list_xTest_allClasses[iclass])[1][:,inode-1])
                else:
                    list_predictions_train_class.append(model.predict(list_xTrain_allClasses[iclass])[:,inode])
                    list_predictions_test_class.append(model.predict(list_xTest_allClasses[iclass])[:,inode])

        list_predictions_train_allNodes_allClasses.append(list_predictions_train_class)
        list_predictions_test_allNodes_allClasses.append(list_predictions_test_class)

    #-- Modify all list elements at once
    # list_yTest_allClasses[:] = [[0.5] * len(inner) for inner in list_yTest_allClasses[:]]

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

    #-- NEW: also append to NN info file the min/max output values of the NN based on the train+test evaluation datasets (--> Can later adapt accordingly the template histogram boundaries to avoid empty bins)
    text_file = open(weightDir + "NN_info.txt", "a+") #Append mode
    for inode in range(opts["nofOutputNodes"]):
        # print(np.concatenate((list_predictions_train_allNodes_allClasses[inode],list_predictions_test_allNodes_allClasses[inode])).shape)
        min_prediction = np.concatenate((np.concatenate(list_predictions_train_allNodes_allClasses[inode]),np.concatenate(list_predictions_test_allNodes_allClasses[inode]))).min(axis=0)
        max_prediction = np.concatenate((np.concatenate(list_predictions_train_allNodes_allClasses[inode]),np.concatenate(list_predictions_test_allNodes_allClasses[inode]))).max(axis=0)
        # min_prediction = np.concatenate(np.concatenate((list_predictions_train_allNodes_allClasses[inode],list_predictions_test_allNodes_allClasses[inode]))).min(axis=0)
        # max_prediction = np.concatenate(np.concatenate((list_predictions_train_allNodes_allClasses[inode],list_predictions_test_allNodes_allClasses[inode]))).max(axis=0)
        # print('min_prediction =', min_prediction)
        # print('max_prediction =', max_prediction)
        text_file.write('bounds'); text_file.write(' ' + str(min_prediction)); text_file.write(' ' + str(max_prediction)+'\n');
    text_file.close()

    return list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, list_yTrain_allClasses, list_yTest_allClasses, list_xTrain_allClasses, list_xTest_allClasses
