'''
Create validation plots: overtraining distributions, loss, metrics, ROCs, input correlations, etc.

NB: ROC curves don't take event weights into account --> Not correct ! (should get ROC from scan of fine-bined TH1...)
'''

import os
import ROOT
import numpy as np
# import keras
import math
import json
import shap
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from tensorflow.keras.models import load_model, model_from_json
from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
from root_numpy import fill_hist
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score, classification_report, confusion_matrix, multilabel_confusion_matrix
# from sklearn.feature_selection import mutual_info_classif
from Utils.Helper import close_event, KS_test, Anderson_Darling_test, ChiSquare_test, s_from_r, r_from_s
from Utils.ColoredPrintout import colors
from Utils.RegressorValidation import *
from pandas.plotting import scatter_matrix
# from ann_visualizer.visualize import ann_viz #Requires installation
from tensorflow.keras.utils import plot_model
from scipy import optimize

# //--------------------------------------------
# //--------------------------------------------
 ######  ########  #######  ########  ########
##    ##    ##    ##     ## ##     ## ##
##          ##    ##     ## ##     ## ##
 ######     ##    ##     ## ########  ######
      ##    ##    ##     ## ##   ##   ##
##    ##    ##    ##     ## ##    ##  ##
 ######     ##     #######  ##     ## ########

########  ########  ######## ########  ####  ######  ######## ####  #######  ##    ##
##     ## ##     ## ##       ##     ##  ##  ##    ##    ##     ##  ##     ## ###   ##
##     ## ##     ## ##       ##     ##  ##  ##          ##     ##  ##     ## ####  ##
########  ########  ######   ##     ##  ##  ##          ##     ##  ##     ## ## ## ##
##        ##   ##   ##       ##     ##  ##  ##          ##     ##  ##     ## ##  ####
##        ##    ##  ##       ##     ##  ##  ##    ##    ##     ##  ##     ## ##   ###
##        ##     ## ######## ########  ####  ######     ##    ####  #######  ##    ##
# //--------------------------------------------
# //--------------------------------------------

def Store_TrainTestPrediction_Histograms(opts, lumiName, list_features, list_labels, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTest_allClasses, list_xTest_allClasses, list_predictions_train_allNodes_allClasses=[], list_PhysicalWeightsTrain_allClasses=[], list_xTrain_allClasses=[], scan=False, op='', WC=''):
    '''
    Apply NN model on train/test datasets to produce ROOT histograms which can later be used to plot ROC curves.

    NB: normalize to avoid effect from different nof entries, ... ?
    '''

    if opts["strategy"] in ["CASCAL"]: return #Not supported yet

    print(colors.fg.lightblue, "\n--- Create & store ROC histos...", colors.reset)

    compare_ROC_inputFeature = 'recoZ_Pt' #If not '', will also store histogram of corresponding feature, so that its ROC curve can be compareed (will only work if the feature displays a left/right separation)
    # compare_ROC_inputFeature = 'recoZ_Eta'
    # compare_ROC_inputFeature = 'lep3_phi'
    # compare_ROC_inputFeature = 'dEta_tjprime'

    xmin_feature = -1; xmax_feature = -1 #Can hard-code feature histogram range here #-1 : use ranges predefined below

    nbins = 100
    if xmin_feature==-1 and xmax_feature==-1:
        if compare_ROC_inputFeature=="recoZ_Pt":
            xmin_feature = 0; xmax_feature=500
        elif compare_ROC_inputFeature=="recoZ_Eta":
            xmin_feature = -5; xmax_feature=5
        else: #Default
            xmin_feature = 0; xmax_feature=500

    maxEvents = 500000 #Upper limit on nof events per class, else validation too slow (problematic for parameterized NN with huge training stat.)

    store_trainHisto = True
    if len(list_predictions_train_allNodes_allClasses)==0: store_trainHisto = False #Only considering 'test' events

    idx_compare_ROC_inputFeature = -1
    for idx,feature in enumerate(list_features):
        if feature==compare_ROC_inputFeature: idx_compare_ROC_inputFeature = idx
        # print(feature, idx_compare_ROC_inputFeature)

    # Fill a ROOT histogram from NumPy arrays, with fine binning (no loss of info)
    signame_tmp = list_labels[0]
    if signame_tmp == "SM": signame_tmp = "EFT" #Special case: in StdVal code, SM point is first... but by convention, want to have signal=EFT (first position)
    rootfile_outname = "../outputs/NN_"+signame_tmp+"_"+lumiName+".root"
    if scan: rootfile_outname = "../outputs/NN_"+signame_tmp+"_"+lumiName+"_"+op+"_"+WC+".root"
    fout = ROOT.TFile(rootfile_outname, "RECREATE")

    # Comparison with input feature: store in separate file
    if idx_compare_ROC_inputFeature>=0:
        rootfile_outname_feature = "../outputs/NN_"+signame_tmp+"_"+lumiName+'_'+list_features[idx_compare_ROC_inputFeature]+'.root'
        if scan: rootfile_outname_feature = "../outputs/NN_"+signame_tmp+"_"+lumiName+'_'+list_features[idx_compare_ROC_inputFeature]+"_"+op+WC+".root"
        if path.exists(rootfile_outname_feature): os.remove(rootfile_outname_feature) #Remove existing file (to avoid amconfusion)
        fout_feature = ROOT.TFile(rootfile_outname_feature, "RECREATE")

    nofOutputNodes = opts["nofOutputNodes"]

    nodes_labels = list_labels
    if opts["strategy"] is "RASCAL":
        nodes_labels = ["LR"]
        for op in opts["listOperatorsParam"]: nodes_labels.append(str("Score "+op))
    elif opts["strategy"] is "regressor":
        if opts["targetVarIdx"][0] >=0:
            nodes_labels = []; [nodes_labels.append(list_features[v]) for v in opts["targetVarIdx"]] #1 node per target variable
        else: nodes_labels = ['target'] #Default target
    elif nofOutputNodes==1 and nodes_labels[0]=="SM": nodes_labels = ["EFT"] #Special case, want to have EFT as signal (not SM) by convention

    #-- For each NN output node, fill+write train/test histos for each class
    # print('nodes_labels', nodes_labels); print('list_labels', list_labels)
    for inode in range(nofOutputNodes):
        if (len(list_predictions_test_allNodes_allClasses)==1 or opts["strategy"] in ["CASCAL"]) and inode>0: break

        fout.cd()
        outname = 'hist_train_NODE_'+nodes_labels[inode]+'_allClasses'
        if store_trainHisto: hist_TrainingEvents_allClasses = TH1F(outname, '', nbins, 0, 1); hist_TrainingEvents_allClasses.Sumw2(); hist_TrainingEvents_allClasses.SetDirectory(0)
        outname = 'hist_test_NODE_'+nodes_labels[inode]+'_allClasses'
        hist_TestingEvents_allClasses = TH1F(outname, '', nbins, 0, 1); hist_TestingEvents_allClasses.Sumw2(); hist_TestingEvents_allClasses.SetDirectory(0)

        for iclass in range(len(list_labels)):
            # print('inode', inode, 'iclass', iclass)

            fout.cd()

            if store_trainHisto:
                hist_TrainingEvents_class = TH1F('hist_train_NODE_'+nodes_labels[inode]+'_CLASS_'+list_labels[iclass], '', nbins, 0, 1); hist_TrainingEvents_class.Sumw2(); hist_TrainingEvents_class.SetDirectory(0)
                fill_hist(hist_TrainingEvents_class, list_predictions_train_allNodes_allClasses[inode][iclass][:maxEvents], weights=list_PhysicalWeightsTrain_allClasses[iclass][:maxEvents])
                hist_TrainingEvents_class.Scale(1./hist_TrainingEvents_class.Integral())
                hist_TrainingEvents_class.Write()

            hist_TestingEvents_class = TH1F('hist_test_NODE_'+nodes_labels[inode]+'_CLASS_'+list_labels[iclass], '', nbins, 0, 1); hist_TestingEvents_class.Sumw2(); hist_TestingEvents_class.SetDirectory(0)
            fill_hist(hist_TestingEvents_class, list_predictions_test_allNodes_allClasses[inode][iclass][:maxEvents], weights=list_PhysicalWeightsTest_allClasses[iclass][:maxEvents])
            hist_TestingEvents_class.Scale(1./hist_TestingEvents_class.Integral())
            hist_TestingEvents_class.Write()

            #-- Sum all classes in a given node #Useless ?
            # if store_trainHisto: fill_hist(hist_TrainingEvents_allClasses, list_predictions_train_allNodes_allClasses[inode][iclass][:maxEvents], weights=list_PhysicalWeightsTrain_allClasses[iclass][:maxEvents])
            # fill_hist(hist_TestingEvents_allClasses, list_predictions_test_allNodes_allClasses[inode][iclass][:maxEvents], weights=list_PhysicalWeightsTest_allClasses[iclass][:maxEvents])

            # Also store histogram for selected input feature, for ROC comparison
            if idx_compare_ROC_inputFeature>=0 and inode==0:

                fout_feature.cd()
                if store_trainHisto:
                    hist_TrainingEvents_class_feature = TH1F('hist_train_NODE_'+nodes_labels[inode]+'_CLASS_'+list_labels[iclass], '', nbins, xmin_feature, xmax_feature); hist_TrainingEvents_class_feature.Sumw2(); hist_TrainingEvents_class_feature.SetDirectory(0)
                    fill_hist(hist_TrainingEvents_class_feature, list_xTrain_allClasses[iclass][:maxEvents,idx_compare_ROC_inputFeature], weights=list_PhysicalWeightsTrain_allClasses[iclass][:maxEvents])
                    hist_TrainingEvents_class_feature.Scale(1./hist_TrainingEvents_class_feature.Integral())
                    hist_TrainingEvents_class_feature.Write()

                hist_TestingEvents_class_feature = TH1F('hist_test_NODE_'+nodes_labels[inode]+'_CLASS_'+list_labels[iclass], '', nbins, xmin_feature, xmax_feature); hist_TestingEvents_class_feature.Sumw2(); hist_TestingEvents_class_feature.SetDirectory(0)
                fill_hist(hist_TestingEvents_class_feature, list_xTest_allClasses[iclass][:maxEvents,idx_compare_ROC_inputFeature], weights=list_PhysicalWeightsTest_allClasses[iclass][:maxEvents])
                hist_TestingEvents_class_feature.Scale(1./hist_TestingEvents_class_feature.Integral())
                hist_TestingEvents_class_feature.Write()

    fout.Close()
    print(colors.fg.lightgrey, "\nSaved output ROOT file containing Keras Predictions as histograms :", colors.reset, rootfile_outname, '\n')
    if idx_compare_ROC_inputFeature>=0:
        fout_feature.Close()
        print(colors.fg.lightgrey, "Saved output ROOT file containing Keras Predictions as histograms (for comparison with single feature) :", colors.reset, rootfile_outname_feature, '\n')

    return

# //--------------------------------------------
# //--------------------------------------------






# //--------------------------------------------
# //--------------------------------------------
##     ##    ###    ##       #### ########     ###    ######## ####  #######  ##    ##
##     ##   ## ##   ##        ##  ##     ##   ## ##      ##     ##  ##     ## ###   ##
##     ##  ##   ##  ##        ##  ##     ##  ##   ##     ##     ##  ##     ## ####  ##
##     ## ##     ## ##        ##  ##     ## ##     ##    ##     ##  ##     ## ## ## ##
 ##   ##  ######### ##        ##  ##     ## #########    ##     ##  ##     ## ##  ####
  ## ##   ##     ## ##        ##  ##     ## ##     ##    ##     ##  ##     ## ##   ###
   ###    ##     ## ######## #### ########  ##     ##    ##    ####  #######  ##    ##
# //--------------------------------------------
# //--------------------------------------------

def Make_Default_Validation_Plots(opts, list_features, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, PhysicalWeights_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, x, y_train, y_test, y_process, y_process_train, y_process_test, list_yTrain_allClasses, list_yTest_allClasses, list_xTrain_allClasses, list_xTest_allClasses, metrics, weight_dir, scores_allClasses_eachOperator, model, score=None, history=None):
    '''
    Call all the relevant sub-functions to create validation plots.
    '''

    #-- Useless, but anyway...
    plt.style.use('default')
    matplotlib.rc_file_defaults() #Restore matplotlib default settings

    print('\n'); print(colors.fg.lightblue, "--- Create control plots...", colors.reset); print('\n')

    if opts["testToy1D"]: Make_Test1D_Plot(opts, model)

    Control_Printouts(opts, list_labels, y_test, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTest_allClasses, score)

    Make_Loss_Plot(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, weight_dir, history)

    Make_Metrics_Plot(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, metrics, weight_dir, history)

    Make_ROC_plots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, list_xTrain_allClasses, list_xTest_allClasses, weight_dir)

    Make_Overtraining_plots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_xTrain_allClasses, list_xTest_allClasses, weight_dir)

    if opts["regress"] == True: Make_Regressor_ControlPlots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, y_test, y_process_test, list_yTest_allClasses, weight_dir, list_xTest_allClasses)

    Create_Correlation_Plots(opts, x, y_process, list_features, weight_dir)

    Plot_Input_Features(opts, x, y_process, PhysicalWeights_allClasses, list_features, weight_dir, False)

    if opts["strategy"] in ["regressor", "ROLR", "RASCAL"]:
        Plot_LR_Pred_vs_Truth(opts, list_features, list_labels, list_yTrain_allClasses, list_yTest_allClasses, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_truth_Test_allClasses, list_xTrain_allClasses, list_xTest_allClasses, weight_dir)
        Make_Pull_Plot(opts, weight_dir, list_yTest_allClasses, list_predictions_test_allNodes_allClasses, list_truth_Test_allClasses, list_PhysicalWeightsTest_allClasses, list_xTest_allClasses)
        doEvaluationPlots(list_yTest_allClasses[0], list_predictions_test_allNodes_allClasses[0][0], list_PhysicalWeightsTest_allClasses[0], weight_dir)

    if opts["shapPlots"]: Make_SHAP_Plots(opts, model, weight_dir, list_xTrain_allClasses, list_xTest_allClasses, list_features)

    if opts["strategy"] in ["RASCAL", "CASCAL"]:Test_Make_Score_Plot(weight_dir, scores_allClasses_eachOperator, y_process, x)

    return

# //--------------------------------------------
# //--------------------------------------------

 #####  #####  # #    # #####  ####  #    # #####
 #    # #    # # ##   #   #   #    # #    #   #
 #    # #    # # # #  #   #   #    # #    #   #
 #####  #####  # #  # #   #   #    # #    #   #
 #      #   #  # #   ##   #   #    # #    #   #
 #      #    # # #    #   #    ####   ####    #

def Control_Printouts(opts, list_labels, y_test, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTest_allClasses, score=None):
    '''
    Printout some information related to NN performance.
    '''

    if opts["nofOutputNodes"] == 1 and score is not None:
        loss = score[0]
        accuracy = score[1]
        print(colors.fg.lightgrey, '** Loss :', float('%.4g' % loss), colors.reset)
        print(colors.fg.lightgrey, '** Metrics :', float('%.4g' % accuracy), colors.reset)

    #-- Apply 2-sided KS test to train/test distributions, for first node / first process class
    KS_test(list_predictions_train_allNodes_allClasses[0][0], list_predictions_test_allNodes_allClasses[0][0])
    # KS_test(list_predictions_train_allNodes_allClasses[0][0][:,0], list_predictions_test_allNodes_allClasses[0][0][:,0])

    #-- Other tests : Anderson-Darling, Chi-2, ...
    # Anderson_Darling_test(list_predictions_train_allNodes_allClasses[0][0], list_predictions_train_allNodes_allClasses[0][0])
    # ChiSquare_test(list_predictions_train_allNodes_allClasses[0][0], list_predictions_train_allNodes_allClasses[0][0])

    #-- Classification report (problem with shapes)
    # if opts["nofOutputNodes"] == 1 and opts["strategy"] == "classifier":
    #     print('\n Classification report :')
    #     print(classification_report(y_true=y_test, y_pred=np.rint(np.concatenate(list_predictions_test_allNodes_allClasses[0])), sample_weight=np.concatenate(list_PhysicalWeightsTest_allClasses), labels=list_labels, target_names=list_labels) ) #Use rint to round to nearest integer (predict a class, not proba)
    #     print(confusion_matrix(y_true=y_test, y_pred=np.rint(np.concatenate(list_predictions_test_allNodes_allClasses[0])), sample_weight=np.concatenate(list_PhysicalWeightsTest_allClasses), labels=list_labels) )
    # print('\n')

    return

# //--------------------------------------------
# //--------------------------------------------

 #       ####   ####   ####
 #      #    # #      #
 #      #    #  ####   ####
 #      #    #      #      #
 #      #    # #    # #    #
 ######  ####   ####   ####

def Make_Loss_Plot(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, weight_dir, history=None):
    '''
    Plot the evolution of the loss versus epochs.

    NB : Possible solution to show 3 y-axes (train, test, lr) : https://stackoverflow.com/questions/48618992/matplotlib-graph-with-more-than-2-y-axes
    '''

    if history == None:
        print('history == None')
        return

    matplotlib.rc_file_defaults() #Restore matplotlib default settings

    # Plotting the loss with the number of iterations
    fig = plt.figure('loss')
    ax1 = fig.gca()
    timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
    timer.add_callback(close_event)

    plt.title('Loss VS Epoch')
    plt.ylabel('Training loss', color='darkorange')
    plt.xlabel('Epoch')
    lns1 = plt.plot(history.history['loss'], color='darkorange', label='Train') #Training loss

    ax2 = ax1.twinx()  # instantiate a second axis that shares the same x-axis
    ax2.set_ylabel('Test loss', color='cornflowerblue')  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', color='cornflowerblue')
    lns2 = ax2.plot(history.history['val_loss'], color='cornflowerblue', label='Test')
    #Validation loss starts usually at much larger values (because no regularization, different event weights, etc.) --> Modify the range of y-axis, zoom on converging part
    # ax2_ymin, ax2_ymax = ax2.get_ylim()
    # ax2.set_ylim(bottom=0, top=1)

    ax3 = ax1.twinx()  # instantiate a 3rd axis that shares the same x-axis
    lns3 = ax3.plot(history.history['lr'], color='dimgrey', linestyle='--', label='lr')
    ax3.get_yaxis().set_ticks([]) #invisible y axis (not important)

    #Legend
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc='best')
    # plt.legend(lns, labs, loc='upper right')

    if opts["displayImages"]==True:
        timer.start()
        plt.show()

    plotname = weight_dir + 'Loss_NN.png'
    fig.savefig(plotname, bbox_inches='tight') #bbox_inches='tight' ensures that second y-axis is visible
    print(colors.fg.lightgrey, "\nSaved Loss plot as :", colors.reset, plotname)
    fig.clear()
    plt.close('loss')

    return

# //--------------------------------------------
# //--------------------------------------------

 #    # ###### ##### #####  #  ####   ####
 ##  ## #        #   #    # # #    # #
 # ## # #####    #   #    # # #       ####
 #    # #        #   #####  # #           #
 #    # #        #   #   #  # #    # #    #
 #    # ######   #   #    # #  ####   ####

def Make_Metrics_Plot(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, metrics, weight_dir, history=None):
    '''
    Plot the evolution of the chosen metrics versus epochs.
    '''

    if history is None: return

    if opts["strategy"] == "RASCAL": #Metrics name depend on the output (e.g. 'val_likelihood_ratio_mean_squared_error' or 'val_score_mean_squared_error')
        metrics = "likelihood_ratio_" + str(metrics)

    # Plotting the error with the number of iterations
    fig = plt.figure('metric')
    ax1 = fig.gca()
    timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
    timer.add_callback(close_event)

    # plt.plot(history.history['val_'+metrics], color='cornflowerblue')
    plt.title('Accuracy VS Epoch')
    plt.ylabel('Train '+str(metrics), color='darkorange')
    plt.xlabel('Epoch')
    lns1 = plt.plot(history.history[metrics], color='darkorange', label='Train') #metrics name
    # print(history.history.keys())

    ax2 = ax1.twinx()  # instantiate a second axis that shares the same x-axis
    ax2.set_ylabel('Test '+str(metrics), color='cornflowerblue')  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', color='cornflowerblue')
    lns2 = ax2.plot(history.history['val_'+metrics], color='cornflowerblue', label='Test')

    ax3 = ax1.twinx()  # instantiate a 3rd axis that shares the same x-axis
    ax3.get_yaxis().set_ticks([]) #invisible y axis (not important)
    lns3 = ax3.plot(history.history['lr'], color='dimgrey', linestyle='--', label='lr')

    #Legend
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc='best')

    if opts["displayImages"]==True:
        timer.start()
        plt.show()

    plotname = weight_dir + 'Accuracy_NN.png'
    fig.savefig(plotname, bbox_inches='tight')
    # print("Saved Accuracy plot as : " + plotname)
    print(colors.fg.lightgrey, "\nSaved Accuracy plot as :", colors.reset, plotname)
    fig.clear()
    plt.close('metric')

    return

# //--------------------------------------------
# //--------------------------------------------

 #####   ####   ####   ####
 #    # #    # #    # #
 #    # #    # #       ####
 #####  #    # #           #
 #   #  #    # #    # #    #
 #    #  ####   ####   ####

def Make_ROC_plots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, list_xTrain_allClasses, list_xTest_allClasses, weight_dir):
    '''
    Plot ROC curves based on the separation between different processes/classes.

    See: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_curve.html
    /!\ NB: can not handle negative event weights ! Hence ROC curves may be non-physical...
    '''

    if opts["strategy"] is "regressor": return #No ROC to plot

    nofOutputNodes = opts["nofOutputNodes"]

    # Caveat: ROC curve can not handle negative weights. Either discard weights, or use absolute...
    list_PhysicalWeightsTest_allClasses_abs = np.absolute(list_PhysicalWeightsTest_allClasses)
    list_PhysicalWeightsTrain_allClasses_abs = np.absolute(list_PhysicalWeightsTrain_allClasses)

    lw = 2 #linewidth

    # print('MEAN: ', np.mean(list_predictions_train_allNodes_allClasses[0][0]))
    # print('MEAN: ', np.mean(list_predictions_test_allNodes_allClasses[0][0]))

    for inode in range(nofOutputNodes):

#-- Binary ROCs

        if opts["strategy"] in ["ROLR", "RASCAL"]: #Transformation r->s : s = 1/(r+1) #No sample weights (unweighted samples)

            if inode != 0: continue #Only for r node

            fpr, tpr, _ = roc_curve(np.concatenate(list_truth_Test_allClasses), 1./(np.concatenate(list_predictions_test_allNodes_allClasses[inode])+1))
            roc_auc = auc(fpr, tpr)

            fpr_train, tpr_train, _ = roc_curve(np.concatenate(list_truth_Train_allClasses), 1./(np.concatenate(list_predictions_train_allNodes_allClasses[inode])+1))
            roc_auc_train = auc(fpr_train, tpr_train)

        elif opts["nofOutputNodes"] == 1 or opts["strategy"] == "CASCAL": #Special case: CASCAL strategy has >1 outputs, but truth is binary (0/1)

            if opts["strategy"] == "CASCAL" and inode > 0: break #Only consider the classification node ?

            if opts["trainAtManyEFTpoints"] is False:
                fpr, tpr, _ = roc_curve(np.concatenate(list_truth_Test_allClasses), np.concatenate(list_predictions_test_allNodes_allClasses[inode]), sample_weight=np.concatenate(list_PhysicalWeightsTest_allClasses_abs) )
                fpr_train, tpr_train, _ = roc_curve(np.concatenate(list_truth_Train_allClasses), np.concatenate(list_predictions_train_allNodes_allClasses[inode]), sample_weight=np.concatenate(list_PhysicalWeightsTrain_allClasses_abs) )
            else:
                fpr, tpr, _ = roc_curve(np.concatenate(list_truth_Test_allClasses), np.concatenate(list_predictions_test_allNodes_allClasses[inode]))
                fpr_train, tpr_train, _ = roc_curve(np.concatenate(list_truth_Train_allClasses), np.concatenate(list_predictions_train_allNodes_allClasses[inode]))

            roc_auc = auc(fpr, tpr)
            roc_auc_train = auc(fpr_train, tpr_train)

        elif opts["nofOutputNodes"] > 1: #Multiclass: 1 ROC (signal vs All) per node

            fpr = dict(); tpr = dict(); roc_auc = dict(); fpr_train = dict(); tpr_train = dict(); roc_auc_train = dict()
            if opts["trainAtManyEFTpoints"] is False:
                fpr[inode], tpr[inode], _ = roc_curve(np.concatenate(list_truth_Test_allClasses)[:,inode], np.concatenate(list_predictions_test_allNodes_allClasses[inode]), sample_weight=np.concatenate(list_PhysicalWeightsTest_allClasses_abs) )
                fpr_train[inode], tpr_train[inode], _ = roc_curve(np.concatenate(list_truth_Train_allClasses)[:,inode], np.concatenate(list_predictions_train_allNodes_allClasses[inode]), sample_weight=np.concatenate(list_PhysicalWeightsTrain_allClasses_abs) )
            else:
                if opts["strategy"] is "CARL_multiclass" and inode > 0: #SM events are used for training each node w/ corresponding WC values; but for evaluation of EFT nodes, need to only consider SM events corresponding to the current node (conservation of proba, meaningless for other nodes)

                    indices_SMevents_test = np.where(list_xTest_allClasses[0][:,-len(opts["listOperatorsParam"])-1+inode]!=0); indices_SMevents_test = np.squeeze(indices_SMevents_test) #2D (current_node,indices) --> 1D (indices)
                    # print('np.shape(indices_SMevents_test)', np.shape(indices_SMevents_test))
                    truth_test_tmp = np.concatenate((list_truth_Test_allClasses[inode][:,inode], list_truth_Test_allClasses[0][indices_SMevents_test,inode]))
                    pred_test_tmp = np.concatenate((list_predictions_test_allNodes_allClasses[inode][inode], list_predictions_test_allNodes_allClasses[inode][0][indices_SMevents_test]))

                    indices_SMevents_train = np.where(list_xTrain_allClasses[0][:,-len(opts["listOperatorsParam"])-1+inode]!=0); indices_SMevents_train = np.squeeze(indices_SMevents_train)
                    truth_train_tmp = np.concatenate((list_truth_Train_allClasses[inode][:,inode], list_truth_Train_allClasses[0][indices_SMevents_train,inode]))
                    pred_train_tmp = np.concatenate((list_predictions_train_allNodes_allClasses[inode][inode], list_predictions_train_allNodes_allClasses[inode][0][indices_SMevents_train]))

                    fpr[inode], tpr[inode], _ = roc_curve(truth_test_tmp, pred_test_tmp)
                    fpr_train[inode], tpr_train[inode], _ = roc_curve(truth_train_tmp, pred_train_tmp)

                else:
                    fpr[inode], tpr[inode], _ = roc_curve(np.concatenate(list_truth_Test_allClasses)[:,inode], np.concatenate(list_predictions_test_allNodes_allClasses[inode]))
                    fpr_train[inode], tpr_train[inode], _ = roc_curve(np.concatenate(list_truth_Train_allClasses)[:,inode], np.concatenate(list_predictions_train_allNodes_allClasses[inode]))

            roc_auc[inode] = auc(fpr[inode], tpr[inode])
            roc_auc_train[inode] = auc(fpr_train[inode], tpr_train[inode])

        # Plot ROC curves
        fig = plt.figure()
        timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
        timer.add_callback(close_event)

        if opts["nofOutputNodes"] == 1 or opts["strategy"] in ["ROLR", "RASCAL", "CASCAL"]: #Consider single node
            if opts["makeValPlotsOnly"] is False: plt.plot(tpr_train, 1-fpr_train, color='darkorange', lw=lw, label='ROC NN (train) (AUC = {1:0.2f})' ''.format(0, roc_auc_train)) #If making validation plots only --> No 'training data'
            plt.plot(tpr, 1-fpr, color='cornflowerblue', lw=lw, label='ROC NN (test) (AUC = {1:0.2f})' ''.format(0, roc_auc))
        else: #for each node
            if opts["makeValPlotsOnly"] is False: plt.plot(tpr_train[inode], 1-fpr_train[inode], color='darkorange', lw=lw, label='ROC NN (train) (AUC = {1:0.2f})' ''.format(0, roc_auc_train[inode]))
            plt.plot(tpr[inode], 1-fpr[inode], color='cornflowerblue', lw=lw, label='ROC NN (test) (AUC = {1:0.2f})' ''.format(0, roc_auc[inode]))

        ax = fig.gca()
        ax.set_xticks(np.arange(0, 1, 0.1))
        ax.set_yticks(np.arange(0, 1., 0.1))
        plt.grid()
        plt.plot([1, 0], [0, 1], 'k--', lw=lw)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.0])
        plt.xlabel('Signal efficiency')
        plt.ylabel('Background rejection')
        plt.title('')
        plt.legend(loc='best')

        #Display plot in terminal for quick check (only for first node)
        if inode == 0 and opts["displayImages"]==True:
            timer.start()
            plt.show()

        plotname = weight_dir + 'ROC_NN_' + list_labels[inode] + '.png'
        fig.savefig(plotname)
        print(colors.fg.lightgrey, "\nSaved ROC plot as :", colors.reset, plotname)
        fig.clear()
        plt.close(fig)

#-- Multiclass ROCs (using test data only)
        if nofOutputNodes > 1 and (opts["strategy"] is "classifier" or (opts["strategy"] is "CARL_multiclass" and inode==0)):

            # Compute ROC curve and ROC area for each class
            fpr = dict(); tpr = dict(); roc_auc = dict(); fpr_train = dict(); tpr_train = dict(); roc_auc_train = dict()

            # Plot ROC curves
            fig = plt.figure('multiclass')
            for iclass in range(len(list_labels)): #iclass <-> iterate over all classes (may be signal or bkg); inode <-> current 'signal'

                if iclass == inode: #current 'signal' (<-> node) vs all
                    if opts["trainAtManyEFTpoints"] is False: fpr[iclass], tpr[iclass], _ = roc_curve(np.concatenate(list_truth_Test_allClasses)[:,inode], np.concatenate(list_predictions_test_allNodes_allClasses[inode]), sample_weight=np.concatenate(list_PhysicalWeightsTest_allClasses_abs) )
                    elif opts["strategy"] is "CARL_multiclass": continue #For this strategy, 'SM vs all' makes no sense (will never evaluate it on mixture of several operators, only single operators at a time)
                    else: fpr[iclass], tpr[iclass], _ = roc_curve(np.concatenate(list_truth_Test_allClasses)[:,inode], np.concatenate(list_predictions_test_allNodes_allClasses[inode]))

                    roc_auc[iclass] = auc(fpr[iclass], tpr[iclass])
                    mylabel = 'vs All (AUC = {1:0.2f})' ''.format(0, roc_auc[iclass])

                else: #current 'signal' (<-> node) vs specific process
                    if opts["trainAtManyEFTpoints"] is False:
                        fpr[iclass], tpr[iclass], _ = roc_curve(np.concatenate((list_truth_Test_allClasses[inode],list_truth_Test_allClasses[iclass]))[:,inode], np.concatenate((list_predictions_test_allNodes_allClasses[inode][inode],list_predictions_test_allNodes_allClasses[inode][iclass])), sample_weight=np.concatenate((list_PhysicalWeightsTest_allClasses_abs[inode],list_PhysicalWeightsTest_allClasses_abs[iclass])) )
                    elif opts["strategy"] is "CARL_multiclass": #SM events are used for training each node w/ corresponding WC values; but for evaluation of EFT nodes, need to only consider SM events corresponding to the current node (conservation of proba, meaningless for other nodes)
                        indices_SMevents_test = np.where(list_xTest_allClasses[0][:,-len(opts["listOperatorsParam"])-1+inode]!=0); indices_SMevents_test = np.squeeze(indices_SMevents_test) #2D (current_node,indices) --> 1D (indices)
                        truth_test_tmp = np.concatenate((list_truth_Test_allClasses[iclass][:,inode], list_truth_Test_allClasses[0][indices_SMevents_test,inode]))
                        pred_test_tmp = np.concatenate((list_predictions_test_allNodes_allClasses[inode][iclass], list_predictions_test_allNodes_allClasses[inode][0][indices_SMevents_test]))
                        fpr[iclass], tpr[iclass], _ = roc_curve(truth_test_tmp, pred_test_tmp)
                    else:
                        fpr[iclass], tpr[iclass], _ = roc_curve(np.concatenate((list_truth_Test_allClasses[inode],list_truth_Test_allClasses[iclass]))[:,inode], np.concatenate((list_predictions_test_allNodes_allClasses[inode][inode],list_predictions_test_allNodes_allClasses[inode][iclass])) )
                    roc_auc[iclass] = auc(fpr[iclass], tpr[iclass])
                    mylabel = 'vs ' + list_labels[iclass]+' (AUC = {1:0.2f})' ''.format(0, roc_auc[iclass])

                plt.plot(tpr[iclass], 1-fpr[iclass], lw=lw, label=mylabel)

            ax = fig.gca()
            ax.set_xticks(np.arange(0, 1, 0.1))
            ax.set_yticks(np.arange(0, 1., 0.1))
            plt.grid()
            plt.plot([1, 0], [0, 1], 'k--', lw=lw)
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.0])
            plt.xlabel('Signal efficiency')
            plt.ylabel('Background rejection')
            plt.title('ROC curves for '+list_labels[inode]+ ' (test data)')
            plt.legend(loc='best')

            plotname = weight_dir + 'ROC_NN_'+list_labels[inode]+'_allClasses.png'
            fig.savefig(plotname)
            print(colors.fg.lightgrey, "\nSaved ROC plot as :", colors.reset, plotname)
            fig.clear()
            plt.close('multiclass')

    return

# //--------------------------------------------
# //--------------------------------------------

  ####  #    # ###### #####  ##### #####    ##   # #    #
 #    # #    # #      #    #   #   #    #  #  #  # ##   #
 #    # #    # #####  #    #   #   #    # #    # # # #  #
 #    # #    # #      #####    #   #####  ###### # #  # #
 #    #  #  #  #      #   #    #   #   #  #    # # #   ##
  ####    ##   ###### #    #   #   #    # #    # # #    #

def Make_Overtraining_plots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_xTrain_allClasses, list_xTest_allClasses, weight_dir):
    '''
    Superimpose the distributions of the NN output for 'signal' and 'backgrounds'.

    NB: follow TMVA conventions: red=bkg, blue=sig; dots=train, filled/hatched=test
    '''

    nofOutputNodes = opts["nofOutputNodes"]

    use_simplified_plotLayout = False #True <-> simplify the plot layout (no test data, ...), e.g. to make illustrative xkcd-style plots -- hard-coded !

    for inode in range(nofOutputNodes):

        for iter in [0,1]: #Dummy loop: 0 <-> plot node vs all process classes; 1 <-> plot node vs each process class

            if iter > 0 and (nofOutputNodes == 1 or opts["strategy"] != "classifier"): break #Only makes sense for multiclassifier
            if opts["strategy"] in ["ROLR", "RASCAL", "CASCAL"] and inode > 0: break #Only for r node

            nbins = 30
            rmin = 0.; rmax = 1.

            fig = plt.figure('overtrain') #figsize=(30,15), dpi=200
            timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
            timer.add_callback(close_event)

            ax = plt.axes()
            ax.set_xlim([rmin,rmax])
            # plt.gcf().subplots_adjust(bottom=0.10)

            #--- COSMETICS
            #grey=#E6E6E6 #white=#FFFFFF
            ax.patch.set_edgecolor('black')
            ax.patch.set_facecolor('#E6E6E6') #inner bkg color
            # ax.patch.set_facecolor('white')
            plt.grid(color='w', linestyle='solid') # draw solid white grid lines
            # hide axis spines
            for spine in ax.spines.values():
                spine.set_visible(False)

                # hide top and right ticks
                ax.xaxis.tick_bottom()
                # ax.yaxis.tick_left()

                # lighten ticks and labels
                ax.tick_params(colors='gray', direction='out')
                for tick in ax.get_xticklabels():
                    tick.set_color('gray')
                    # tick.set_weight('bold')
                    for tick in ax.get_yticklabels():
                        tick.set_color('gray')

            if iter > 0: #Plot current node vs each other process class

                cmap = [plt.cm.Set1(i) for i in range(20)] #Qualitative colormap

                for iclass in range(len(list_labels)):

                    col = cmap[iclass] #Same color for train/test
                    plt.hist(list_predictions_train_allNodes_allClasses[inode][iclass], bins=nbins, range=(rmin,rmax), weights=list_PhysicalWeightsTrain_allClasses[iclass], density=True, histtype='step', log=False, label=list_labels[iclass]+' (Train)', color=col) #Training data for current class --> Line
                    plt.hist(list_predictions_test_allNodes_allClasses[inode][iclass], bins=nbins, range=(rmin,rmax), weights=list_PhysicalWeightsTest_allClasses[iclass], density=True, histtype='step', log=False, label=list_labels[iclass]+' (Test)', color=col, fill=True, alpha=0.30) #Testing data for current class --> Filled

            else: #Plot current node vs all other process classes

                #Class labels (for legend)
                label_class0 = "Signal"; label_class1 = "Backgrounds"
                if opts["trainAtManyEFTpoints"] == True and inode == 0: label_class0 = "EFT"; label_class1 = "SM" #Only have 'SM vs EFT' scenario when considering SM node vs the rest (EFT)
                # if opts["trainAtManyEFTpoints"] == True and inode == 0: label_class0 = "SM"; label_class1 = "EFT" #Only have 'SM vs EFT' scenario when considering SM node vs the rest (EFT)

                #-- Plot TEST sig/bkg histos, normalized (no errors displayed <-> don't need TH1Fs)
                suffix = " (Test)"
                if use_simplified_plotLayout: suffix = ""

                if opts["strategy"] in ["ROLR", "RASCAL"]: tmp = 1/(list_predictions_test_allNodes_allClasses[inode][inode]+1) #Transform r -> s
                else: tmp = list_predictions_test_allNodes_allClasses[inode][inode]

                plt.hist(tmp, bins=nbins, range=(rmin,rmax), color= 'cornflowerblue', alpha=0.50, weights=list_PhysicalWeightsTest_allClasses[inode], density=True, histtype='step', log=False, label=list_labels[inode]+' (Test)', edgecolor='cornflowerblue',fill=True) #Testing data for current class --> Filled blue

                #-- Trick : want to get arrays of predictions/weights for all test events *which do not belong to class of current node* => Loop on classes, check if matches node
                lists_predictions_bkgs = []; lists_weights_bkg = []
                for iclass in range(0, len(list_labels)):
                    if iclass != inode:
                        # print(list_labels[iclass], list_PhysicalWeightsTest_allClasses[iclass][:20])
                        if opts["strategy"] is "CARL_multiclass" and inode>0:
                            if iclass != 0: continue #For EFT nodes, only consider SM as background, not the other EFT operators
                            indices_SMevents = np.where(list_xTest_allClasses[0][:,-len(opts["listOperatorsParam"])-1+inode]!=0) #SM events are used for training each node w/ corresponding WC values; but for evaluation of EFT nodes, need to only consider SM events corresponding to the current node (conservation of proba, meaningless for other nodes)
                            pred_tmp = list_predictions_test_allNodes_allClasses[inode][0][indices_SMevents]
                            w_tmp = list_PhysicalWeightsTest_allClasses[0][indices_SMevents]
                        elif opts["strategy"] in ["ROLR", "RASCAL"]:
                            pred_tmp = 1/(list_predictions_test_allNodes_allClasses[inode][iclass]+1) #Transform r -> s
                            w_tmp = list_PhysicalWeightsTest_allClasses[iclass]
                        else:
                            pred_tmp = list_predictions_test_allNodes_allClasses[inode][iclass]
                            w_tmp = list_PhysicalWeightsTest_allClasses[iclass]
                        lists_predictions_bkgs.append(pred_tmp)
                        lists_weights_bkg.append(w_tmp)

                #Plot test bkg --> Red hatched
                if len(lists_predictions_bkgs) > 0: #Protection
                    predictions_bkgs = np.concatenate(lists_predictions_bkgs); weights_bkgs = np.concatenate(lists_weights_bkg)
                    if use_simplified_plotLayout: plt.hist(predictions_bkgs, bins=nbins, range=(rmin,rmax), color='orangered', alpha=0.50, weights=weights_bkgs, density=True, histtype='step', log=False, label=label_class1+suffix, edgecolor='orangered',fill=False)
                    else: plt.hist(predictions_bkgs, bins=nbins, range=(rmin,rmax), color='orangered', alpha=0.50, weights=weights_bkgs, density=True, histtype='step', log=False, label=label_class1+suffix, hatch='/', edgecolor='orangered',fill=False)

                #-- For TRAIN histos, we want to compute/display the bin errors correctly; to do this we first fill TH1Fs, then read their bin contents/errors
                hist_overtrain_sig = TH1F('hist_overtrain_sig', '', nbins, rmin, rmax); hist_overtrain_sig.Sumw2(); hist_overtrain_sig.SetDirectory(0)
                hist_overtrain_bkg = TH1F('hist_overtrain_bkg', '', nbins, rmin, rmax); hist_overtrain_bkg.Sumw2(); hist_overtrain_bkg.SetDirectory(0)

                #-- Only signal process
                if opts["strategy"] is "regressor" and opts["nofOutputNodes"] > len(list_labels): continue
                for val, w in zip((list_predictions_train_allNodes_allClasses[inode][inode]), list_PhysicalWeightsTrain_allClasses[inode]): #'zip' stops when the shorter of the lists stops
                    # if opts["strategy"] in ["ROLR", "RASCAL"]: val = 1/(val+1) #Transform r -> s
                    if opts["strategy"] in ["ROLR", "RASCAL"]: val = s_from_r(val) #Transform r -> s
                    hist_overtrain_sig.Fill(val, w)

                #-- Only background processes
                for iclass in range(0, len(list_labels)):
                    if iclass != inode:
                        if opts["strategy"] is "CARL_multiclass" and inode > 0 and iclass != 0: continue #For EFT nodes, only consider SM as background, not the other EFT operators

                        for val, w, x in zip(list_predictions_train_allNodes_allClasses[inode][iclass], list_PhysicalWeightsTrain_allClasses[iclass], list_xTrain_allClasses[iclass]):
                            if opts["strategy"] is "CARL_multiclass" and inode > 0 and x[-len(opts["listOperatorsParam"])-1+inode]==0: continue #SM events are used for training each node w/ corresponding WC values; but for evaluation of EFT nodes, need to only consider SM events corresponding to the current node (conservation of proba, meaningless for other nodes)
                            if opts["strategy"] in ["ROLR", "RASCAL"]: val = s_from_r(val) #Transform r -> s
                            hist_overtrain_bkg.Fill(val, w)
                            # print(inode, iclass, val)

                #-- Normalize
                int_sig = hist_overtrain_sig.Integral(0,hist_overtrain_sig.GetNbinsX()+1)
                int_bkg = hist_overtrain_bkg.Integral(0,hist_overtrain_bkg.GetNbinsX()+1)
                if int_sig <= 0: int_sig = 1
                if int_bkg <= 0: int_bkg = 1
                sf_integral = abs(rmax - rmin) / nbins #h.Scale(1/integral) makes the sum of contents equal to 1, but does not account for the bin width
                hist_overtrain_sig.Scale(1./(int_sig*sf_integral))
                hist_overtrain_bkg.Scale(1./(int_bkg*sf_integral))

                #-- Read bin contents/errors
                bin_centres = []; counts_sig = []; err_sig = []; counts_bkg = []; err_bkg = []
                for ibin in range(1, hist_overtrain_sig.GetNbinsX()+1):
                    bin_centres.append(hist_overtrain_sig.GetBinCenter(ibin))
                    counts_sig.append(hist_overtrain_sig.GetBinContent(ibin)); counts_bkg.append(hist_overtrain_bkg.GetBinContent(ibin))
                    err_sig.append(hist_overtrain_sig.GetBinError(ibin)); err_bkg.append(hist_overtrain_bkg.GetBinError(ibin))

                #-- Plot TRAIN sig/bkg histos, normalized, with errorbars
                if use_simplified_plotLayout is False:
                    plt.errorbar(bin_centres, counts_sig, marker='o', yerr=err_sig, linestyle='None', markersize=6, color='blue', alpha=0.90, label=label_class0+' (Train)')
                    plt.errorbar(bin_centres, counts_bkg, marker='o', yerr=err_bkg, linestyle='None', markersize=6, color='red', alpha=0.90, label=label_class1+' (Train)')

            #-- Labels
            plt.legend(loc='best', numpoints=1)
            plt.title("Output distributions for Signal & Background")
            plt.grid(axis='y', alpha=0.75)
            plt.grid(axis='x', alpha=0.75)
            plt.xlabel("Classifier output")
            plt.ylabel('PDF')

            if use_simplified_plotLayout:
                plt.title("")
                plt.xlabel("")
                # plt.xlabel("Variable from Machine-Learning")

            if inode == 0 and opts["displayImages"]==True:
                timer.start()
                plt.show()

            suffix = ''
            if iter == 1: suffix = '_eachClass'
            plotname = weight_dir + 'Overtraining_NN_' + list_labels[inode] + suffix + '.png'
            fig.savefig(plotname, bbox_inches = "tight")
            # fig.savefig(plotname)
            # print("Saved Overtraining plot as : " + plotname)
            print(colors.fg.lightgrey, "\nSaved Overtraining plot as :", colors.reset, plotname)
            fig.clear()
            plt.close('overtrain')

    return

# //--------------------------------------------
# //--------------------------------------------

 #####  ######  ####  #####  ######  ####   ####   ####  #####
 #    # #      #    # #    # #      #      #      #    # #    #
 #    # #####  #      #    # #####   ####   ####  #    # #    #
 #####  #      #  ### #####  #           #      # #    # #####
 #   #  #      #    # #   #  #      #    # #    # #    # #   #
 #    # ######  ####  #    # ######  ####   ####   ####  #    #

def Make_Regressor_ControlPlots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, y_test, y_process_test, list_yTest_allClasses, weight_dir, list_xTest_allClasses):
    '''
    Create control plots specific to regressors: truth/reco correlation, etc.
    '''

    plot_truth_histos = True #True <-> Display truth histos on plot
    norm = False #True <-> Normalize histogram to unity

    nbins = 20
    rmin = -1.; rmax = 5
    nofOutputNodes = opts["nofOutputNodes"]

    label_class0 = "EFT"; label_class1 = "SM"
    if opts["strategy"] is "regressor": label_class0 = "Signal"; label_class1 = "Backgrounds"

    plotSingleEventClass = False #By default, plot sig+bkg
    if len(list_predictions_train_allNodes_allClasses[0]) == 1 or len(list_predictions_train_allNodes_allClasses[0][1]) == 0 or len(list_predictions_test_allNodes_allClasses[0][1]) == 0: plotSingleEventClass = True #e.g. for regressor on signle class --> no 'bkg class' to plot

    for inode in range(nofOutputNodes):

        #Auto-adjust plot range
        # '''
        for z in range(nofOutputNodes):
            tmp = np.concatenate(list_predictions_test_allNodes_allClasses[z])
            minElement = np.amin(tmp)
            rmin = minElement if minElement < rmin else rmin
            # maxElement = np.amax(tmp)
            # rmax = maxElement if maxElement > rmax else rmax
            quantile = np.quantile(tmp, 0.95)
            rmax = quantile if quantile > rmax else rmax #Define rmax so as to display ~95% of events (better in case few predictions are very off, axis range will stay reasonable)
        rmin = math.floor(rmin) #Round down
        rmax = math.ceil(rmax) #Round up
        # '''

        if opts["strategy"] in ["ROLR", "RASCAL"]: rmin=-2; rmax=2

        fig = plt.figure('regressor')
        # timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
        # timer.add_callback(close_event)

        ax = plt.axes()
        ax.set_xlim([rmin,rmax])

        #--- COSMETICS
        #grey=#E6E6E6 #white=#FFFFFF
        ax.patch.set_edgecolor('black')
        ax.patch.set_facecolor('#E6E6E6') #inner bkg color
        # ax.patch.set_facecolor('white')
        plt.grid(color='w', linestyle='solid') # draw solid white grid lines
        # hide axis spines
        for spine in ax.spines.values():
            spine.set_visible(False)

            # hide top and right ticks
            ax.xaxis.tick_bottom()
            # ax.yaxis.tick_left()

            # lighten ticks and labels
            ax.tick_params(colors='gray', direction='out')
            for tick in ax.get_xticklabels():
                tick.set_color('gray')
                for tick in ax.get_yticklabels():
                    tick.set_color('gray')

        """
        #-- Trick : for truth histos, we want to compute the bin errors correctly ; to do this we first fill TH1Fs, then read their bin contents/errors
        hist_truth_SM = TH1F('hist_truth_SM', '', nbins, rmin, rmax); hist_truth_SM.Sumw2(); hist_truth_SM.SetDirectory(0)
        hist_truth_EFT = TH1F('hist_truth_EFT', '', nbins, rmin, rmax); hist_truth_EFT.Sumw2(); hist_truth_EFT.SetDirectory(0)

        #Only 'SM' events
        for val in y_test[y_process_test==1]: #'zip' stops when the shorter of the lists stops
            if opts["strategy"] is "RASCAL": val = val[inode]
            hist_truth_SM.Fill(val, 1.)

        #Only 'EFT' events
        for val in y_test[y_process_test==0]: #'zip' stops when the shorter of the lists stops
            if opts["strategy"] is "RASCAL": val = val[inode]
            hist_truth_EFT.Fill(val, 1.)

        #Normalize
        int_sm = hist_truth_SM.Integral(0,hist_truth_SM.GetNbinsX()+1)
        int_eft = hist_truth_EFT.Integral(0,hist_truth_EFT.GetNbinsX()+1)
        if int_sm <= 0: int_sm = 1
        if int_eft <= 0: int_eft = 1
        sf_integral = abs(rmax - rmin) / nbins #h.Scale(1/integral) makes the sum of contents equal to 1, but does not account for the bin width
        hist_truth_SM.Scale(1./(int_sm*sf_integral))
        hist_truth_EFT.Scale(1./(int_eft*sf_integral))

        #Read bin contents/errors
        bin_centres = []; counts_SM = []; err_SM = []; counts_EFT = []; err_EFT = []
        for ibin in range(1, hist_truth_SM.GetNbinsX()+1):
            bin_centres.append(hist_truth_SM.GetBinCenter(ibin))
            counts_SM.append(hist_truth_SM.GetBinContent(ibin)); counts_EFT.append(hist_truth_EFT.GetBinContent(ibin))
            err_SM.append(hist_truth_SM.GetBinError(ibin)); err_EFT.append(hist_truth_EFT.GetBinError(ibin))

        counts_SM[0]+= hist_truth_SM.GetBinContent(0); counts_EFT[0]+= hist_truth_EFT.GetBinContent(0) # Add underflow
        counts_SM[-1]+= hist_truth_SM.GetBinContent(hist_truth_SM.GetNbinsX()+1); counts_EFT[-1]+= hist_truth_EFT.GetBinContent(hist_truth_EFT.GetNbinsX()+1) # Add overflow

        #Plot truth histos, normalized, with errorbars
        if plot_truth_histos == True:
            plt.errorbar(bin_centres, counts_SM, marker='o', yerr=err_SM, linestyle='None', markersize=6, color='blue', alpha=0.90, label=label_class0+' (Truth)')
            plt.errorbar(bin_centres, counts_EFT, marker='o', yerr=err_EFT, linestyle='None', markersize=6, color='red', alpha=0.90, label=label_class1+' (Truth)')
        """

        #Prediction and truth histos, normalized
        #NB: apply weights, but they are all equal to 1 for parameterized strategies
        if plot_truth_histos == True:
            if nofOutputNodes == 1:
                # plt.hist(y_test[y_process_test==1], weights=list_PhysicalWeightsTest_allClasses[0], bins=nbins, range=(rmin,rmax), color='cornflowerblue', density=norm, alpha=0.50, histtype='step', log=False, label=label_class0+" (Truth)", edgecolor='cornflowerblue',fill=False, linewidth=4.)
                # plt.hist(y_test[y_process_test==0], weights=list_PhysicalWeightsTest_allClasses[1], bins=nbins, range=(rmin,rmax), color='orangered', density=norm, alpha=0.50, histtype='step', log=False, label=label_class1+" (Truth)", edgecolor='orangered',fill=False, linewidth=4.)
                plt.hist(list_yTest_allClasses[0], weights=list_PhysicalWeightsTest_allClasses[0], bins=nbins, range=(rmin,rmax), color='cornflowerblue', density=norm, alpha=0.50, histtype='step', log=False, label=label_class0+" (Truth)", edgecolor='cornflowerblue',fill=False, linewidth=3.)
                if not plotSingleEventClass: plt.hist(list_yTest_allClasses[1], weights=list_PhysicalWeightsTest_allClasses[1], bins=nbins, range=(rmin,rmax), color='orangered', density=norm, alpha=0.50, histtype='step', log=False, label=label_class1+" (Truth)", edgecolor='orangered',fill=False, linewidth=3.)

                if opts["comparVarIdx"] >= 0:
                    plt.hist(list_xTest_allClasses[0][:,opts["comparVarIdx"]], weights=list_PhysicalWeightsTest_allClasses[0], bins=nbins, range=(rmin,rmax), density=norm, alpha=0.50, histtype='step', log=False, label=label_class0+" (Comparison)", edgecolor='orangered',fill=False, linewidth=3.)
                    # print('truth', list_yTest_allClasses[0][:10])
                    # print('pred', list_predictions_test_allNodes_allClasses[inode][0][:10])
                    # print('comp', list_xTest_allClasses[0][:10,opts["comparVarIdx"]])

            else:
                # plt.hist(y_test[y_process_test==1][:,inode], weights=list_PhysicalWeightsTest_allClasses[0], bins=nbins, range=(rmin,rmax), color='cornflowerblue', density=norm, alpha=0.50, histtype='step', log=False, label=label_class0+" (Truth)", edgecolor='cornflowerblue',fill=False, linewidth=4.)
                # plt.hist(y_test[y_process_test==0][:,inode], weights=list_PhysicalWeightsTest_allClasses[1], bins=nbins, range=(rmin,rmax), color='orangered', density=norm, alpha=0.50, histtype='step', log=False, label=label_class1+" (Truth)", edgecolor='orangered',fill=False, linewidth=4.)
                plt.hist(list_yTest_allClasses[0][:,inode], weights=list_PhysicalWeightsTest_allClasses[0], bins=nbins, range=(rmin,rmax), color='cornflowerblue', density=norm, alpha=0.50, histtype='step', log=False, label=label_class0+" (Truth)", edgecolor='cornflowerblue',fill=False, linewidth=3.)
                if not plotSingleEventClass: plt.hist(list_yTest_allClasses[1][:,inode], weights=list_PhysicalWeightsTest_allClasses[1], bins=nbins, range=(rmin,rmax), color='orangered', density=norm, alpha=0.50, histtype='step', log=False, label=label_class1+" (Truth)", edgecolor='orangered',fill=False, linewidth=3.)

        plt.hist(list_predictions_test_allNodes_allClasses[inode][0], weights=list_PhysicalWeightsTest_allClasses[0], bins=nbins, range=(rmin,rmax), color='cornflowerblue', density=norm, histtype='step', log=False, label=label_class0+" (Pred.)", alpha=0.50, edgecolor='cornflowerblue',fill=True)
        if not plotSingleEventClass: plt.hist(list_predictions_test_allNodes_allClasses[inode][1], weights=list_PhysicalWeightsTest_allClasses[1], bins=nbins, range=(rmin,rmax), color='orangered', density=norm, histtype='step', log=False, label=label_class1+" (Pred.)", alpha=0.50, edgecolor='orangered',fill=True)

        nodename = 'node' + str(inode)
        if opts["strategy"] in ["ROLR", "RASCAL"] and inode == 0: nodename = 'LR'
        elif opts["strategy"] is "RASCAL": nodename = 'Score_' + opts["listOperatorsParam"][inode-1] #First node is r

        plt.legend(loc='best', numpoints=1)
        plt.title("Regressor output for SM & EFT (test data)")
        plt.grid(axis='y', alpha=0.75)
        plt.grid(axis='x', alpha=0.75)
        plt.xlabel(nodename)
        plt.ylabel('PDF')

        plotname = weight_dir + 'Regressor_NN_' + nodename + '.png'
        fig.savefig(plotname)
        print(colors.fg.lightgrey, "\nSaved Regressor plot as :", colors.reset, plotname)
        fig.clear()
        plt.close('regressor')

    return

# //--------------------------------------------
# //--------------------------------------------

  ####   ####  #####  #####  ###### #        ##   ##### #  ####  #    #
 #    # #    # #    # #    # #      #       #  #    #   # #    # ##   #
 #      #    # #    # #    # #####  #      #    #   #   # #    # # #  #
 #      #    # #####  #####  #      #      ######   #   # #    # #  # #
 #    # #    # #   #  #   #  #      #      #    #   #   # #    # #   ##
  ####   ####  #    # #    # ###### ###### #    #   #   #  ####  #    #

def Create_Correlation_Plots(opts, x, y_process, list_features, weight_dir):
    '''
    Create the relevant input feature correlation plots, depending on strategy.
    '''

    if opts["strategy"] is "classifier": #Separate between 'sig' (first process) and bkg (other processes)
        if opts["nofOutputNodes"] == 1: x_sig = x[y_process==1]; x_bkg = x[y_process==0] #Binary
        else: x_sig = x[y_process[:,0]==1]; x_bkg = x[y_process[:,0]==0] #Multiclass
        Make_Correlation_Plot(opts, x_sig, list_features, weight_dir + 'CorrelMatrix_sig.png')
        Make_Correlation_Plot(opts, x_bkg, list_features, weight_dir + 'CorrelMatrix_bkg.png')

    elif opts["trainAtManyEFTpoints"] is True: #Separate between 'SM' and 'EFT'
        x_SM = x[y_process==0]; x_EFT = x[y_process>0]
        Make_Correlation_Plot(opts, x_SM, list_features, weight_dir + 'CorrelMatrix_SM.png')
        Make_Correlation_Plot(opts, x_EFT, list_features, weight_dir + 'CorrelMatrix_EFT.png')

    else: Make_Correlation_Plot(opts, x, list_features, weight_dir + 'CorrelMatrix.png') #Use all events

    return

def Make_Correlation_Plot(opts, x, list_features, outname):
    '''
    Make a correlation plot of input features.
    '''

# //--------------------------------------------
    doNotPlotP4 = False #Can choose to only consider high-level variables (not p4 variables) to improve readability
# //--------------------------------------------

    # sns.set(font_scale=1.4) #Scale up label font size #NB : also sets plotting options to seaborn's default

    #Can avoid plotting theory parameter inputs (by default), and p4 variables
    list_features = np.array(list_features)
    indices = []
    for ivar in range(len(list_features)):
        if doNotPlotP4 and len(list_features) > 5 and list_features[ivar].endswith(('_pt','_eta','_phi','_DeepCSV','_DeepJet')): indices.append(ivar) #Substrings corresponding to p4 vars (hardcoded)
        if opts["parameterizedNN"] == True and ivar >= len(list_features)-len(opts["listOperatorsParam"]): indices.append(ivar)

    mask = np.ones(len(list_features), np.bool) #Default
    if len(indices) < len(list_features): #Protection: ignore doNotPlotP4 if only low-level vars are included
        indices = np.array(indices, dtype=int)
        mask[indices] = 0
        # print(indices)

    maxEvents=50000 #Don't need all events

    #-- Convert np array to pd dataframe
    df = pd.DataFrame(data=x[:maxEvents,:][:,mask], columns=list_features[mask]) #x = (events, vars) ; colums names are var names
    # print(df); print(df.describe())

    #-- Get correlation matrix
    corr = df.corr()

    mask_diag = np.triu(np.ones_like(corr, dtype=np.bool)) #Mask upper right triangle
    # mask_diag = np.tril(np.ones_like(corr, dtype=np.bool)) #Mask bottom left triangle

    #-- Take abs(values)
    # corr = abs(corr)

    #Set small values to 0
    corr[np.abs(corr)<.01] = 0

    #-- Color palette
    # palette = sns.diverging_palette(240, 10, n=9)
    # palette = sns.diverging_palette(20, 220, n=256)
    palette = 'coolwarm'
    # palette = 'YlGn' #From yellow to green

    # fig, ax = plt.subplots()
    fig, ax = plt.subplots(figsize=(10, 10))
    # plt.title('Input features correlations', y=-0.01)
    ax.set_title('Input features correlations',fontsize=30)

    # Draw the heatmap -- see : https://seaborn.pydata.org/generated/seaborn.heatmap.html
    if len(list_features[mask]) > 15: #Lots of features to display, remove values
        hm = sns.heatmap(corr, mask=mask_diag, cmap=palette, vmin=-1., vmax=1., center=0, square=True, linewidths=0.5, annot = False, fmt='.1g', cbar_kws={"shrink": .5},)
        hm.set_xticklabels(hm.get_xticklabels(), horizontalalignment='right', fontsize = 14) #bottom labels
        hm.set_yticklabels(hm.get_yticklabels(), fontsize = 14)
    else:
        hm = sns.heatmap(corr, mask=mask_diag, cmap=palette, vmin=-1., vmax=1., center=0, square=True, linewidths=0.5, annot = True, fmt='.1g', cbar_kws={"shrink": .5},)
        hm.set_xticklabels(hm.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize = 14) #bottom labels
        hm.set_yticklabels(hm.get_yticklabels(), fontsize = 14)

    # ax.set_ylim(bottom=-0.5, top=len(list_features)+0.5)
    # ax.set_xlim(left=-0.5, right=len(list_features)+0.5)
    # ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False) # x labels appear on top
    # hm.set_xticklabels(hm.get_xticklabels(), rotation=45, horizontalalignment='left', va='bottom', fontsize = 16)
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True) # x labels appear at bottom
    fig.tight_layout()

    # plotname = weight_dir + 'CorrelMatrix.png'
    plt.savefig(outname)
    print(colors.fg.lightgrey, "\nSaved correlation matrix plot as :", colors.reset, outname)
    plt.close('all')

    return

# //--------------------------------------------
# //--------------------------------------------

 ###### ######   ##   ##### #    # #####  ######  ####
 #      #       #  #    #   #    # #    # #      #
 #####  #####  #    #   #   #    # #    # #####   ####
 #      #      ######   #   #    # #####  #           #
 #      #      #    #   #   #    # #   #  #      #    #
 #      ###### #    #   #    ####  #    # ######  ####

def Plot_Input_Features(opts, x, y_process, weights, list_features, weight_dir, isControlNorm=False):
    '''
    Plot distributions of input features (train/etst/both depending on dataset given in arg).

    NB: 'isControlNorm' <-> rescaled features, training data (for control plots)
    '''

    useMostExtremeWCvaluesOnly = True #True <-> for 'EFT' class, will only consider points generated at the most extreme WC values included during training (not all the intermediate points) #Use this to create more "representative" val plots, in which only a few specific WC values are included instead of all points
    plot_eachSingleFeature = False #True <-> also save 1 single plot per feature
    doNotPlotP4 = True #Can choose to only consider high-level variables (not p4 variables) to improve readability
    nbins = 30 #Binning for plots

# //--------------------------------------------

    if opts["makeValPlotsOnly"] is True and isControlNorm is True: return
    if x is None: print('Error, can\'t produce input features plots : x=None !'); return

    if isControlNorm: useMostExtremeWCvaluesOnly = False #Can't have both options at the same time, because 'isControlNorm' corresponds to transformed inputs, so we can not condition upon the values of some features (the input WC values) to select events !

    #NB: don't call plt functions here (mess with layout of subsequent plots...)
    list_features = np.array(list_features)
    sns.set(palette='coolwarm', font_scale=1.4) #Scale up label font size #NB : this also sets plotting options to seaborn's default

    #-- First, shuffle the input feature array (necessary e.g. for param. NN, since x array is ordered by theta values)
    shuf = np.arange(len(x))
    np.random.shuffle(shuf)
    x = np.copy(x[shuf]); y_process = np.copy(y_process[shuf]); weights = np.copy(weights[shuf])

    nMax = 500000 # >0 <-> Don't use more events (slow) #Warning : for parameterized NN, biases distributions of WCs (will most likely not cover all EFT points)
    if nMax>0 and len(x) > nMax: x = x[:nMax]; y_process = y_process[:nMax]; weights = weights[:nMax] #Else, too slow

    if opts["parameterizedNN"] == True: #Only retain EFT events whose operator values are at boundaries #NB: for SM can select all events (input features don't depend on WC, only output prediction)
        if useMostExtremeWCvaluesOnly is True:
            sl = np.s_[x.shape[1]-len(opts["listOperatorsParam"]):] #Specify slice corresponding to indices of theory parameters features (placed last)
            # indices_class0 = np.where(y_process==0) #Class 0
            indices_class0 = np.where(np.logical_and(y_process==0, np.logical_or(np.all(x[:,sl]==opts['minWC'],axis=1),np.all(x[:,sl]==opts['maxWC'],axis=1))) ) #Class 1
            indices_class1 = np.where(np.logical_and(y_process==1, np.logical_or(np.all(x[:,sl]==opts['minWC'],axis=1),np.all(x[:,sl]==opts['maxWC'],axis=1))) ) #Class 1
            # print(x[indices_class0].shape); print(x[indices_class1].shape)
            # print(np.mean(x[indices_class0,0].T, axis=0)); print(np.mean(x[indices_class1,0].T, axis=0))
            # print(x[-10:])
            x = np.append(x[indices_class0],x[indices_class1],axis=0); y_process = np.append(y_process[indices_class0],y_process[indices_class1],axis=0); weights = np.append(weights[indices_class0],weights[indices_class1],axis=0)

    #Can avoid plotting theory parameter inputs (by default), and p4 variables
    list_features = np.array(list_features)
    indices = []
    for ivar in range(len(list_features)):
        if doNotPlotP4 and len(list_features) > 5 and list_features[ivar].endswith(('_pt','_eta','_phi','_DeepCSV','_DeepJet')): indices.append(ivar) #Substrings corresponding to p4 vars (hardcoded)
        # if opts["parameterizedNN"] == True and ivar >= len(list_features)-len(opts["listOperatorsParam"]): indices.append(ivar) #Don't show input WCs

    mask = np.ones(len(list_features), np.bool) #Default
    if len(indices) < len(list_features): #Protection: ignore doNotPlotP4 if only low-level vars are included
        indices = np.array(indices, dtype=int)
        mask[indices] = 0

    #-- Convert np array to pd dataframe
    df = pd.DataFrame(data=x[:,:][:,mask], columns=list_features[mask]) #x = (events, vars) ; colums names are var names

    #-- Can tweak panda options to print full df instead of partial
    # pd.set_option('display.max_rows', None)
    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.width', None)
    # pd.set_option('display.max_colwidth', -1)
    # print(df['maxDelRbL'][:])

    #-- Insert a column corresponding to the class label
    if y_process.ndim==1: df.insert(loc=0, column='class', value=y_process[:len(x)], allow_duplicates=False)
    else: df.insert(loc=0, column='class', value=y_process[:len(x),0], allow_duplicates=False) #Only care about first column=main signal (rest -> bkg)

    #-- Insert a column corresponding to physical event weights
    df.insert(loc=0, column='weight', value=weights[:len(x),]) #Only care about first column=main signal (rest -> bkg)
    # print(df); print(df.describe())

    #-- Create multiplot #NB: only process columns corresponding to phy vars
    df[df['class']==1].hist(figsize=(30,30), label='Signal', column=list_features[mask], weights=df['weight'][df['class']==1], bins=nbins, alpha=0.4, density=True, color='b') #signal
    df[df['class']==0].hist(figsize=(30,30), label='Backgrounds', column=list_features[mask], weights=df['weight'][df['class']==0], bins=nbins, alpha=0.4, density=True, color='r', ax=plt.gcf().axes[:len(list_features[mask])]) #bkgs

    #-- debug
    # print(len(df[df['class']==0])); print(len(df[df['class']==1]))
    # print(df[df['class']==0].mean(axis=0)); print(df[df['class']==1].mean(axis=0))
    # print(df[df['class']==0].iloc[:3]); print(df[df['class']==1].iloc[:3])

    if isControlNorm == True: #Control plot, different name, general plot only
        plotname = weight_dir + 'InputFeatures_normTrain.png'
        print(colors.fg.lightgrey, "\nSaved input features plot [using normalized training set] as :", colors.reset, plotname)

    else:
        plotname = weight_dir + 'InputFeatures.png'
        print(colors.fg.lightgrey, "\nSaved input features plot as :", colors.reset, plotname)

    plt.savefig(plotname)
    plt.close()

    #-- Pairplot -- Plot pairwise relationships in a dataset. #Too much data, not relevant
    # sns.pairplot(df, vars=["maxDijetDelR", "dEtaFwdJetBJet"], diag_kind="kde", hue='class')
    # plotname = weight_dir + 'PairPlot.png'
    # plt.savefig(plotname)
    # print(colors.fg.lightgrey, "\nSaved pairwise relationship plot as :", colors.reset, plotname)

    #-- Also create individual plots
    # if plot_eachSingleFeature == True:
    if plot_eachSingleFeature == True and isControlNorm==False:
        for feature in list_features[mask]:
            fig = plt.figure('input_feature')
            _, ax = plt.subplots(figsize=(10,10))
            plt.ylabel("normalized", fontsize=20)
            plt.xlabel(feature, fontsize=20)

            x = [df[df['class'] == 1][feature].to_numpy(), df[df['class'] == 0][feature].to_numpy()]
            plt.hist(x[0], label='Signal', color='b', density=True, histtype='stepfilled', linewidth=2, bins=nbins, alpha=0.3, weights=df['weight'][df['class']==1])
            plt.hist(x[1], label='Backgrounds', color='r', density=True, histtype='stepfilled', linewidth=2, bins=nbins, alpha=0.3, weights=df['weight'][df['class']==0])
            plt.legend(loc='best')

            os.makedirs(weight_dir + 'features/', exist_ok=True)
            plotname = weight_dir + 'features/'+feature+'.png'
            if isControlNorm == True: plotname = weight_dir + 'features/'+feature+'_normTrain.png' #Control plot, different name, general plot only
            plt.savefig(plotname)
            print(colors.fg.lightgrey, "\nSaved input features plot as :", colors.reset, plotname)
            plt.close('input_feature')

    matplotlib.rc_file_defaults() #Restore matplotlib default settings

    return

# //--------------------------------------------
# //--------------------------------------------

   ##   #####   ####  #    # # ##### ######  ####  ##### #    # #####  ######
  #  #  #    # #    # #    # #   #   #      #    #   #   #    # #    # #
 #    # #    # #      ###### #   #   #####  #        #   #    # #    # #####
 ###### #####  #      #    # #   #   #      #        #   #    # #####  #
 #    # #   #  #    # #    # #   #   #      #    #   #   #    # #   #  #
 #    # #    #  ####  #    # #   #   ######  ####    #    ####  #    # ######

def Visualize_NN_architecture(weight_dir):
    '''
    Visualize the NN architecture using external packages.

    NB: Does not work yet...
    '''

    return

    '''
    with open(str(weight_dir+'arch_NN.json'), 'r') as json_file:
        model = model_from_json(json_file.read()) #model_from_json() expects a JSON string as its first parameter, not the name of a JSON file

    outname = weight_dir+'graphviz_NN.png'
    plot_model(model, to_file=outname, show_shapes=True, show_layer_names=True, dpi=96)
    print("\n-- Created NN arch plot with graphviz : " + outname)

    outname = weight_dir+'annviz_NN.gv'
    ann_viz(model, title="Neural network architecture", filename=outname, view=True) #cant handle batchnorm?
    print("\n-- Created NN arch plot with annviz : " + outname)

    return
    '''

# //--------------------------------------------
# //--------------------------------------------

def Make_SHAP_Plots(opts, model, weight_dir, list_xTrain_allClasses, list_xTest_allClasses, list_features):
    '''
    Use SHAP library to produce validation and control plots related to neural network.
    See: https://github.com/slundberg/shap

    NB: to truncate digits in plot, had to modify ~/Documents/Programmes/Anaconda3/lib/python3.7/site-packages/shap/plots/force_matplotlib.py !
    '''

# //--------------------------------------------

    plot_importance_allFeatures = True #True <-> create importance plot for each single input feature (very slow !)
    nmax=1000 #Max nof events *per process class* #None <-> use all events #SLOW !

# //--------------------------------------------

    if opts["nofOutputNodes"] > 1: return #Don't make SHAP plots for multiclass for now (less interpretable)
    # shap.initjs() #Load Javascript library #Obsolete
    shap.explainers.deep.deep_tf.op_handlers["AddV2"] = shap.explainers.deep.deep_tf.passthrough #Fix -- does not work on NAF ?
    if not isinstance(list_features, list): list_features = list_features.tolist() #Fix
    if opts['activInputLayer'] == 'lrelu' or opts['activHiddenLayers'] == 'lrelu': return #Not compatible

    #== Use 'Kernel SHAP' to explain test set predictions #Meant to approximate SHAP values for deep learning models.
    #-- model: The model to be explained. The output of the model can be a vector of size n_samples or a matrix of size [n_samples x n_output] (for a classification model).
    #-- data: Background dataset to generate the perturbed dataset required for training surrogate models. We simulate missing data by replacing the feature with the values it takes in the background dataset. So if the background dataset is a simple sample of all zeros, then we would approximate a feature being missing by setting it to zero. For small problems this background dataset can be the whole training set, but for larger problems consider using a single reference value or using the kmeans function to summarize the dataset.
    #-- link: A function to connect feature contribution values to the model output. For a classification model, we generally explain the logit of the predicted probability as a sum of feature contributions. Hence, if the output of the model (the first argument) is a probability, we set link = 'logit' to get the feature contributions in logit form.
    explainer = shap.DeepExplainer(model, data=np.concatenate([list[:nmax] for list in list_xTrain_allClasses]))
    # print('[SHAP] Number of classes: ', len(explainer.expected_value.numpy()))
    # print('[SHAP] Base value for first class:', explainer.expected_value.numpy()[0])

    #== Get SHAP values
    #-- X: Dataset on which to explain the model output.
    #-- nsamples: Nof samples to draw to build the surrogate model for explaining each prediction.
    # shap_values = explainer.shap_values(X=list_xTest_allClasses[0][:nmax]) #Returns a list of size n_classes. For binary classification model, n_classes=2. Each object of this list is an array of size [n_samples, n_features] and corresponds to the SHAP values for the respective class. For regression models, we get a single set of shap values of size [n_samples, n_features].
    shap_values = explainer.shap_values(X=np.concatenate([list[:nmax] for list in list_xTest_allClasses])) #Returns a list of size n_classes. For binary classification model, n_classes=2. Each object of this list is an array of size [n_samples, n_features] and corresponds to the SHAP values for the respective class. For regression models, we get a single set of shap values of size [n_samples, n_features].
    # shap_values = explainer.shap_values(X=list_xTest_allClasses[0][:100])#Returns a list of size n_classes. For binar classification model, n_classes=2. Each object of this list is an array of size [n_samples, n_features] and corresponds to the SHAP values for the respective class. For regression models, we get a single set of shap values of size [n_samples, n_features].
    # shap_values = explainer.shap_values(X=list_xTest_allClasses[0]) #Compute 'shap values'
    # print('[SHAP] Shape of shap value for each class: ', shap_values[0].shape)

    #== Plot SHAP values for first event, for specific investigation
    #-- See: https://github.com/slundberg/shap/blob/06c9d18f3dd014e9ed037a084f48bfaf1bc8f75a/shap/plots/force.py#L31
    #-- link: "identity" or "logit": the transformation used when drawing the tick mark labels. Using logit will change log-odds numbers into probabilities.
    #-- SHAP values are all relative to some base value. By default, the base value is explainer.expected_value: the mean of the raw model predictions for the training data. Predictions will be made relative to the base value, which must be provided
    # fig = plt.figure('shap_singleEvent')
    # shap_singleEvent = shap.force_plot(explainer.expected_value[0].numpy(), shap_values[0][0,:], list_xTest_allClasses[0][0,:], matplotlib=True, show=False, feature_names=list_features, link='logit') #matplotlib=True <-> avoid javascript #Need show=False to save plot
    # plt.savefig(weight_dir+"shap_singleEvent.png", bbox_inches='tight', dpi=1000)
    # plt.close('shap_singleEvent')
    # print(colors.fg.lightgrey, "\nSaved shap_singleEvent plot as :", colors.reset, weight_dir+"shap_singleEvent.png")

    # Plot SHAP values for all events
    # fig = plt.figure('all')
    shap_all = shap.force_plot(explainer.expected_value[0].numpy(), shap_values[0], np.concatenate([list[:nmax] for list in list_xTest_allClasses]).round(2), show=False, feature_names=list_features, link='logit', text_rotation=90) #matplotlib=True not supported for multiple events
    shap.save_html(weight_dir+'shap_all.html', shap_all)
    # plt.savefig("./shap_all.png", bbox_inches='tight', dpi=1000)
    # plt.close('all')
    print(colors.fg.lightgrey, "Saved shap_all plot as :", colors.reset, weight_dir+"shap_all.html")

    #-- Feature importance plot
    # See: https://github.com/slundberg/shap/blob/master/shap/plots/_summary.py#L18
    fig = plt.figure('shap_summary')
    if opts["nofOutputNodes"]==1: shap_summary = shap.summary_plot(shap_values=np.concatenate(shap_values), features=np.concatenate([list[:nmax] for list in list_xTest_allClasses]), show=False, feature_names=list_features, max_display=25)
    else: shap_summary = shap.summary_plot(shap_values=shap_values, features=np.concatenate([list[:nmax] for list in list_xTest_allClasses]), show=False, feature_names=list_features, max_display=25)
    plt.savefig(weight_dir+"shap_summary.png", bbox_inches='tight', dpi=1000)
    plt.close('shap_summary')
    print(colors.fg.lightgrey, "Saved shap_summary plot as :", colors.reset, weight_dir+"shap_summary.png")

    fig = plt.figure('shap_summary_bar')
    shap_summary_bar = shap.summary_plot(plot_type="bar", shap_values=np.concatenate(shap_values), features=np.concatenate([list[:nmax] for list in list_xTest_allClasses]), show=False, feature_names=list_features, max_display=30)
    plt.savefig(weight_dir+"shap_summary_bar.png", bbox_inches='tight', dpi=1000)
    plt.close('shap_summary_bar')
    print(colors.fg.lightgrey, "Saved shap_summary_bar plot as :", colors.reset, weight_dir+"shap_summary_bar.png")

    #-- Identical to 'shap_summary'
    # fig = plt.figure('shap_summary_singleClass')
    # shap_summary_singleClass = shap.summary_plot(shap_values[0], np.concatenate([list[:nmax] for list in list_xTest_allClasses]), show=False, feature_names=list_features, max_display=20)
    # plt.savefig(weight_dir+"shap_summary_singleClass.png", bbox_inches='tight', dpi=600)
    # plt.close('shap_summary_singleClass')
    # print(colors.fg.lightgrey, "Saved shap_summary_singleClass plot as :", colors.reset, weight_dir+"shap_summary_singleClass.png")

    #See: https://github.com/slundberg/shap/blob/master/shap/plots/decision.py#L217 #Very slow, not very useful
    # fig = plt.figure('decision_plot')
    # shap.decision_plot(explainer.expected_value[0].numpy(), shap_values[0], np.concatenate([list[:nmax] for list in list_xTest_allClasses]), feature_names=list_features, feature_order="importance", link='logit', show=False)
    # plt.savefig(weight_dir+"decision_plot.png", bbox_inches='tight', dpi=600)
    # plt.close('decision_plot')
    # print(colors.fg.lightgrey, "Saved decision_plot plot as :", colors.reset, weight_dir+"decision_plot.png")

    #== Multi-output decision plot #Not very useful (plot for 1 event only... ?)
    #-- See: https://github.com/slundberg/shap/blob/a4fd466193c7f9602e948e7f9fd65d49249ba4bd/shap/plots/decision.py#L553
    # row_index = 0
    # fig = plt.figure('multioutput_decision_plot')
    # shap.multioutput_decision_plot(explainer.expected_value.numpy().tolist(), shap_values, row_index=row_index, feature_names=list_features, highlight=None, legend_location='lower right', show=False)
    # plt.savefig(weight_dir+"multioutput_decision_plot.png", bbox_inches='tight', dpi=600)
    # plt.close('multioutput_decision_plot')
    # print(colors.fg.lightgrey, "Saved multioutput_decision_plot plot as :", colors.reset, weight_dir+"multioutput_decision_plot.png")

    #== Dependence of decision on single feature
    #-- See: https://github.com/slundberg/shap/blob/master/shap/plots/dependence.py#L15
    if plot_importance_allFeatures is True: #Make importance plots for all features (slow)
        for feature in list_features:
            fig = plt.figure('dependence_plot')
            shap.dependence_plot(feature, shap_values[0], np.concatenate([list[:nmax] for list in list_xTest_allClasses]), feature_names=list_features, alpha=0.5, interaction_index=None, show=False) #'interaction_index=inds[i]' shows interaction with 2nd variable on z-axis
            plt.savefig(weight_dir+"dependence_plot_"+feature+".png", bbox_inches='tight', dpi=600)
            plt.close('dependence_plot')
            print(colors.fg.lightgrey, "Saved dependence_plot plot as :", colors.reset, weight_dir+"dependence_plot_"+feature+".png")

    else: #Only plot few specific features (hard-coded names)
        feature1 = 'ctw'
        if feature1 in list_features:
            fig = plt.figure('dependence_plot')
            shap.dependence_plot(feature1, shap_values[0], np.concatenate([list[:nmax] for list in list_xTest_allClasses]), feature_names=list_features, alpha=0.5, interaction_index=None, show=False) #'interaction_index=inds[i]' shows interaction with 2nd variable on z-axis
            plt.savefig(weight_dir+"dependence_plot.png", bbox_inches='tight', dpi=600)
            plt.close('dependence_plot')
            print(colors.fg.lightgrey, "Saved dependence_plot plot as :", colors.reset, weight_dir+"dependence_plot.png")

        #Second dependence plot
        feature2 = "recoZ_Pt"
        if feature2 in list_features:
            fig = plt.figure('dependence_plot2')
            shap.dependence_plot(feature2, shap_values[0], np.concatenate([list[:nmax] for list in list_xTest_allClasses]), feature_names=list_features, alpha=0.5, interaction_index=None, show=False) #'interaction_index=inds[i]' shows interaction with 2nd variable on z-axis
            plt.savefig(weight_dir+"dependence_plot2.png", bbox_inches='tight', dpi=600)
            plt.close('dependence_plot2')
            print(colors.fg.lightgrey, "Saved dependence_plot2 plot as :", colors.reset, weight_dir+"dependence_plot2.png")

        #Third dependence plot
        feature3 = "recoZ_Eta"
        if feature3 in list_features:
            fig = plt.figure('dependence_plot3')
            shap.dependence_plot(feature3, shap_values[0], np.concatenate([list[:nmax] for list in list_xTest_allClasses]), feature_names=list_features, alpha=0.5, interaction_index=None, show=False)
            plt.savefig(weight_dir+"dependence_plot2.png", bbox_inches='tight', dpi=600)
            plt.close('dependence_plot3')
            print(colors.fg.lightgrey, "Saved dependence_plot3 plot as :", colors.reset, weight_dir+"dependence_plot3.png")

    return

# //--------------------------------------------
# //--------------------------------------------

def Make_Test1D_Plot(opts, model):
    '''
    Plot NN output as a function of dummy input feature x, for several values of parameter theta.
    Similar to what is done in reference paramNN article.
    '''

    nPts = 500
    fig = plt.figure('TEST')
    cm = plt.cm.get_cmap('RdYlBu_r')

    #Best gaus fit
    def gaus(x,a,b): return np.exp(-((x-a)/b)**2)

    for theta in [-4, 0, 4]:
        x = [xx for xx in np.linspace(-10, 10, num=nPts)] #Feature x
        y = np.squeeze(model.predict(np.append(np.atleast_2d(x).T, np.atleast_2d(np.full(nPts,theta)).T, axis=1))) #NN response (depends on x and theta)

        plt.scatter(x,y, alpha=0.50, label=str(theta))
        # params, params_covariance = optimize.curve_fit(gaus, x, y, p0=[theta,0.5])
        # plt.plot(x, gaus(x, params[0], params[1]), label='theta='+str(theta))
        # print(params[0]); print(params[1])

    plt.legend(loc='best', numpoints=1, title='Theta param.')
    # plt.title('xxx')
    plt.ylabel('NN output')
    plt.xlabel('1D input feature')

    plotname = './TEST2.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved TEST2 plot as :", colors.reset, plotname)
    fig.clear()
    plt.close('TEST')

    return

# //--------------------------------------------
# //--------------------------------------------

# Plot the score... testing
def Test_Make_Score_Plot(weight_dir, scores_allClasses_eachOperator, y_process, x):

    fig = plt.figure('loss')

    # print(x[:,-1])

    nbins = 50
    rmin = -1; rmax = 1
    plt.hist(scores_allClasses_eachOperator[y_process==0], bins=nbins, range=(rmin,rmax), color= 'cornflowerblue', alpha=0.50, density=True, histtype='step')
    plt.hist(scores_allClasses_eachOperator[y_process==1], bins=nbins, range=(rmin,rmax), color= 'darkorange', alpha=0.50, density=True, histtype='step')
    # plt.hist(scores_allClasses_eachOperator[np.logical_and(y_process==0, x[:,-1]==5)], bins=nbins, range=(rmin,rmax), color= 'cornflowerblue', alpha=0.50, density=True, histtype='step')
    # plt.hist(scores_allClasses_eachOperator[np.logical_and(y_process==1, x[:,-1]==5)], bins=nbins, range=(rmin,rmax), color= 'darkorange', alpha=0.50, density=True, histtype='step')

    plotname = weight_dir + 'TEST.png'
    fig.savefig(plotname, bbox_inches='tight') #bbox_inches='tight' ensures that second y-axis is visible
    print(colors.fg.lightgrey, "\nSaved Loss plot as :", colors.reset, plotname)
    fig.clear()
    plt.close('loss')

    return

# //--------------------------------------------
# //--------------------------------------------
