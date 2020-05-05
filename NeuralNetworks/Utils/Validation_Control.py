#Create control plots, ROC histos, etc.

import os
import ROOT
import numpy as np
import keras
import math
import json
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from tensorflow.keras.models import load_model, model_from_json
from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
from root_numpy import fill_hist
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score, classification_report, confusion_matrix, multilabel_confusion_matrix
# from sklearn.feature_selection import mutual_info_classif
from Utils.Helper import close_event, Printout_Outputs_FirstLayer, KS_test, Anderson_Darling_test, ChiSquare_test
from Utils.ColoredPrintout import colors
from pandas.plotting import scatter_matrix
from ann_visualizer.visualize import ann_viz
from tensorflow.keras.utils import plot_model

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

#Apply NN model on train/test datasets to produce ROOT histograms which can later be used to plot ROC curves
def Store_TrainTestPrediction_Histograms(opts, lumiName, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses):

    print(colors.fg.lightblue, "--- Create & store ROC histos...", colors.reset); print('\n')

    maxEvents = 500000 #Upper limit on nof events per class, else validation too slow (problematic for parameterized NN with huge training stat.)

    # Fill a ROOT histogram from NumPy arrays, with fine binning (no loss of info)
    rootfile_outname = "../outputs/NN_"+list_labels[0]+"_"+lumiName+".root"
    fout = ROOT.TFile(rootfile_outname, "RECREATE")

    nofOutputNodes = opts["nofOutputNodes"]

    nodes_labels = list_labels
    if opts["strategy"] == "RASCAL":
        nodes_labels = ["LR"]
        for op in opts["listOperatorsParam"]: nodes_labels.append(str("Score "+op))

    for inode in range(nofOutputNodes):

        hist_TrainingEvents_allClasses = TH1F('hist_train_NODE_'+nodes_labels[inode]+'_allClasses', '', 1000, 0, 1); hist_TrainingEvents_allClasses.Sumw2()
        hist_TestingEvents_allClasses = TH1F('hist_test_NODE_'+nodes_labels[inode]+'_allClasses', '', 1000, 0, 1); hist_TestingEvents_allClasses.Sumw2()

        for iclass in range(len(list_labels)):
            # print('inode', inode, 'iclass', iclass)

            hist_TrainingEvents_class = TH1F('hist_train_NODE_'+nodes_labels[inode]+'_CLASS_'+list_labels[iclass], '', 1000, 0, 1); hist_TrainingEvents_class.Sumw2()
            fill_hist(hist_TrainingEvents_class, list_predictions_train_allNodes_allClasses[inode][iclass][:maxEvents], weights=list_PhysicalWeightsTrain_allClasses[iclass][:maxEvents])
            hist_TrainingEvents_class.Write()

            hist_TestingEvents_class = TH1F('hist_test_NODE_'+nodes_labels[inode]+'_CLASS_'+list_labels[iclass], '', 1000, 0, 1); hist_TestingEvents_class.Sumw2()
            fill_hist(hist_TestingEvents_class, list_predictions_test_allNodes_allClasses[inode][iclass][:maxEvents], weights=list_PhysicalWeightsTest_allClasses[iclass][:maxEvents])
            hist_TestingEvents_class.Write()

            fill_hist(hist_TrainingEvents_allClasses, list_predictions_train_allNodes_allClasses[inode][iclass][:maxEvents], weights=list_PhysicalWeightsTrain_allClasses[iclass][:maxEvents])
            fill_hist(hist_TestingEvents_allClasses, list_predictions_test_allNodes_allClasses[inode][iclass][:maxEvents], weights=list_PhysicalWeightsTest_allClasses[iclass][:maxEvents])

        hist_TrainingEvents_allClasses.Write()
        hist_TestingEvents_allClasses.Write()

    fout.Close()
    # print("Saved output ROOT file containing Keras Predictions as histograms : " + rootfile_outname)
    print(colors.fg.lightgrey, "\nSaved output ROOT file containing Keras Predictions as histograms :", colors.reset, rootfile_outname, '\n')

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

# Call all sub-functions
def Make_Default_Validation_Plots(opts, list_features, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, PhysicalWeights_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, x, y_train, y_test, y_process, y_process_train, y_process_test, metrics, weight_dir, score, history=None):

    print('\n'); print(colors.fg.lightblue, "--- Create control plots...", colors.reset); print('\n')

    #NB: styling not working properly ? related to figsize... ?
    # sns.set(palette='coolwarm', font_scale=1.4) #Scale up label font size #NB: this also sets plotting options to seaborn's default
    # plt.tight_layout()

    Control_Printouts(opts, score, list_labels, y_test, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTest_allClasses)

    Make_Loss_Plot(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, weight_dir, history)

    Make_Metrics_Plot(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, metrics, weight_dir, history)

    Make_ROC_plots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, weight_dir)

    Make_Overtraining_plots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, weight_dir)

    if opts["regress"] == True: Make_Regressor_ControlPlots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, y_train, y_test, y_process_train, y_process_test, weight_dir)

    Create_Correlation_Plot(opts, x, list_features, weight_dir)

    Plot_Input_Features(opts, x, y_process, PhysicalWeights_allClasses, list_features, weight_dir, False)

    return

# //--------------------------------------------
# //--------------------------------------------

 #####  #####  # #    # #####  ####  #    # #####
 #    # #    # # ##   #   #   #    # #    #   #
 #    # #    # # # #  #   #   #    # #    #   #
 #####  #####  # #  # #   #   #    # #    #   #
 #      #   #  # #   ##   #   #    # #    #   #
 #      #    # # #    #   #    ####   ####    #

#Printout some information related to NN performance
def Control_Printouts(opts, score, list_labels, y_test, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTest_allClasses):

    if opts["nofOutputNodes"] == 1:
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

def Make_Loss_Plot(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, weight_dir, history=None):

    # NB : Possible solution to show 3 y-axes (train, test, lr) : https://stackoverflow.com/questions/48618992/matplotlib-graph-with-more-than-2-y-axes

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
    plt.legend(lns, labs, loc='upper right')

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

def Make_Metrics_Plot(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, metrics, weight_dir, history=None):

    if opts["strategy"] == "RASCAL": #Metrics name depend on the output (e.g. 'val_likelihood_ratio_mean_squared_error' or 'val_score_mean_squared_error')
        metrics = "likelihood_ratio_" + metrics

    # Plotting the error with the number of iterations
    fig = plt.figure('metric')
    ax1 = fig.gca()
    timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
    timer.add_callback(close_event)

    # plt.plot(history.history['val_'+metrics], color='cornflowerblue')
    plt.title('Accuracy VS Epoch')
    plt.ylabel('Train '+metrics, color='darkorange')
    plt.xlabel('Epoch')
    lns1 = plt.plot(history.history[metrics], color='darkorange', label='Train') #metrics name
    # print(history.history.keys())

    ax2 = ax1.twinx()  # instantiate a second axis that shares the same x-axis
    ax2.set_ylabel('Test '+metrics, color='cornflowerblue')  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', color='cornflowerblue')
    lns2 = ax2.plot(history.history['val_'+metrics], color='cornflowerblue', label='Test')

    ax3 = ax1.twinx()  # instantiate a 3rd axis that shares the same x-axis
    ax3.get_yaxis().set_ticks([]) #invisible y axis (not important)
    lns3 = ax3.plot(history.history['lr'], color='dimgrey', linestyle='--', label='lr')

    #Legend
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc='lower center')

    timer.start()
    plt.show()

    plotname = weight_dir + 'Accuracy_NN_.png'
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

def Make_ROC_plots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, weight_dir):

    nofOutputNodes = opts["nofOutputNodes"]

    for inode in range(nofOutputNodes):

        lw = 2 #linewidth

        if opts["strategy"] in ["ROLR", "RASCAL"]: #Transformation r->s : s = 1/(r+1)

            if inode != 0: continue #Only for r node

            fpr, tpr, _ = roc_curve(np.concatenate(list_truth_Test_allClasses), 1./(np.concatenate(list_predictions_test_allNodes_allClasses[inode])+1))
            roc_auc = auc(fpr, tpr)

            fpr_train, tpr_train, _ = roc_curve(np.concatenate(list_truth_Train_allClasses), 1./(np.concatenate(list_predictions_train_allNodes_allClasses[inode])+1))
            roc_auc_train = auc(fpr_train, tpr_train)

        elif opts["nofOutputNodes"] == 1:

            fpr, tpr, _ = roc_curve(np.concatenate(list_truth_Test_allClasses), np.concatenate(list_predictions_test_allNodes_allClasses[inode]))
            roc_auc = auc(fpr, tpr)

            fpr_train, tpr_train, _ = roc_curve(np.concatenate(list_truth_Train_allClasses), np.concatenate(list_predictions_train_allNodes_allClasses[inode]))
            roc_auc_train = auc(fpr_train, tpr_train)

        elif opts["nofOutputNodes"] > 1: # Multiclass: 1 ROC (signal vs All) per node

            fpr = dict()
            tpr = dict()
            roc_auc = dict()
            fpr_train = dict()
            tpr_train = dict()
            roc_auc_train = dict()

            fpr[inode], tpr[inode], _ = roc_curve(np.concatenate(list_truth_Test_allClasses)[:,inode], np.concatenate(list_predictions_test_allNodes_allClasses[inode]))
            roc_auc[inode] = auc(fpr[inode], tpr[inode])

            fpr_train[inode], tpr_train[inode], _ = roc_curve(np.concatenate(list_truth_Train_allClasses)[:,inode], np.concatenate(list_predictions_train_allNodes_allClasses[inode]))
            roc_auc_train[inode] = auc(fpr_train[inode], tpr_train[inode])

        # Plot ROC curves
        fig = plt.figure()
        timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
        timer.add_callback(close_event)

        if opts["nofOutputNodes"] == 1 or opts["strategy"] in ["ROLR", "RASCAL"]: #single node
            plt.plot(tpr_train, 1-fpr_train, color='darkorange', lw=lw, label='ROC NN (train) (AUC = {1:0.2f})' ''.format(0, roc_auc_train))
            plt.plot(tpr, 1-fpr, color='cornflowerblue', lw=lw, label='ROC NN (test) (AUC = {1:0.2f})' ''.format(0, roc_auc))
        else: #for each node
            plt.plot(tpr_train[inode], 1-fpr_train[inode], color='darkorange', lw=lw, label='ROC NN (train) (AUC = {1:0.2f})' ''.format(0, roc_auc_train[inode]))
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
        plt.legend(loc="lower left")

        #Display plot in terminal for quick check (only for first node)
        if inode == 0:
            timer.start()
            plt.show()

        plotname = weight_dir + 'ROC_NN_' + list_labels[inode] + '.png'
        fig.savefig(plotname)
        print(colors.fg.lightgrey, "\nSaved ROC plot as :", colors.reset, plotname)
        fig.clear()
        plt.close(fig)

    #Multiclass ROC (first node only)
    if (opts["strategy"] is "classifier" or opts["strategy"] is "CARL_multiclass") and nofOutputNodes > 1:

        # Compute ROC curve and ROC area for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        fpr_train = dict()
        tpr_train = dict()
        roc_auc_train = dict()

        # Plot ROC curves
        fig = plt.figure('multiclass')
        timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
        timer.add_callback(close_event)

        for iclass in range(len(list_labels)):

            # print(list_predictions_test_allNodes_allClasses[0][iclass])
            # print(list_truth_Test_allClasses[iclass][:,0])
            # print(list_predictions_test_allNodes_allClasses[0][iclass].shape)
            # print(list_truth_Test_allClasses[iclass][:,0].shape)

            if iclass == 0: #first signal vs all
                fpr[iclass], tpr[iclass], _ = roc_curve(np.concatenate(list_truth_Test_allClasses)[:,0], np.concatenate(list_predictions_test_allNodes_allClasses[0]))
                roc_auc[iclass] = auc(fpr[iclass], tpr[iclass])
                mylabel = 'vs All (AUC = {1:0.2f})' ''.format(0, roc_auc[iclass])
            else: #first signal vs specific process
                fpr[iclass], tpr[iclass], _ = roc_curve(np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[iclass]))[:,0], np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][iclass])) )
                roc_auc[iclass] = auc(fpr[iclass], tpr[iclass])
                mylabel = 'vs ' + list_labels[iclass]+' (AUC = {1:0.2f})' ''.format(0, roc_auc[iclass])

            # fpr_train[iclass], tpr_train[iclass], _ = roc_curve(np.concatenate(list_truth_Train_allClasses)[:,0], np.concatenate(list_predictions_train_allNodes_allClasses[0][iclass]))
            # roc_auc_train[iclass] = auc(fpr_train[iclass], tpr_train[iclass])

            # plt.plot(tpr_train[iclass], 1-fpr_train[iclass], color='darkorange', lw=lw, label='ROC NN (train) (AUC = {1:0.2f})' ''.format(0, roc_auc_train[iclass]))
            # plt.plot(tpr[iclass], 1-fpr[iclass], color='cornflowerblue', lw=lw, label='ROC NN (test) (AUC = {1:0.2f})' ''.format(0, roc_auc[iclass]))

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
        plt.title('ROC curves for '+list_labels[0]+ ' (test data)')
        plt.legend(loc="lower left")

        timer.start()
        plt.show()

        plotname = weight_dir + 'ROC_NN_node0_allClasses.png'
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

def Make_Overtraining_plots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, weight_dir):

    nofOutputNodes = opts["nofOutputNodes"]

    for inode in range(nofOutputNodes):

        if opts["strategy"] in ["ROLR", "RASCAL"] and inode > 0: continue #Only for r node

        nbins = 20
        rmin = 0.
        rmax = 1.

        fig = plt.figure('overtrain')
        timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
        timer.add_callback(close_event)

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
                # tick.set_weight('bold')
                for tick in ax.get_yticklabels():
                    tick.set_color('gray')

        #-- Trick : for training histos, we want to compute the bin errors correctly ; to do this we first fill TH1Fs, then read their bin contents/errors
        hist_overtrain_train_sig = TH1F('hist_overtrain_train_sig', '', nbins, rmin, rmax); hist_overtrain_train_sig.Sumw2(); hist_overtrain_train_sig.SetDirectory(0)
        hist_overtrain_train_bkg = TH1F('hist_overtrain_train_bkg', '', nbins, rmin, rmax); hist_overtrain_train_bkg.Sumw2(); hist_overtrain_train_bkg.SetDirectory(0)

        #only signal process
        for val, w in zip((list_predictions_train_allNodes_allClasses[inode][inode]), list_PhysicalWeightsTrain_allClasses[inode]): #'zip' stops when the shorter of the lists stops
            if opts["strategy"] in ["ROLR", "RASCAL"]: val = 1/(val+1) #Transform r -> s
            hist_overtrain_train_sig.Fill(val, w)

        #only background processes
        for iclass in range(0, len(list_labels)):
            if iclass != inode:
                for val, w in zip(list_predictions_train_allNodes_allClasses[inode][iclass], list_PhysicalWeightsTrain_allClasses[iclass]):
                    if opts["strategy"] in ["ROLR", "RASCAL"]: val = 1/(val+1) #Transform r -> s
                    hist_overtrain_train_bkg.Fill(val, w)

        #Normalize
        sf_integral = abs(rmax - rmin) / nbins #h.Scale(1/integral) makes the sum of contents equal to 1, but does not account for the bin width
        hist_overtrain_train_sig.Scale(1./(hist_overtrain_train_sig.Integral()*sf_integral))
        hist_overtrain_train_bkg.Scale(1./(hist_overtrain_train_bkg.Integral()*sf_integral))

        #Plot testing sig/bkg histos, normalized (no errors displayed <-> don't need TH1Fs)
        if opts["strategy"] in ["ROLR", "RASCAL"]: tmp = 1/(list_predictions_test_allNodes_allClasses[inode][inode]+1) #Transform r -> s
        else: tmp = list_predictions_test_allNodes_allClasses[inode][inode]
        plt.hist(tmp, bins=nbins, range=(rmin,rmax), color= 'cornflowerblue', alpha=0.50, weights=list_PhysicalWeightsTest_allClasses[inode], density=True, histtype='step', log=False, label="Signal (Test)", edgecolor='cornflowerblue',fill=True)

        #Trick : want to get arrays of predictions/weights for all events *which do not belong to class of current node* => Loop on classes, check if matches node
        lists_predictions_bkgs = []; lists_weights_bkg = []
        for iclass in range(0, len(list_labels)):
            if iclass != inode:
                if opts["strategy"] in ["ROLR", "RASCAL"]: tmp = 1/(list_predictions_test_allNodes_allClasses[inode][iclass]+1) #Transform r -> s
                else: tmp = list_predictions_test_allNodes_allClasses[inode][iclass]
                lists_predictions_bkgs.append(tmp)
                # lists_predictions_bkgs.append(list_predictions_test_allNodes_allClasses[inode][iclass])
                lists_weights_bkg.append(list_PhysicalWeightsTest_allClasses[iclass])
        predictions_bkgs = np.concatenate(lists_predictions_bkgs); weights_bkgs = np.concatenate(lists_weights_bkg)

        plt.hist(predictions_bkgs, bins=nbins, range=(rmin,rmax), color='orangered', alpha=0.50, weights=weights_bkgs, density=True, histtype='step', log=False, label="Background (Test)", hatch='/', edgecolor='orangered',fill=False)

        #Read bin contents/errors
        bin_centres = []; counts_sig = []; err_sig = []; counts_bkg = []; err_bkg = []
        for ibin in range(1, hist_overtrain_train_sig.GetNbinsX()+1):
            bin_centres.append(hist_overtrain_train_sig.GetBinCenter(ibin))
            counts_sig.append(hist_overtrain_train_sig.GetBinContent(ibin)); counts_bkg.append(hist_overtrain_train_bkg.GetBinContent(ibin))
            err_sig.append(hist_overtrain_train_sig.GetBinError(ibin)); err_bkg.append(hist_overtrain_train_bkg.GetBinError(ibin))

        #Plot training sig/bkg histos, normalized, with errorbars
        plt.errorbar(bin_centres, counts_sig, marker='o', yerr=err_sig, linestyle='None', markersize=6, color='blue', alpha=0.90, label='Signal (Train)')
        plt.errorbar(bin_centres, counts_bkg, marker='o', yerr=err_bkg, linestyle='None', markersize=6, color='red', alpha=0.90, label='Background (Train)')

        myxlabel = "Classifier output"

        plt.legend(loc='upper center', numpoints=1)
        plt.title("Output distributions for Signal & Background")
        plt.grid(axis='y', alpha=0.75)
        plt.grid(axis='x', alpha=0.75)
        plt.xlabel(myxlabel)
        plt.ylabel('PDF')

        if inode == 0:
            timer.start()
            plt.show()

        plotname = weight_dir + 'Overtraining_NN_' + list_labels[inode] + '.png'
        fig.savefig(plotname)
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

def Make_Regressor_ControlPlots(opts, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_PhysicalWeightsTest_allClasses, list_truth_Train_allClasses, list_truth_Test_allClasses, y_train, y_test, y_process_train, y_process_test, weight_dir):

    nbins = 20
    rmin = 0; rmax = 1
    nofOutputNodes = opts["nofOutputNodes"]

    for inode in range(nofOutputNodes):

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

        fig = plt.figure('regressor')
        timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
        timer.add_callback(close_event)

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
                # tick.set_weight('bold')
                for tick in ax.get_yticklabels():
                    tick.set_color('gray')

        #-- Trick : for truth histos, we want to compute the bin errors correctly ; to do this we first fill TH1Fs, then read their bin contents/errors
        hist_truth_SM = TH1F('hist_truth_SM', '', nbins, rmin, rmax); hist_truth_SM.Sumw2(); hist_truth_SM.SetDirectory(0)
        hist_truth_EFT = TH1F('hist_truth_EFT', '', nbins, rmin, rmax); hist_truth_EFT.Sumw2(); hist_truth_EFT.SetDirectory(0)

        for val in y_test[y_process_test==0]: #'zip' stops when the shorter of the lists stops
            if opts["strategy"] is "RASCAL": val = val[inode]
            hist_truth_SM.Fill(val, 1.)

        #only background processes
        for val in y_test[y_process_test==1]: #'zip' stops when the shorter of the lists stops
            if opts["strategy"] is "RASCAL": val = val[inode]
            hist_truth_EFT.Fill(val, 1.)

        #Normalize
        sf_integral = abs(rmax - rmin) / nbins #h.Scale(1/integral) makes the sum of contents equal to 1, but does not account for the bin width
        hist_truth_SM.Scale(1./(hist_truth_SM.Integral()*sf_integral))
        hist_truth_EFT.Scale(1./(hist_truth_EFT.Integral()*sf_integral))

        plt.hist(np.concatenate(list_predictions_test_allNodes_allClasses[inode])[y_process_test==0], bins=nbins, range=(rmin,rmax), color= 'cornflowerblue', alpha=0.50, density=True, histtype='step', log=False, label="SM (Pred.)", edgecolor='cornflowerblue',fill=True)
        plt.hist(np.concatenate(list_predictions_test_allNodes_allClasses[inode])[y_process_test==1], bins=nbins, range=(rmin,rmax), color= 'orangered', alpha=0.50, density=True, histtype='step', log=False, label="EFT (Pred.)", edgecolor='cornflowerblue',fill=True)

        #Read bin contents/errors
        bin_centres = []; counts_SM = []; err_SM = []; counts_EFT = []; err_EFT = []
        for ibin in range(1, hist_truth_SM.GetNbinsX()+1):
            bin_centres.append(hist_truth_SM.GetBinCenter(ibin))
            counts_SM.append(hist_truth_SM.GetBinContent(ibin)); counts_EFT.append(hist_truth_EFT.GetBinContent(ibin))
            err_SM.append(hist_truth_SM.GetBinError(ibin)); err_EFT.append(hist_truth_EFT.GetBinError(ibin))

        #Plot truth histos, normalized, with errorbars
        plt.errorbar(bin_centres, counts_SM, marker='o', yerr=err_SM, linestyle='None', markersize=6, color='blue', alpha=0.90, label='SM (Truth)')
        plt.errorbar(bin_centres, counts_EFT, marker='o', yerr=err_EFT, linestyle='None', markersize=6, color='red', alpha=0.90, label='EFT (Truth)')

        nodename = 'node' + str(inode)
        if opts["strategy"] in ["ROLR", "RASCAL"] and inode == 0: nodename = 'LR'
        elif opts["strategy"] is "RASCAL": nodename = 'Score_' + opts["listOperatorsParam"][inode-1] #First node is r

        # myxlabel = "Regressor output"
        plt.legend(loc='upper center', numpoints=1)
        plt.title("Regressor output for SM & EFT")
        plt.grid(axis='y', alpha=0.75)
        plt.grid(axis='x', alpha=0.75)
        plt.xlabel(nodename)
        plt.ylabel('PDF')

        if inode == 0:
            timer.start()
            plt.show()

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

#Plot correlations between input features
def Create_Correlation_Plot(opts, x, list_features, weight_dir):

    #Use this variables to avoid plotting also Wilson Coefficients given as inputs
    iVarMax = None
    if opts["parameterizedNN"] == True: iVarMax = -len(opts["listOperatorsParam"])

    #-- Convert np array to pd dataframe
    df = pd.DataFrame(data=x[0:,0:iVarMax], columns=list_features[:iVarMax]) #x = (events, vars) ; colums names are var names
    # print(df); print(df.describe())

    #-- Get correlation matrix
    corr = df.corr()

    mask = np.triu(np.ones_like(corr, dtype=np.bool)) #Mask upper right triangle
    # mask = np.tril(np.ones_like(corr, dtype=np.bool)) #Mask bottom left triangle

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
    hm = sns.heatmap(corr, mask=mask, cmap=palette, vmin=-1., vmax=1., center=0, square=True, linewidths=0.5, annot = True, fmt='.1g', cbar_kws={"shrink": .5},)
    sns.set(font_scale=1.4) #Scale up label font size #NB : also sets plotting options to seaborn's default
    hm.set_yticklabels(hm.get_yticklabels(), fontsize = 16)
    # ax.set_ylim(bottom=-0.5, top=len(list_features)+0.5)
    # ax.set_xlim(left=-0.5, right=len(list_features)+0.5)
    # ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False) # x labels appear on top
    # hm.set_xticklabels(hm.get_xticklabels(), rotation=45, horizontalalignment='left', va='bottom', fontsize = 16)
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True) # x labels appear at bottom
    hm.set_xticklabels(hm.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize = 16) #bottom labels
    fig.tight_layout()

    plotname = weight_dir + 'CorrelMatrix.png'
    plt.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved correlation matrix plot as :", colors.reset, plotname)
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

#Plot distributions of input features (according to dataset in arg)
def Plot_Input_Features(opts, x, y_process, weights, list_features, weight_dir, isControlNorm=False):

    plot_eachSingleFeature = False #True <-> 1 single plot per feature
# //--------------------------------------------

    # if opts["parameterizedNN"] is True and isControlNorm is True: return #Gain time
    if x is None: print('Error, can\'t produce input features plots : x=None !'); return

    sns.set(palette='coolwarm', font_scale=1.4) #Scale up label font size #NB : this also sets plotting options to seaborn's default
    plt.tight_layout()

    #Use this variables to avoid plotting also Wilson Coefficients given as inputs
    iVarMax = None
    n_features = len(list_features)
    if opts["parameterizedNN"] == True:
        iVarMax = -len(opts["listOperatorsParam"])
        n_features+= iVarMax
        nMax = 50000 #Don't use more events (slow) #Warning : for parameterized NN, biases distributions of WCs (will most likely not cover all EFT points)
        if len(x) > nMax: x = x[:nMax]; y_process = y_process[:nMax]; weights = weights[:nMax] #Else, too slow

    #-- Convert np array to pd dataframe
    df = pd.DataFrame(data=x[0:,0:iVarMax], columns=list_features[:iVarMax]) #x = (events, vars) ; colums names are var names

    #-- Insert a column corresponding to the class label
    if y_process.ndim==1: df.insert(loc=0, column='class', value=y_process[:], allow_duplicates=False)
    else: df.insert(loc=0, column='class', value=y_process[:,0], allow_duplicates=False) #Only care about first column=main signal (rest -> bkg)

    #-- Insert a column corresponding to physical event weights
    df.insert(loc=0, column='weight', value=weights[:,]) #Only care about first column=main signal (rest -> bkg)
    # print(df); print(df.describe())

    #-- Create multiplot #NB: only process columns corresponding to phy vars
    df[df['class']==1].hist(figsize=(15,15), label='Signal', column=list_features[:iVarMax], weights=df['weight'][df['class']==1], bins=20, alpha=0.4, density=True, color='r') #signal
    df[df['class']==0].hist(figsize=(15,15), label='Backgrounds', column=list_features[:iVarMax], weights=df['weight'][df['class']==0], bins=20, alpha=0.4, density=True, color='b', ax=plt.gcf().axes[:n_features]) #bkgs

    if isControlNorm == True: #Control plot, different name, general plot only
        plotname = weight_dir + 'InputFeatures_normTrain.png'
        print(colors.fg.lightgrey, "\nSaved input features plot [using normalized training set] as :", colors.reset, plotname)

    else:
        plotname = weight_dir + 'InputFeatures.png'
        print(colors.fg.lightgrey, "\nSaved input features plot as :", colors.reset, plotname)

    plt.savefig(plotname)
    plt.close('all')

    #-- Pairplot -- Plot pairwise relationships in a dataset. #Too much data, not relevant
    # sns.pairplot(df, vars=["maxDijetDelR", "dEtaFwdJetBJet"], diag_kind="kde", hue='class')
    # plotname = weight_dir + 'PairPlot.png'
    # plt.savefig(plotname)
    # print(colors.fg.lightgrey, "\nSaved pairwise relationship plot as :", colors.reset, plotname)

    #-- Also create individual plots
    # if plot_eachSingleFeature == True:
    if plot_eachSingleFeature == True and isControlNorm==False:
        for feature in list_features:
            # fig, ax = plt.subplots()
            fig, ax = plt.subplots(figsize=(10,10))
            plt.ylabel("normalized", fontsize=20)
            plt.xlabel(feature, fontsize=20)

            xx = [df[df['class'] == 1][feature].to_numpy(), df[df['class'] == 0][feature].to_numpy()]
            plt.hist(xx[0], label='Signal', color='r', density=True, histtype='stepfilled', linewidth=2, bins=20, alpha=0.4, weights=df['weight'][df['class']==1])
            plt.hist(xx[1], label='Backgrounds', color='b', density=True, histtype='stepfilled', linewidth=2, bins=20, alpha=0.4, weights=df['weight'][df['class']==0])
            plt.legend(loc="upper center")

            os.makedirs(weight_dir + 'features/', exist_ok=True)
            plotname = weight_dir + 'features/'+feature+'.png'
            if isControlNorm == True: #Control plot, different name, general plot only
                plotname = weight_dir + 'features/'+feature+'_normTrain.png'
            plt.savefig(plotname)
            plt.close('all')
            print(colors.fg.lightgrey, "\nSaved input features plot as :", colors.reset, plotname)

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

#Does not work yet...
def Visualize_NN_architecture(weight_dir):

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
