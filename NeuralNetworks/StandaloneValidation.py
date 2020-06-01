#xxx

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
import matplotlib.cm as cm
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score, classification_report, confusion_matrix, multilabel_confusion_matrix
from Utils.Helper import *
from Utils.ColoredPrintout import colors
from Utils.RegressorValidation import *
from pandas.plotting import scatter_matrix
from ann_visualizer.visualize import ann_viz
from tensorflow.keras.utils import plot_model
from Utils.GetData import Get_Data
from Train_Neural_Network import optsTrain,_list_lumiYears,_list_labels,_list_features,_list_processClasses
from Utils.Validation_Control import *
from Utils.Predictions import *
from Utils.RegressorValidation import *

# //--------------------------------------------
# //--------------------------------------------

nEventsStandaloneVal = 500 #Nof events to sample/display per point

#== SINGLE POINT AT WHICH TO EVALUATE EVENTS #NB: i.e. 'rwgt_ctW_3' corresponds to asking the NN 'are these events more EFT(ctW=3)-like, or more reference-like (<-> SM-like)'.
# evalPoint = "SM"
evalPoint = "rwgt_ctZ_3"
# evalPoint = "rwgt_ctZ_5"

#== LIST OF POINTS FROM WHICH TO SAMPLE EVENTS  #NB: order of operators should be the same as used for training #NB: for CARL_multiclass, only 1 operator can be activated per point !
list_points_sampling = []
# list_points_sampling.append("SM") #Keep this
list_points_sampling.append("rwgt_ctZ_1")
list_points_sampling.append("rwgt_ctZ_3")
list_points_sampling.append("rwgt_ctZ_5")
# list_points_sampling.append("rwgt_ctW_1")
# list_points_sampling.append("rwgt_ctW_2")
# list_points_sampling.append("rwgt_ctW_3")
# list_points_sampling.append("rwgt_ctW_4")
# list_points_sampling.append("rwgt_ctW_2_cpQ3_4.5")
# list_points_sampling.append("rwgt_ctZ_3_ctW_0_cpQM_0_cpQ3_0")
# list_points_sampling.append("rwgt_ctZ_3_ctW_0_cpQM_0_cpQ3_0_cpt_0")
# list_points_sampling.append("rwgt_ctZ_3_ctW_3_cpQM_3_cpQ3_3_cpt_3")

# //--------------------------------------------
# //--------------------------------------------

def Standalone_Validation(optsTrain, _list_lumiYears, _list_labels, _list_features, _list_processClasses, list_points_sampling, evalPoint, nEventsStandaloneVal):

    optsTrain["makeValPlotsOnly"] = True #Don't delete existing files
    _lumiName, _weightDir, _h5modelName, _ntuplesDir, _batchSize = Initialization_And_SanityChecks(optsTrain, _list_lumiYears, _list_processClasses, _list_labels)

    if optsTrain["parameterizedNN"] is False:
        print(colors.fg.red, "Error: strategy =", optsTrain["strategy"], ". Standalone validation not available for non-parameterized strategies (check validation plots produced by main training code)", colors.reset)
        return

    #-- Add option to control nof events to sample per hypothesis
    optsTrain["nEventsStandaloneVal"] = nEventsStandaloneVal
    #-- Add option to set point (WC values) at which DNN is to be evaluated
    optsTrain["evalPoint"] = np.squeeze(AddMissingOperatorsToValPointsNames(optsTrain, evalPoint) )

    #-- If want to validate over a single operator ( e.g. 'rwgt_ctZ_3') but the DNN was trained over more operators (e.g. 'rwgt_ctZ_3_ctW_0_cpQM_0_cpQ3_0_cpt_0', etc.), need to include missing operators into points names
    AddMissingOperatorsToValPointsNames(optsTrain, list_points_sampling)

    #-- Create output dir.
    standaloneValDir = _weightDir + 'StandaloneVal/'
    os.makedirs(standaloneValDir, exist_ok=True)

    #--- Load model
    tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    model = load_model(_h5modelName, compile=False) #compile=False <-> does not need to define any custom loss, since not needed for testing

    #-- Get data
    x_all=[]; y_all=[]; y_process_all=[]; PhysicalWeights_all=[]
    print(colors.fg.lightblue, "\n\n--- Get the data...\n", colors.reset)
    for idx, point in enumerate(list_points_sampling):
        print(colors.fg.lightblue, "=== POINT: ", point, " ===\n", colors.reset)
        x_tmp, y_tmp, y_process_tmp, PhysicalWeights_tmp, list_labels = Get_Data(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features, _weightDir, _ntuplesDir, _lumiName, singleThetaName=point)
        y_process_tmp = np.squeeze(y_process_tmp)

        #?? Keep that for multiclass ?
        y_process_tmp = np.zeros(len(y_process_tmp)) #1D
        y_process_tmp[:] = idx #Arbitrary identifier

        x_all.append(x_tmp); y_all.append(y_tmp); y_process_all.append(y_process_tmp); PhysicalWeights_all.append(PhysicalWeights_tmp)
        # print(x_tmp[:5]); print(y_tmp[:5]); print(y_process_tmp[:5]); print(model.predict(x_tmp[:5]))
    x=np.concatenate(x_all)
    y=np.concatenate(y_all)
    y_process=np.concatenate(y_process_all)
    PhysicalWeights=np.concatenate(PhysicalWeights_all)

    #-- Get model predictions
    predictions = np.squeeze(model.predict(x))
    # print(predictions)

    #-- Alter data for testing/debugging
    # x[y_process==0][:,-len(optsTrain["listOperatorsParam"]):] = 0.5

    # print(y[-10:])
    # print(y_process[-10:])
    # print(model.predict(x[-10:]))

#-- Create validation plots

    #For regressors
    Make_ScatterPlot_TrueVSPred(optsTrain, standaloneValDir, y, predictions, y_process, list_points_sampling)
    Make_Pull_Plot(optsTrain, standaloneValDir, y, predictions, list_points_sampling)

    #For classifiers
    Make_OvertrainingPlot_SinglePoints(optsTrain, standaloneValDir, list_labels, predictions, y_process, list_points_sampling)
    Make_ScatterPlot_2Dvars(optsTrain, _list_features, standaloneValDir, x, predictions, y_process, list_points_sampling)

    return

# //--------------------------------------------
# //--------------------------------------------

  ####   ####    ##   ##### ##### ###### #####  #####  #       ####  #####
 #      #    #  #  #    #     #   #      #    # #    # #      #    #   #
  ####  #      #    #   #     #   #####  #    # #    # #      #    #   #
      # #      ######   #     #   #      #####  #####  #      #    #   #
 #    # #    # #    #   #     #   #      #   #  #      #      #    #   #
  ####   ####  #    #   #     #   ###### #    # #      ######  ####    #

def Make_ScatterPlot_TrueVSPred(opts, standaloneValDir, truth, pred, procClass, list_points_sampling):
#See https://seaborn.ppred_data.org/tutorial/distributions.html
    """
    Make validation plots for regressor NN using test data. Compare predictions to true target values.
    """

    if opts["strategy"] not in ["ROLR", "RASCAL"]: return

    nodename='r'
    # if opts['regress_onLogr'] == True: nodename='log(r)'

    # xmin=-1.; xmax=20
    # ymin=-1.; ymax=20
    mycol = 'g'

    # print(truth.shape); print(pred.shape); print(procClass.shape)
    # print(truth[:10]); print(pred[:10]); print(procClass[:10])

# //--------------------------------------------
#DATA

    if opts["nofOutputNodes"] > 1:
        truth = truth[:,0]
        pred = pred[:,0]

    truth = np.squeeze(truth); pred = np.squeeze(pred); procClass = np.squeeze(procClass) #Make data 1D

    truth, pred, procClass = unison_shuffled_copies(truth, pred, procClass)
    # print(truth.shape); print(pred.shape)

    #Auto-adjust plot range
    quantileMin = np.quantile(truth, 0.10)
    quantileMax = np.quantile(truth, 0.90)
    xmin = math.floor(quantileMin) #Round down
    xmax = math.ceil(quantileMax) #Round up
    quantileMin = np.quantile(pred, 0.10)
    quantileMax = np.quantile(pred, 0.90)
    ymin = math.floor(quantileMin) #Round down
    ymax = math.ceil(quantileMax) #Round up

# //--------------------------------------------

    fig = plt.figure('splot', figsize=(10, 10))
    plt.title('Predicted VS True '+nodename+r'(x|$\theta_0,\theta_1$)')
    plt.xlabel(r'True '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string
    plt.ylabel(r'Learned '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string #color='darkorange'

    splot = sns.scatterplot(x=truth, y=pred, hue=procClass, palette='muted', s=15) #s <-> point size
    # splot = sns.scatterplot(x=truth[:,0], y=pred[:,0], hue=procClass[:], palette='muted')
    leg_handles = splot.get_legend_handles_labels()[0]
    splot.legend(leg_handles, list_points_sampling)
    ax = fig.gca()
    ax.set(xlim=(xmin, xmax))
    ax.set(ylim=(ymin, ymax))

    diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

    plotname = standaloneValDir + 'StandaloneScatterPlotLR_PredvsTruth.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved  LR scatter plot as :", colors.reset, plotname)
    fig.clear(); plt.close('splot')

# //--------------------------------------------
#Jointplot (scatterplot+density)

    #NB: 'jointplot' func creates its own figure (as do: FacetGrid, factorplot, lmplot, PairGrid, pairplot, JointGrid, jointplot.)
    jplot = sns.jointplot(x=truth, y=pred, color=mycol, alpha=0.5, xlim=(xmin,xmax), ylim=(ymin,ymax))
    plotname = standaloneValDir + 'JointPlot_PredvsTruth.png'
    jplot.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR joint plot as :", colors.reset, plotname)

    return

# //--------------------------------------------
# //--------------------------------------------

def Make_ScatterPlot_2Dvars(opts, list_features, standaloneValDir, x, pred, procClass, list_points_sampling):
#See https://seaborn.ppred_data.org/tutorial/distributions.html
    """
    2D scatterplot representing 2 input features in XY, and the DNN response as the z-axis colorbar.
    Should use this function to represent how the DNN response varies depending on some powerful variables, for events sampled according to a single scenario.
    """

    nEvents = 150 #Nof events to display

    # cmap = sns.cubehelix_palette(as_cmap=True)
    cmap = sns.color_palette('coolwarm', 7)

    var1 = "recoZ_Pt"
    var2 = "Mass_3l"

#Mass_3l, maxDelPhiLL, recoZ_Pt, mTW, recoTop_Eta, recoTop_Pt, ...

    idx1=-1; idx2=-1
    for i, feature in enumerate(list_features):
        if feature is var1: idx1=i
        elif feature is var2: idx2=i
    if idx1 is -1 or idx2 is -1: print('ERROR : feature not found !'); return

# //--------------------------------------------
#DATA

    x1 = x[:,idx1]
    x2 = x[:,idx2]

    x1, x2, pred, procClass = unison_shuffled_copies(x1, x2, pred, procClass)
    for array in [x1, x2, pred, procClass]:
        array = array[:nEvents]

    # x1, x2, pred, procClass = unison_shuffled_copies(x1, x2, pred, procClass)
    # print(x1.shape); print(x2.shape); print(pred.shape); print(procClass.shape)

    #Auto-adjust plot range
    # minElement = np.amin(x1)
    quantileMin = np.quantile(x1, 0.10)
    quantileMax = np.quantile(x1, 0.90)
    xmin = math.floor(quantileMin) #Round down
    xmax = math.ceil(quantileMax) #Round up
    quantileMin = np.quantile(x2, 0.10)
    quantileMax = np.quantile(x2, 0.90)
    ymin = math.floor(quantileMin) #Round down
    ymax = math.ceil(quantileMax) #Round up
    quantileMin = np.quantile(pred, 0.10)
    quantileMax = np.quantile(pred, 0.90)
    zmin = math.floor(quantileMin) #Round down
    zmax = math.ceil(quantileMax) #Round up

# //--------------------------------------------

    fs = 16 #fontsize
    pd = 20 #title padding
    fig = plt.figure('splot', figsize=(15, 10))
    # plt.tight_layout()
    plt.title('NN prediction VS ('+var1+','+var2+')', fontsize = 18)
    plt.xlabel(var1, fontsize=fs, labelpad=pd)
    plt.ylabel(var2, fontsize=fs, labelpad=pd)

    # splot = sns.scatterplot(x=x1, y=x2, hue=procClass, palette='muted')
    # splot = sns.scatterplot(x=x1, y=x2, hue=pred, s=50, cmap=cmap)
    splot = plt.scatter(x1, x2, c=pred, cmap="coolwarm", s=15) #s <-> point size
    # plt.colorbar(splot)
    cbar = plt.colorbar(splot)
    cbar.set_label('NN output value', fontsize=fs, labelpad=pd)
    # leg_handles = splot.get_legend_handles_labels()[0]
    # splot.legend(leg_handles, list_points_sampling)
    ax = fig.gca()
    ax.set(xlim=(xmin, xmax))
    ax.set(ylim=(ymin, ymax))
    plt.clim(zmin, zmax)
    ax.tick_params(labelsize=20)

    plotname = standaloneValDir + 'StandaloneScatterPlot2Dvar.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved  2Dvar scatter plot as :", colors.reset, plotname)
    fig.clear(); plt.close('splot')

    return

# //--------------------------------------------
# //--------------------------------------------

 #####  #    # #      #
 #    # #    # #      #
 #    # #    # #      #
 #####  #    # #      #
 #      #    # #      #
 #       ####  ###### ######

def Make_Pull_Plot(opts, standaloneValDir, truth, pred, list_points_sampling):
    """
    Plot difference between prediction and truth (error).
    """

    # if opts["strategy"] not in ["regressor", "ROLR", "RASCAL"]: return #Only useful for regressors

    #Transform classifier -> LR
    if opts["strategy"] is "classifier" or pts["strategy"] is "CARL":
        truth = r_from_s(truth)
        pred = r_from_s(pred)

    hpull = TH1F('Pull', 'Pull', 30, 0, 3); hpull.Sumw2(); hpull.SetDirectory(0)
    for idx in range(len(pred)):
        # tmp = pred[idx] - truth[idx]
        # tmp = (pred[idx] - truth[idx]) / truth[idx]
        # tmp = (pred[idx] - truth[idx]) / min(pred[idx], truth[idx])
        tmp = (pred[idx]) / truth[idx]
        hpull.Fill(tmp, 1.)
        # print('pred ', pred[idx], 'truth ', truth[idx], ' => ', tmp)

    hpull.SetFillColor(18)
    hpull.SetLineColor(1)

    c1 = ROOT.TCanvas()
    hpull.Draw("hist")
    plotname = standaloneValDir + "Pull_plot.png"
    c1.SaveAs(plotname)
    print(colors.fg.lightgrey, "\nSaved pull plot as :", colors.reset, plotname)

    return

# //--------------------------------------------
# //--------------------------------------------

  ####  #    # ###### #####  ##### #####    ##   # #    #
 #    # #    # #      #    #   #   #    #  #  #  # ##   #
 #    # #    # #####  #    #   #   #    # #    # # # #  #
 #    # #    # #      #####    #   #####  ###### # #  # #
 #    #  #  #  #      #   #    #   #   #  #    # # #   ##
  ####    ##   ###### #    #   #   #    # #    # # #    #

def Make_OvertrainingPlot_SinglePoints(opts, standaloneValDir, list_labels, predictions, y_process, list_points_sampling):

    nofOutputNodes = opts["nofOutputNodes"]

    #-- Define colors for all validation points
    #See: https://matplotlib.org/3.2.1/gallery/color/colormap_reference.html
    cols = cm.rainbow(np.linspace(0, 1, len(list_points_sampling))) #Rainbow
    # cols = cm.Oranges(np.linspace(0, 1, len(list_points_sampling))) #Nuances of orange
    # cols = cm.Pastel1(np.linspace(0, 1, len(list_points_sampling))) #Pastel
    # cols = cm.Set1(np.linspace(0, 1, len(list_points_sampling))) #Qualitative

    #-- Define legend names for all validation points
    legendNames = GetLegendNameEFTpoint(list_points_sampling)

    for inode in range(nofOutputNodes): #For each output node

        if opts["strategy"] in ["ROLR", "RASCAL"] and inode > 0: continue #Only for r node

        nbins = 20
        rmin = 0.; rmax = 1.

        fig = plt.figure('overtrain')
        timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
        timer.add_callback(close_event)

        ax = plt.axes()
        ax.set_xlim([rmin,rmax])

        #--- COSMETICS
        ax.patch.set_edgecolor('black')
        ax.patch.set_facecolor('#E6E6E6') #inner bkg color
        plt.grid(color='w', linestyle='solid') # draw solid white grid lines
        for spine in ax.spines.values():
            spine.set_visible(False)
            ax.xaxis.tick_bottom()
            ax.tick_params(colors='gray', direction='out')
            for tick in ax.get_xticklabels():
                tick.set_color('gray')
                for tick in ax.get_yticklabels():
                    tick.set_color('gray')


        for ipt, point in enumerate(list_points_sampling): #For each validation point

            #For each point, get corresponding color and legend name
            col = cols[ipt]
            if(point is "SM"): col = 'dimgrey'
            leg = legendNames[ipt]

            #-- Trick : for training histos, we want to compute the bin errors correctly ; to do this we first fill TH1Fs, then read their bin contents/errors
            # h = TH1F('h', '', nbins, rmin, rmax); h.Sumw2(); h.SetDirectory(0)
            # for pred in predictions:
            #     if opts["strategy"] in ["ROLR", "RASCAL"]: pred = 1/(pred+1) #Transform r -> s
            #     h.Fill(pred, 1.)
            # integ = h.Integral(0,h.GetNbinsX()+1)
            # if integ <= 0: integ = 1
            # sf_integral = abs(rmax - rmin) / nbins #h.Scale(1/integral) makes the sum of contents equal to 1, but does not account for the bin width
            # h.Scale(1./(integ*sf_integral))

            #Plot testing sig/bkg histos, normalized (no errors displayed <-> don't need TH1Fs)
            if nofOutputNodes == 1:
                if opts["strategy"] in ["ROLR", "RASCAL"]: tmp = 1./(predictions[y_process==ipt]+1) #Transform r -> s
                else: tmp = predictions[y_process==ipt]
            else:
                if opts["strategy"] in ["ROLR", "RASCAL"]: tmp = 1./(predictions[y_process==ipt][:,inode]+1) #Transform r -> s
                else: tmp = predictions[y_process==ipt][:,inode]

            if point is "SM":
                plt.hist(tmp, bins=nbins, range=(rmin,rmax), color=col, alpha=0.50, density=True, histtype='step', log=False, label=leg, edgecolor=col,fill=True)
            else:
                plt.hist(tmp, bins=nbins, range=(rmin,rmax), color=col, density=True, histtype='step', log=False, label=leg, edgecolor=col,fill=False, linewidth=2.5)

            # plt.hist(tmp, bins=nbins, range=(rmin,rmax), color=col, density=True, histtype='step', log=False, label=leg, edgecolor=col,fill=False, linewidth=2.5)

        myxlabel = "Classifier output"

        # plt.legend(loc='best', fontsize=14.5)
        plt.legend(fontsize=13, bbox_to_anchor=(1.1, 1.1), loc="upper right")
        # plt.legend(bbox_to_anchor=(1.1, 1.05))
        # plt.legend(loc='upper center', numpoints=1)
        # plt.title("Output distributions")
        plt.grid(axis='y', alpha=0.75)
        plt.grid(axis='x', alpha=0.75)
        plt.xlabel(myxlabel)
        plt.ylabel('PDF')

        if inode == 0:
            timer.start()
            plt.show()

        plotname = standaloneValDir + 'Overtraining_NN_' + list_labels[inode] + '.png'
        fig.savefig(plotname)
        # print("Saved Overtraining plot as : " + plotname)
        print(colors.fg.lightgrey, "\nSaved Overtraining plot as :", colors.reset, plotname)
        fig.clear()
        plt.close('overtrain')

    return

# //--------------------------------------------
# //--------------------------------------------


# //--------------------------------------------
# //--------------------------------------------

if __name__ == "__main__":

    Standalone_Validation(optsTrain,_list_lumiYears,_list_labels,_list_features,_list_processClasses, list_points_sampling, evalPoint, nEventsStandaloneVal)