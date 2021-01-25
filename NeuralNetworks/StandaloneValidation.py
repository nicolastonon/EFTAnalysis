#Perform validation independently from the main training code
#Useful to validate parameterized NN strategies. Sample 'new' data from given SM/EFT points, and evaluate at desired SM/EFT point.

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
from tensorflow.keras.utils import plot_model
from Utils.GetData import Get_Data
from Train_Neural_Network import optsTrain,_list_lumiYears,_list_labels,_list_features,_list_processClasses
from Utils.Validation_Control import *
from Utils.Predictions import *
from Utils.RegressorValidation import *
# from ann_visualizer.visualize import ann_viz
from pathlib import Path


# //--------------------------------------------
# //--------------------------------------------

use_xkcd_style = False #True <-> use pyplot's xkcd style for plotting

nEventsStandaloneVal = 50000 #Nof events to sample/display per point

#== CHOOSE SINGLE POINT AT WHICH TO EVALUATE EVENTS #NB: i.e. 'rwgt_ctW_3' corresponds to asking the NN 'are these events more EFT(ctW=3)-like, or more reference-like (<-> SM-like)'
#== NB: evalPoint=='' <-> evaluation point corresponds to the point to which each sample is drawn (<-> WC input values set accordingly)
#== NB: evalPoint=='' <-> ROC/... don't make sense (sig/bkg evaluated at different points)
evalPoint = ''
# evalPoint = "SM"
# evalPoint = "rwgt_ctz_1"
# evalPoint = "rwgt_ctz_3"
# evalPoint = "rwgt_ctw_1"
# evalPoint = "rwgt_ctw_2"
# evalPoint = "rwgt_ctw_3"
# evalPoint = "rwgt_ctw_5"
# evalPoint = "rwgt_cpqm_10"
# evalPoint = "rwgt_cpq3_5"
# evalPoint = "rwgt_cpt_15"
# evalPoint = "rwgt_ctz_5_ctw_5"
# evalPoint = "rwgt_ctZ_3_ctW_3_cpQM_3_cpQ3_3_cpt_3"

#== LIST OF POINTS FROM WHICH TO SAMPLE EVENTS  #NB: order of operators should be the same as used for training #NB: for CARL_multiclass, only 1 operator can be activated per point !
list_points_sampling = ["SM"] #Keep this !
# list_points_sampling.append("rwgt_ctz_1")
# list_points_sampling.append("rwgt_ctz_3")
list_points_sampling.append("rwgt_ctz_5")
# list_points_sampling.append("rwgt_ctw_1")
# list_points_sampling.append("rwgt_ctw_2")
# list_points_sampling.append("rwgt_ctw_3")
# list_points_sampling.append("rwgt_ctw_4")
list_points_sampling.append("rwgt_ctw_5")
# list_points_sampling.append("rwgt_cpqm_5")
# list_points_sampling.append("rwgt_cpq3_10")
# list_points_sampling.append("rwgt_cpt_15")
# list_points_sampling.append("rwgt_ctz_5_ctw_5")
# list_points_sampling.append("rwgt_ctW_2_cpQ3_4.5")

#== SCAN OPTIONS ==#
scan_singleOperator = False #True <-> plot output distributions for several values of a single operator
operator_scan = 'cpqm' #Operator to scan
range_step = [-15, 15, 3] #(range,steps) with which to scan operator
# range_step = [-1, 1, 2] #(range,steps) with which to scan operator
# range_step = [0, 15, 3] #(range,steps) with which to scan operator
only_SM_events = False #True <-> sample same SM events for each input WC value to test

# //--------------------------------------------
# //--------------------------------------------

def Standalone_Validation(optsTrain, _list_lumiYears, _list_labels, _list_features, _list_processClasses, list_points_sampling, evalPoint, nEventsStandaloneVal):

    optsTrain["makeValPlotsOnly"] = True #Don't delete existing files
    _lumiName, _weightDir, _h5modelName, _batchSize, _list_features = Initialization_And_SanityChecks(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features)

    if optsTrain["trainAtManyEFTpoints"] is False: #Trick: standalone val. code works with SMEFT samples; if applying it on classifier training, need to update lists; also change some options so that data is sampled properly
        optsTrain["trainAtManyEFTpoints"] = True;
        if optsTrain["strategy"] != "CARL_singlePoint":
            _list_processClasses = [["PrivMC_tZq"]]; _list_labels = ["PrivMC_tZq"]
            print(colors.fg.orange, "Warning: code [StandaloneValidation] only works for parameterized SMEFT samples. Setting [_list_processClasses = PrivMC_tZq] by default !", colors.reset)
            # print(colors.fg.red, "Error: strategy =", optsTrain["strategy"], ". Standalone validation not available for non-parameterized strategies (check validation plots produced by main training code)", colors.reset); return

    if scan_singleOperator:
        evalPoint = '' #Evaluate each sample at corresponding point (<-> set input WCs accordingly)
        ymax = -1 #Trick to keep same y-axis for each single plot
        idx_opScan = -1
        list_points_sampling, WCs = Get_ListPointsSampling_SingleOp(operator_scan, range_step)
        # print(list_points_sampling, WCs)
        for idx_op,opParam in enumerate(optsTrain["listOperatorsParam"]): #Find index of feature corresponding to the WC that we scan
            if operator_scan == opParam: idx_opScan = idx_op

    #-- If want to validate over a single operator ( e.g. 'rwgt_ctZ_3') but the DNN was trained over more operators (e.g. 'rwgt_ctZ_3_ctW_0_cpQM_0_cpQ3_0_cpt_0', etc.), need to include missing operators into points names
    AddMissingOperatorsToValPointsNames(optsTrain, list_points_sampling)

    optsTrain["nEventsStandaloneVal"] = nEventsStandaloneVal # Add option to control nof events to sample per hypothesis
    optsTrain["evalPoint"] = np.squeeze(AddMissingOperatorsToValPointsNames(optsTrain, evalPoint)) # Add option to set point (WC values) at which DNN is to be evaluated

    #-- Create output dir.
    standaloneValDir = _weightDir + 'StandaloneVal/'
    os.makedirs(standaloneValDir, exist_ok=True)

    #-- Load model
    tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    model = load_model(_h5modelName, compile=False) #compile=False <-> does not need to define any custom loss, since not needed for testing

    #-- Get data
    x_all=[]; y_all=[]; y_process_all=[]; PhysicalWeights_all=[]; pred_all=[]
    print(colors.fg.lightblue, "\n\n--- Get the data...\n", colors.reset)
    for idx, point in enumerate(list_points_sampling):
        print(colors.fg.lightblue, "\n=== POINT: ", point, " ===\n", colors.reset)

        # if scan_singleOperator:
        #     optsTrain["evalPoint"] = point
        #     print('evalPoint = ', point)

        if not scan_singleOperator or not only_SM_events or idx == 0: #-- Special case: may only need to reuse SM distribution at different EFT points (<-> only need to sample events at SM point)
            x_tmp, y_tmp, y_process_tmp, PhysicalWeights_tmp, list_labels, list_features, jointLR_allClasses, scores_allClasses_eachOperator = Get_Data(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features, _weightDir, ntuplesDir, _lumiName, singleThetaName=point)
            # list_labels = list_labels[::-1] #Trick: in main code, EFT is first (sig=1) and SM second (bkg=0); but in this code the first default sample is SM --> Reverse order
            y_process_tmp = np.squeeze(y_process_tmp)
            pred_tmp = np.squeeze(model.predict(x_tmp))

            # (Need to keep that for multiclass ?)
            y_process_tmp = np.zeros(len(y_process_tmp)) #1D
            y_process_tmp[:] = idx #Set to arbitrary identifier

            if optsTrain["strategy"] in ["RASCAL","CASCAL"]: pred_tmp = pred_tmp.T

            x_all.append(x_tmp); y_all.append(y_tmp); y_process_all.append(y_process_tmp); PhysicalWeights_all.append(PhysicalWeights_tmp); pred_all.append(pred_tmp)
            # print(x_tmp[:5]); print(y_tmp[:5]); print(y_process_tmp[:5]); print(model.predict(x_tmp[:5]))

        if scan_singleOperator:
            if point=='SM':
                x_SM=np.copy(x_tmp); y_process_SM=np.copy(y_process_tmp); PhysicalWeights_SM=np.copy(PhysicalWeights_tmp)
                if only_SM_events:
                    if optsTrain["trainAtManyEFTpoints"]: x_SM[:,x_SM.shape[1]-len(optsTrain["listOperatorsParam"])+idx_opScan] = 0 #Keep default SM distribution as reference
                    pred_SM = np.squeeze(model.predict(x_SM))
                # print(x_SM[0,:])
            else:
                if only_SM_events: #Case 1: keep 'default' SM distribution as fix reference --> don't modify 'xxx_SM' arrays ; only modify the 'xxx_tmp' arrays (copy SM distribution, change the WC values according to current point, and re-evaluate NN prediction)
                    x_tmp=np.copy(x_SM) #Start from SM events sampled previously
                    if optsTrain["trainAtManyEFTpoints"]: x_tmp[:,x_tmp.shape[1]-len(optsTrain["listOperatorsParam"])+idx_opScan] = WCs[idx] #Update WC values according to current EFT point
                    pred_tmp = np.squeeze(model.predict(x_tmp)) #Update NN predictions
                else: #Case 2: change both the SM and EFT distributions according to current EFT point --> already obtained the EFt data above, now need to modify the SM events (adapt WC values) and re-evalaute the NN predictions
                    if optsTrain["trainAtManyEFTpoints"]: x_SM[:,x_SM.shape[1]-len(optsTrain["listOperatorsParam"])+idx_opScan] = WCs[idx]
                    pred_SM=np.squeeze(model.predict(x_SM)) #Also need to set the WC input value for SM events according to current EFT point #Update prediction

                Store_TrainTestPrediction_Histograms(optsTrain, _lumiName, _list_features, ['EFT','SM'], [[pred_tmp,pred_SM]], [PhysicalWeights_tmp,PhysicalWeights_SM], [x_tmp,x_SM], [],[],[], True, operator_scan, str(WCs[idx]).replace('.0',''))
                ymax = Make_OvertrainingPlot_SinglePoints(optsTrain, standaloneValDir, list_labels, _list_features, np.concatenate((x_SM,x_tmp)), np.concatenate((pred_SM,pred_tmp)), np.concatenate((y_process_SM,y_process_tmp)), np.concatenate((PhysicalWeights_SM,PhysicalWeights_tmp)), ['SM',point], True, operator_scan, WCs, idx, ymax, weightDir=_weightDir)
                Make_ROCs(optsTrain, standaloneValDir, list_labels, _list_features, np.concatenate((x_SM,x_tmp)), np.concatenate((y_process_SM,y_process_tmp)), np.concatenate((pred_SM,pred_tmp)), np.concatenate((PhysicalWeights_SM,PhysicalWeights_tmp)), list_points_sampling, True, operator_scan, str(WCs[idx]).replace('.0',''), feature_name='recoZ_Pt') #recoZ_Pt/dEta_tjprime

    if scan_singleOperator:
        Make_Animation_fromParamOutputPlots(standaloneValDir, list_labels, list_points_sampling, operator_scan, WCs)
        Make_Animation_fromParamOutputPlots(standaloneValDir, list_labels, list_points_sampling, operator_scan, WCs, ROC=True)
        return

    x=np.concatenate(x_all)
    y=np.concatenate(y_all)
    y_process=np.concatenate(y_process_all)
    PhysicalWeights=np.concatenate(PhysicalWeights_all)
    predictions=np.concatenate(pred_all)

    #-- Alter data for testing/debugging
    # x[y_process==0][:,-len(optsTrain["listOperatorsParam"]):] = 0.5
    # print(y[-10:]); print(y_process[-10:]); print(model.predict(x[-10:]))

#-- Create validation plots

    #For regressors
    Make_ScatterPlot_TrueVSPred(optsTrain, standaloneValDir, y, predictions, y_process, list_points_sampling)
    Make_Pull_Plot(optsTrain, standaloneValDir, y, predictions, list_points_sampling)

    #For classifiers
    Make_OvertrainingPlot_SinglePoints(optsTrain, standaloneValDir, list_labels,_list_features, x, predictions, y_process, PhysicalWeights, list_points_sampling, weightDir=_weightDir)
    Make_OvertrainingPlot_SinglePoints(optsTrain, standaloneValDir, list_labels,_list_features, x, predictions, y_process, PhysicalWeights, list_points_sampling, feature_name="recoZ_Pt", weightDir=_weightDir)
    Make_ROCs(optsTrain, standaloneValDir, list_labels, _list_features, x, y, predictions, PhysicalWeights, list_points_sampling)
    Store_TrainTestPrediction_Histograms(optsTrain, _lumiName, _list_features, list_labels, [pred_all], PhysicalWeights_all, x_all)
    Make_Multiple_ROCs(optsTrain, standaloneValDir, list_labels, _list_features, x, y, predictions, PhysicalWeights, list_points_sampling, y_process)
    # Make_ScatterPlot_2Dvars(optsTrain, _list_features, standaloneValDir, x, predictions, y_process, PhysicalWeights, list_points_sampling, scores_allClasses_eachOperator)

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

    if opts["strategy"] is "ROLR":
        if xmax > 30: xmax = 30
        if ymax > 30: ymax = 30

# //--------------------------------------------

    fig = plt.figure('splot', figsize=(10, 10))
    plt.title('Predicted VS True '+nodename+r'(x|$\theta_0,\theta_1$)')
    plt.xlabel(r'True '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string
    plt.ylabel(r'Learned '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string #color='darkorange'

    splot = sns.scatterplot(x=truth, y=pred, hue=procClass, palette='muted', s=20) #s <-> point size
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

def Make_ScatterPlot_2Dvars(opts, list_features, standaloneValDir, x, pred, procClass, PhysicalWeights, list_points_sampling, scores_allClasses_eachOperator):
#See https://seaborn.ppred_data.org/tutorial/distributions.html
    """
    2D scatterplot representing 2 input features in XY, and the DNN response as the z-axis colorbar.
    Should use this function to represent how the DNN response varies depending on some powerful variables, for events sampled according to a single scenario.
    """

    nEvents = 150 #Nof events to display

    # cmap = sns.cubehelix_palette(as_cmap=True)
    cmap = sns.color_palette('coolwarm', 7)

    var1 = "recoZ_Pt"
    var2 = "recoZ_dPhill" #Mass_3l, recoZ_Eta,

    #Mass_3l, maxDelPhiLL, recoZ_Pt, mTW, recoTop_Eta, recoTop_Pt, ...

    if opts["strategy"] == "CARL_multiclass": return #Not supported anymore
    elif opts["strategy"] == "CASCAL": pred = pred[:,1] #First score component

    idx1=-1; idx2=-1
    for i, feature in enumerate(list_features):
        if feature is var1: idx1=i
        elif feature is var2: idx2=i
    if idx1 is -1 or idx2 is -1: print('ERROR : feature', var1, 'or ', var2, ' not found !'); return

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

    if opts["strategy"] not in ["regressor", "ROLR", "RASCAL"]: return #Only useful for regressors
    # if opts["strategy"] is "CARL_multiclass": return

    #Transform classifier -> LR
    # if opts["strategy"] is "classifier" or opts["strategy"] is "CARL":
    #     truth = r_from_s(truth)
    #     pred = r_from_s(pred)

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

def Make_OvertrainingPlot_SinglePoints(opts, standaloneValDir, list_labels, list_features, x, predictions, y_process, PhysicalWeights, list_points_sampling, scan=False, operator_scan='', WCs=[], idx_WC=-1, ymax=-1, feature_name="", weightDir=""):
    '''
    Plot output distributions for points in 'list_points_sampling'.

    scan: True <-> only consider a single EFT point at a time. Save separate plots.
    '''

    nofOutputNodes = opts["nofOutputNodes"]
    # print(list_points_sampling)

    if feature_name != "" and scan: return #Not compatible

    #-- Define colors for all validation points
    #See: https://matplotlib.org/3.2.1/gallery/color/colormap_reference.html
    cols = cm.rainbow(np.linspace(0, 1, len(list_points_sampling))) #Rainbow
    # cols = cm.Set1(np.linspace(0, 1, 10)) #Qualitative

    # cols = cm.viridis(np.linspace(0, 1, len(list_points_sampling))) #Perceptually uniform 1
    # cols = cm.plasma(np.linspace(0, 1, len(list_points_sampling))) #Perceptually uniform 2
    # cols = cm.cividis(np.linspace(0, 1, len(list_points_sampling))) #Perceptually uniform 3
    # cols = cm.Oranges(np.linspace(0, 1, len(list_points_sampling))) #Nuances of orange
    # cols = cm.Dark2(np.linspace(0, 1, len(list_points_sampling))) #Qualitative

    #-- Define legend names for all validation points
    legendNames = GetLegendNameEFTpoint(list_points_sampling)

    #-- Option: instead of NN prediction, may represent the distribution of some input feature
    for idx_feature, feature in enumerate(list_features):
        if feature_name == feature: predictions = x[:,idx_feature]

    for inode in range(nofOutputNodes): #For each output node

        # print('inode', inode)

        if opts["strategy"] in ["ROLR", "RASCAL", "CASCAL"] and inode > 0: continue #Only for r node

        nbins = 20
        # nbins = 70

        xrange = None #Hist x-range

        fig = plt.figure('overtrain')
        timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
        timer.add_callback(close_event)

        ax = plt.axes()
        if feature_name == "":
            # rmin = 0.3; rmax = 0.7; xrange = (rmin,rmax)
            rmin = 0.; rmax = 1.; xrange = (rmin,rmax)
            # ax.set_xlim([rmin,rmax])
            myxlabel = "Classifier output"

            #-- Read contents of the card we just produced #Set automatic bounds
            if Path(weightDir+"NN_info.txt").is_file():
                f = open(weightDir+"NN_info.txt", "r")
                contents = f.readlines()
                f.close()
                for line in contents:
                    if 'bounds' in line:
                        rmin = float(line.split()[1]); rmax = float(line.split()[2])

                        #-- Adjust to 0.05 below/above
                        rmin = int(rmin*100); rmin-= int(rmin%5); rmin/= 100.; rmin=0. if rmin < 0 else rmin
                        rmax = int(rmax*100)+5; rmax-= int(rmax%5); rmax/= 100.; rmax=1. if rmax > 1 else rmax

                        print('-- Set auto bounds: rmin = ', rmin, ' / rmax = ', rmax)

        else:
            if feature_name == "recoZ_Pt": #Hardcode xrange for easier comparisons b/w histos
                rmin = 0.; rmax = 750
                # ax.set_xlim([0,750])
            myxlabel = feature_name

        xrange = (rmin,rmax)

        #--- COSMETICS
        ax.patch.set_edgecolor('black')
        ax.patch.set_facecolor('#E6E6E6') #inner bkg color
        if not use_xkcd_style:
            plt.grid(color='w', linestyle='solid') # draw solid white grid lines
            plt.grid(axis='y', alpha=0.75)
            plt.grid(axis='x', alpha=0.75)

        for spine in ax.spines.values():
            spine.set_visible(False)
            ax.xaxis.tick_bottom()
            ax.tick_params(colors='gray', direction='out')
            for tick in ax.get_xticklabels():
                tick.set_color('gray')
                for tick in ax.get_yticklabels(): tick.set_color('gray')

        for ipt, point in enumerate(list_points_sampling): #For each validation point
            # print('Point: ', point)

            if point=='' or (scan==True and ipt>0 and point=='SM'): continue #Don't plot SM twice

            #For each point, get corresponding color and legend name
            col = cols[ipt]
            if(point is "SM"): col = 'dimgrey'
            leg = legendNames[ipt]

            #-- Plot normalized TEST sig/bkg histos (no errors displayed -> don't need TH1Fs)
            if nofOutputNodes == 1:
                if opts["strategy"] in ["ROLR", "RASCAL"] and feature_name == "":
                    tmp = 1./(predictions[y_process==ipt]+1) #Transform r -> s
                    weights_tmp = np.ones(len(tmp)) #Unweighted events
                else:
                    if scan and ipt>0:
                        tmp = predictions[y_process>0] #Scan single point at a time: no relation between ipt and y_process, only care about SM/EFT
                        weights_tmp = PhysicalWeights[y_process>0]
                    else:
                        tmp = predictions[y_process==ipt]
                        weights_tmp = PhysicalWeights[y_process==ipt]
            else:
                if opts["strategy"] in ["ROLR", "RASCAL"] and feature_name == "":
                    tmp = 1./(predictions[y_process==ipt][:,inode]+1) #Transform r -> s
                    weights_tmp = np.ones(len(tmp)) #Unweighted events
                elif opts["strategy"] == 'CASCAL' and feature_name == "":
                    # print('predictions[y_process==ipt].shape', predictions[y_process==ipt].shape)
                    tmp = predictions[y_process==ipt][:,inode]
                    weights_tmp = np.ones(len(tmp)) #Unweighted events
                else:
                    if scan and ipt>0:
                        tmp = predictions[y_process>0] #Scan single point at a time: no relation between ipt and y_process, only care about SM/EFT
                        weights_tmp = PhysicalWeights[y_process>0]
                    else:
                        tmp = predictions[y_process==ipt]
                        weights_tmp = PhysicalWeights[y_process==ipt]

            if point is "SM": plt.hist(tmp, bins=nbins, range=xrange, weights=weights_tmp, color=col, alpha=0.50, density=True, histtype='step', log=False, label=leg, edgecolor=col,fill=True)
            else: plt.hist(tmp, bins=nbins, range=xrange, weights=weights_tmp, color=col, density=True, histtype='step', log=False, label=leg, edgecolor=col,fill=False, linewidth=2.5)

        bottom, top = ax.get_ylim()
        ax.set_ylim([0., top*1.1])

        # plt.legend(loc='best', fontsize=14.5)
        plt.legend(fontsize=13, bbox_to_anchor=(1.1, 1.1), loc="upper right")
        # plt.legend(bbox_to_anchor=(1.1, 1.05))
        # plt.legend(loc='upper center', numpoints=1)
        # plt.title("Output distributions")
        plt.xlabel(myxlabel)
        plt.ylabel('PDF')

        if ymax != -1: plt.ylim([0.,ymax])
        elif scan: ymax = ax.get_ylim()[1] #Save ymax <-> keep same axis for all plots in scan

        # if inode == 0 and scan is False:
        #     timer.start()
        #     plt.show()

        plotname = standaloneValDir + 'Overtraining_NN_' + list_labels[inode]
        if scan: plotname+= '_' + operator_scan + '_' + str(WCs[idx_WC]).replace('.0','')
        elif feature_name != "": plotname+= "_" + feature_name
        plotname+= '.png'
        fig.savefig(plotname)
        print(colors.fg.lightgrey, "\nSaved Overtraining plot as :", colors.reset, plotname)
        fig.clear()
        plt.close('overtrain')

    return ymax

# //--------------------------------------------
# //--------------------------------------------

 #####   ####   ####   ####
 #    # #    # #    # #
 #    # #    # #       ####
 #####  #    # #           #
 #   #  #    # #    # #    #
 #    #  ####   ####   ####

#-- Make single '1vs all' ROC curve plot
def Make_ROCs(opts, standaloneValDir, list_labels, list_features, x, truth, predictions, PhysicalWeights, list_points_sampling, scan=False, op='', WC='', feature_name=''):

    if "CARL" not in opts["strategy"]: return
    if "sm" not in list_points_sampling and "SM" not in list_points_sampling: return #Compare EFT to SM
    lw = 2 #linewidth

    fig = plt.figure('roc')

    truth_tmp = truth[:]
    truth_tmp[truth_tmp>0] = 1 #Positive integers --> Set to 0 or 1

    fpr, tpr, _ = roc_curve(truth, predictions)
    roc_auc = auc(fpr, tpr)
    plt.plot(1-fpr, tpr, color='cornflowerblue', lw=lw, label='ROC NN (test) (AUC = {1:0.2f})' ''.format(0, roc_auc))

    #-- Option: instead of NN prediction, make ROC for some input feature
    if feature_name != '':
        idx_feature = -1
        for idx_feature, feature in enumerate(list_features):
            if feature_name == feature: predictions = x[:,idx_feature]
        if idx_feature != -1:
            fpr, tpr, _ = roc_curve(truth, predictions)
            roc_auc = auc(fpr, tpr)
            plt.plot(1-fpr, tpr, color='darkorange', lw=lw, label='ROC '+ feature_name +' (test) (AUC = {1:0.2f})' ''.format(0, roc_auc))


    ax = fig.gca()
    ax.set_xticks(np.arange(0, 1, 0.1))
    ax.set_yticks(np.arange(0, 1., 0.1))
    plt.grid()
    plt.plot([1, 0], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('Signal efficiency')
    plt.ylabel('Background rejection')
    if scan:
        plt.title(op + '=' + WC + ' vs SM', fontsize=20)
        plt.legend(loc='lower left')
    else:
        plt.title('')
        plt.legend(loc='best')

    plotname = standaloneValDir + 'ROC'
    if scan: plotname+= "_"+op+'_'+WC
    plotname+= ".png"
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved ROC plot as :", colors.reset, plotname)
    fig.clear()
    plt.close('roc')

    return


#-- Make ROC for 1 vs 1 (as many ROC curves as competing hypotheses)
def Make_Multiple_ROCs(opts, standaloneValDir, list_labels, list_features, x, truth, predictions, PhysicalWeights, list_points_sampling, y_process, scan=False, op='', WC='', feature_name=''):

    if "CARL" not in opts["strategy"]: return
    if "sm" not in list_points_sampling and "SM" not in list_points_sampling: return #Compare EFT to SM
    lw = 2 #linewidth

    fig = plt.figure('multiroc')

    nofOutputNodes = opts["nofOutputNodes"]

    truth_tmp = truth[:]
    truth_tmp[truth_tmp>0] = 1 #Positive integers --> Set to 0 or 1

    legendNames = GetLegendNameEFTpoint(list_points_sampling)

    #-- Define colors for all validation points
    #See: https://matplotlib.org/3.2.1/gallery/color/colormap_reference.html
    cols = cm.rainbow(np.linspace(0, 1, len(list_points_sampling))) #Rainbow

    for ipt, point in enumerate(list_points_sampling): #For each validation point
        # print('Point: ', point)

        if ipt == 0: continue #Don't plot SM vs SM...!

        #For each point, get corresponding color and legend name
        col = cols[ipt]
        if(point is "SM"): col = 'dimgrey'
        leg = legendNames[ipt]

        fpr, tpr, _ = roc_curve(np.concatenate((truth[y_process==0],truth[y_process==ipt])), np.concatenate((predictions[y_process==0],predictions[y_process==ipt])))
        roc_auc = auc(fpr, tpr)
        plt.plot(1-fpr, tpr, color=col, lw=lw, label=leg+' (AUC = {1:0.2f})' ''.format(0, roc_auc))

        #-- Option: instead of NN prediction, make ROC for some input feature
        if feature_name != '':
            idx_feature = -1
            for idx_feature, feature in enumerate(list_features):
                if feature_name == feature: predictions = x[:,idx_feature]
            if idx_feature != -1:
                fpr, tpr, _ = roc_curve(truth, predictions)
                roc_auc = auc(fpr, tpr)
                plt.plot(1-fpr, tpr, color='darkorange', lw=lw, label='ROC '+ feature_name +' (test) (AUC = {1:0.2f})' ''.format(0, roc_auc))


    ax = fig.gca()
    ax.set_xticks(np.arange(0, 1, 0.1))
    ax.set_yticks(np.arange(0, 1., 0.1))
    plt.grid()
    plt.plot([1, 0], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('Signal efficiency')
    plt.ylabel('Background rejection')
    if scan:
        plt.title(op + '=' + WC + ' vs SM', fontsize=20)
        plt.legend(loc='lower left')
    else:
        plt.title('')
        plt.legend(loc='best')

    plotname = standaloneValDir + 'MultiROC'
    if scan: plotname+= "_"+op+'_'+WC
    plotname+= ".png"
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved multiROC plot as :", colors.reset, plotname)
    fig.clear()
    plt.close('multiroc')

    return

# //--------------------------------------------
# //--------------------------------------------


# //--------------------------------------------
# //--------------------------------------------

if __name__ == "__main__":

    if use_xkcd_style: plt.xkcd() # XKCD-style plotting

    Standalone_Validation(optsTrain,_list_lumiYears,_list_labels,_list_features,_list_processClasses, list_points_sampling, evalPoint, nEventsStandaloneVal)
