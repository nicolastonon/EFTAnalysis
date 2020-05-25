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
from Utils.Helper import unison_shuffled_copies

# //--------------------------------------------
# //--------------------------------------------

list_points_val = []
# list_points_val.append("SM") #Keep
list_points_val.append("rwgt_ctZ_0.5")
# list_points_val.append("rwgt_ctZ_3")
list_points_val.append("rwgt_ctZ_5")

# //--------------------------------------------
# //--------------------------------------------

def Standalone_Validation(optsTrain, _list_lumiYears, _list_labels, _list_features, _list_processClasses, list_points_val):

    optsTrain["makeValPlotsOnly"] = True #Don't delete existing files
    _lumiName, _weightDir, _h5modelName, _ntuplesDir, _batchSize = Initialization_And_SanityChecks(optsTrain, _list_lumiYears, _list_processClasses, _list_labels)

    #-- Add option to control nof events to sample per hypothesis
    optsTrain["nEventsStandaloneVal"] = 1000

    #-- Create output dir.
    standaloneValDir = _weightDir + 'StandaloneVal/'
    os.makedirs(standaloneValDir, exist_ok=True)

    #--- Load model
    tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    model = load_model(_h5modelName)

    #-- Get data
    x_all=[]; y_all=[]; y_process_all=[]; PhysicalWeights_all=[]
    print(colors.fg.lightblue, "\n\n--- Get the data...\n", colors.reset)
    for idx, point in enumerate(list_points_val):
        x_tmp, y_tmp, y_process_tmp, PhysicalWeights_tmp, _ = Get_Data(optsTrain, _list_lumiYears, _list_processClasses, _list_labels, _list_features, _weightDir, _ntuplesDir, _lumiName, singleThetaName=point)
        y_process_tmp = np.squeeze(y_process_tmp); y_process_tmp[:] = idx
        x_all.append(x_tmp); y_all.append(y_tmp); y_process_all.append(y_process_tmp); PhysicalWeights_all.append(PhysicalWeights_tmp)
        # print(x_tmp[:10]); print(y_tmp[:10]); print(model.predict(x_tmp[:10]))
    x=np.concatenate(x_all); y=np.concatenate(y_all); y_process=np.concatenate(y_process_all); PhysicalWeights=np.concatenate(PhysicalWeights_all)

    #-- Get model predictions
    predictions = np.squeeze(model.predict(x))

    #-- Alter data for testing/debugging
    # x[y_process==0][:,-len(optsTrain["listOperatorsParam"]):] = 0.5

    # print(y[-10:])
    # print(y_process[-10:])
    # print(model.predict(x[-10:]))

    #-- Create validation plots
    Make_ScatterPlot_TrueVSPred(optsTrain, standaloneValDir, y, predictions, y_process)

    Make_Pull_Plot(optsTrain, standaloneValDir, y, predictions)

    return

# //--------------------------------------------
# //--------------------------------------------
#Scatterplot
#See https://seaborn.ppred_data.org/tutorial/distributions.html
def Make_ScatterPlot_TrueVSPred(opts, standaloneValDir, truth, pred, procClass):
    """
    Make validation plots for regressor NN using test data. Compare predictions to true target values.
    """

    if opts["strategy"] not in ["ROLR", "RASCAL"]: return

    nodename='r'
    # if opts['regress_onLogr'] == True: nodename='log(r)'

    xmin=-1.; xmax=20
    ymin=-1.; ymax=20
    mycol = 'g'

    # print(truth.shape); print(pred.shape); print(procClass.shape)
    # print(truth[:10]); print(pred[:10]); print(procClass[:10])

# //--------------------------------------------
#DATA

    if opts["nofOutputNodes"] > 1:
        truth = truth[:,0]
        pred = pred[:,0]

    truth, pred, procClass = unison_shuffled_copies(truth, pred, procClass)
    # print(truth.shape); print(pred.shape)

# //--------------------------------------------

    fig = plt.figure('splot', figsize=(10, 10))
    plt.title('Predicted VS True '+nodename+r'(x|$\theta_0,\theta_1$)')
    plt.xlabel(r'True '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string
    plt.ylabel(r'Learned '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string #color='darkorange'

    splot = sns.scatterplot(x=truth, y=pred, hue=procClass, palette='muted')
    # splot = sns.scatterplot(x=truth[:,0], y=pred[:,0], hue=procClass[:], palette='muted')
    leg_handles = splot.get_legend_handles_labels()[0]
    splot.legend(leg_handles, list_points_val)
    ax = fig.gca()
    ax.set(xlim=(xmin, xmax))
    ax.set(ylim=(ymin, ymax))

    diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

    plotname = standaloneValDir + 'StandaloneScatterPlotLR_PredvsTruth.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved  LR scatter plot as :", colors.reset, plotname)
    fig.clear(); plt.close('splot')

    return

# //--------------------------------------------
# //--------------------------------------------
def Make_Pull_Plot(opts, standaloneValDir, truth, pred):
    """
    Plot difference between prediction and truth (error).
    """

    hpull = TH1F('', '', 50, -1, 1); hpull.Sumw2(); hpull.SetDirectory(0)
    for idx in range(len(pred)):
        # tmp = pred[idx] - truth[idx]
        # tmp = (pred[idx] - truth[idx]) / truth[idx]
        tmp = (pred[idx] - truth[idx]) / min(pred[idx], truth[idx])
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



# //--------------------------------------------
# //--------------------------------------------

if __name__ == "__main__":

    Standalone_Validation(optsTrain,_list_lumiYears,_list_labels,_list_features,_list_processClasses, list_points_val)
