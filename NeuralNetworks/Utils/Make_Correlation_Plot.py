# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import argparse
# //--------------------------------------------
import tensorflow
import keras
import numpy as np
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score

import ROOT
from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
import numpy as np
from numpy import random
from root_numpy import root2array, tree2array, array2root, hist2array
import os

import math
import json
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def Make_Correlation_Plot():
    '''
    Make a correlation plot of input features.
    '''

    list_features = []
    list_features.append('Mass_3l')
    list_features.append('lAsymmetry')
    list_features.append('maxDeepJet')
    list_features.append('nbjets')
    list_features.append('recoZ_Eta')
    list_features.append('recoZ_dPhill')
    list_features.append('dR_lWjprime')
    list_features.append('jprime_Pt')
    list_features.append('metEt')
    list_features.append('mbjMax')

    list_features.append('mTW')
    list_features.append('jPrimeAbsEta')
    list_features.append('dR_blW')
    list_features.append('recoLepTop_Pt')

    # list_features.append('lep1_pt')
    # list_features.append('lep2_pt')
    # list_features.append('lep3_pt')
    # list_features.append('lep1_eta')
    # list_features.append('lep2_eta')
    # list_features.append('lep3_eta')
    # list_features.append('lep1_phi')
    # list_features.append('lep2_phi')
    # list_features.append('lep3_phi')
    # list_features.append('jet1_pt')
    # list_features.append('jet2_pt')
    # list_features.append('jet3_pt')
    # list_features.append('jet1_eta')
    # list_features.append('jet2_eta')
    # list_features.append('jet3_eta')
    # list_features.append('jet1_phi')
    # list_features.append('jet2_phi')
    # list_features.append('jet3_phi')
    # list_features.append('jet1_DeepJet')
    # list_features.append('jet2_DeepJet')
    # list_features.append('jet3_DeepJet')

    # list_features = []
    # list_features.append('mHT')
    # list_features.append('mTW')
    # list_features.append('Mass_3l')
    # list_features.append('recoZ_Pt')
    # list_features.append('recoZ_Eta')
    # list_features.append('recoZ_dPhill')
    # list_features.append('dR_blW')
    # list_features.append('dR_tZ')
    # list_features.append('recoLepTop_Pt')
    # list_features.append('dEta_Zjprime')
    # list_features.append('jPrimeAbsEta')

    filepath = '../../input_ntuples/DATA.root'
    file = TFile.Open(filepath)
    # tree = file.Get('result')
    tree = file.Get('result')

    maxEvents=-1 #Don't need all events
    mask = np.ones(len(list_features), np.bool) #Default

    x = tree2array(tree, branches=list_features)[:maxEvents]
    x = np.column_stack([x[name] for name in x.dtype.names]) #1D --> 2D
    x = x.astype(np.float32) #Convert all to floats
    # print(x)

    #Can avoid plotting theory parameter inputs (by default), and p4 variables
    list_features = np.array(list_features)

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

    outname = './CorrelMatrix.png'
    plt.savefig(outname)
    print("\nSaved correlation matrix plot as :", outname)
    plt.close('all')

    return


# //--------------------------------------------
if __name__ == "__main__":

    Make_Correlation_Plot()
