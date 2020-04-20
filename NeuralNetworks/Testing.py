# Analysis options
# //--------------------------------------------
# -- Choose here what data you want to consider (separate ntuples per year) ; same convention as for main analysis code
# Naming convention enforced : 2016+2017 <-> "201617" ; etc.; 2016+2017+2018 <-> "Run2" # NB : years must be placed in the right order !
_lumi_years = []
# _lumi_years.append("2016")
_lumi_years.append("2017")
# _lumi_years.append("2018")

#Signal process must be first
_processClasses_list = []
_processClasses_list.append(["tZq"])
_processClasses_list.append(["PrivMC_tZq"])
# _processClasses_list.append(["PrivMC_tZq_ctz"])
# _processClasses_list.append(["PrivMC_tZq_ctw"])
# _processClasses_list.append(["ttZ"])
# _processClasses_list.append(["ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep"])
# _processClasses_list.append(["ttZ", "ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep",])

_labels_list = []
_labels_list.append("tZq")
_labels_list.append("PrivMC_tZq")
# _labels_list.append("PrivMC_tZq_ctz")
# _labels_list.append("PrivMC_tZq_ctw")
# _labels_list.append("ttZ")
# _labels_list.append("Backgrounds")

cuts = "passedBJets==1" #Event selection, both for train/test ; "1" <-> no cut
# //--------------------------------------------

#--- Training options
# //--------------------------------------------
_regress = False #True <-> DNN used for regression ; False <-> classification

_nepochs = 5 #Number of training epochs (<-> nof times the full training dataset is shown to the NN)
_batchSize = 500 #Batch size (<-> nof events fed to the network before its parameter get updated)

_maxEvents_perClass = 10 #max nof events to be used for each process ; -1 <-> all events
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
import pandas as pd
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score
from tensorflow.keras.models import load_model

from Utils.FreezeSession import freeze_session
from Utils.Helper import batchOutput, Write_Variables_To_TextFile, TimeHistory, Get_LumiName, SanityChecks_Parameters, Printout_Outputs_FirstLayer
from Utils.Model import Create_Model
from Utils.Callbacks import Get_Callbacks
from Utils.GetData import Get_Data
from Utils.Optimizer import Get_Loss_Optim_Metrics
from Utils.ColoredPrintout import colors
from Utils.Output_Plots_Histos import Create_TrainTest_ROC_Histos, Create_Control_Plots

from Utils.EFTUtilities import *

import matplotlib.pyplot as plt
import seaborn as sns

# //--------------------------------------------
# //--------------------------------------------

def TEST():
        print('\n\n')
        print(colors.bg.orange, colors.bold, "=====================================", colors.reset)
        print('\t', colors.fg.orange, colors.bold, "...TESTING...", colors.reset)
        print(colors.bg.orange, colors.bold, "=====================================", colors.reset, '\n\n')

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

        #Get data
        print(colors.fg.lightblue, "--- Read and shape the data...", colors.reset); print('\n')
        _transfType = 'quantile' #Feature norm. method -- 'range', 'gauss', 'quantile'
        x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, EFTweights_train, EFTweights_test, EFTweightIDs_train, EFTweightIDs_test, EFT_FitCoeffs_train, EFT_FitCoeffs_test, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, EFTweights_allClasses, EFTweightIDs_allClasses, EFT_FitCoeffs_allClasses, shifts, scales, x_control_firstNEvents, xTrainRescaled = Get_Data(_regress, weight_dir, _lumi_years, _ntuples_dir, _processClasses_list, _labels_list, var_list, cuts, _nof_output_nodes, _maxEvents_perClass, _splitTrainEventFrac, _nEventsTot_train, _nEventsTot_test, lumiName, _transfType)

# //--------------------------------------------

    test_points = ['rwgt_ctz_3p15_ctw_min4p3_cpqm_4p62_cpq3_5p83_cpt_min3p48', 'rwgt_ctz_4p15_ctw_min2p3_cpqm_1p62_cpq3_0p83_cpt_min3p48']
    operatorNames, operatorWCs = Parse_EFTpoint_IDs(test_points) #Get the lists of operator names and WC values for this process #NB: assumes that they are identical for all events in this process
    n_components, components = Find_Components(operatorNames) #Determine the components required to parameterize the event weight #NB: assumes that they are identical for all events in this process
    effWC_components = Get_EffectiveWC_eachComponents(n_components, components, operatorWCs) #Determine the 'effective WC' values associated with each component, for each benchmark point
    # fit_coeffs = Get_FitCoefficients(effWC_components, benchmark_weights=list_EFTweights_allClasses[iclass]) #Determine the fit coefficients of the events, based on the benchmark weights and 'effective WC' values

    newWeights = Get_Extrapolated_EFTweight(effWC_components, EFT_FitCoeffs_allClasses)

    print(newWeights)

# //--------------------------------------------
# # //--------------------------------------------
# def xxx():


# //--------------------------------------------
# //--------------------------------------------

TEST()

# str='rwgt_ctZ_3.27_ctW_3.61_cpQM_-4.04_cpQ3_-1.85_cpt_-0.52'
# str=['rwgt_ctZ_5_ctW_3']
# str=['rwgt_ctZ_5_ctW_3', 'rwgt_ctZ_3.27_ctW_3.61']
# list_op, list_WC = Parse_EFTpoint_IDs(str)

# n_components, components = Find_Components(list_op)
# effWC_eachComponent = Get_EffectiveWC_eachComponents(n_components, components, list_WC)
# effWC_eachComponent = np.random.random_sample((6,6))
# benchmark_weights = np.array([0.5, 1, 2, 3, 4, 5])
# fit_coeffs = Get_FitCoefficients(effWC_eachComponent, benchmark_weights)
# str='rwgt_ctZ_5_ctW_7'
# w = Get_Extrapolated_EFTweight(n_components, components, fit_coeffs, str)
