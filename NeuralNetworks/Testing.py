
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
 #######  ########  ######## ####  #######  ##    ##  ######
##     ## ##     ##    ##     ##  ##     ## ###   ## ##    ##
##     ## ##     ##    ##     ##  ##     ## ####  ## ##
##     ## ########     ##     ##  ##     ## ## ## ##  ######
##     ## ##           ##     ##  ##     ## ##  ####       ##
##     ## ##           ##     ##  ##     ## ##   ### ##    ##
 #######  ##           ##    ####  #######  ##    ##  ######
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

# Analysis options
# //--------------------------------------------
# -- Choose here what data you want to consider (separate ntuples per year) ; same convention as for main analysis code
# Naming convention enforced : 2016+2017 <-> "201617" ; etc.; 2016+2017+2018 <-> "Run2" # NB : years must be placed in the right order !
_lumi_years = []
# _lumi_years.append("2016")
_lumi_years.append("2017")
# _lumi_years.append("2018")

#Signal process must be first
_processClasses_list = [
                # ["PrivMC_tZq"],
                ["tZq"],
                ["ttZ"]]
                # ["ttZ"], ["ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep"]]
                # ["ttZ", "ttW", "ttH", "WZ", "ZZ4l", "TTbar_DiLep",]]

_labels_list =  ["tZq",
                # "ttZ", "Backgrounds"]
                "Backgrounds"]

cuts = "passedBJets==1" #Event selection, both for train/test ; "1" <-> no cut
# //--------------------------------------------

#--- Training options
# //--------------------------------------------
_nepochs = 10 #Number of training epochs (<-> nof times the full training dataset is shown to the NN)
_batchSize = 1000 #Batch size (<-> nof events fed to the network before its parameter get updated)

_maxEvents_perClass = 100 #max nof events to be used for each process ; -1 <-> all events
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
from Utils.GetData import Get_Data_For_DNN_Training
from Utils.Optimizer import Get_Loss_Optim_Metrics
from Utils.ColoredPrintout import colors
from Utils.Output_Plots_Histos import Create_TrainTest_ROC_Histos, Create_Control_Plots

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

        x_train, y_train, x_test, y_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, x, y, PhysicalWeights_allClasses, LearningWeights_allClasses, means, stddev, x_control_firstNEvents = Get_Data_For_DNN_Training(weight_dir, _lumi_years, _ntuples_dir, _processClasses_list, _labels_list, var_list, cuts, _nof_output_nodes, _maxEvents_perClass, _splitTrainEventFrac, _nEventsTot_train, _nEventsTot_test, lumiName)

        # print(x.shape); print(x)

        df = pd.DataFrame(data=x[0:,0:], columns=var_list[:]) #Should use separate list for varnames
        # print(df)

        corr = df.corr()

        fig, ax = plt.subplots()

        mask = np.triu(np.ones_like(corr, dtype=np.bool))
        # mask = np.tril(np.ones_like(corr, dtype=np.bool))

        # corr = abs(corr)

        # palette = sns.diverging_palette(240, 10, n=9)
        # palette = sns.diverging_palette(20, 220, n=256)
        palette = 'coolwarm'

        # Draw the heatmap -- see : https://seaborn.pydata.org/generated/seaborn.heatmap.html
        # hm = sns.heatmap(corr, mask=mask, cmap=palette, vmin=-1, vmax=1, center=0, square=True, linewidths=1., linecolor='black', cbar_kws={"shrink": .5}, annot = True, fmt='.1g', )
        hm = sns.heatmap(corr, mask=mask, cmap=palette, vmin=0, vmax=1, center=0, square=True, linewidths=0.5, linecolor='black')
        hm.set_xticklabels(hm.get_xticklabels(), rotation=45, horizontalalignment='right')

        plt.show()

        # df.hist(figsize = (12,10))


# //--------------------------------------------
# //--------------------------------------------

TEST()
