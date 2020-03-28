# Nicolas Tonon (DESY)
# Python helper functions

import time   # time accounting
import ROOT
import numpy as np
from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
import tensorflow
import keras
from matplotlib import pyplot as plt
from Utils.ColoredPrintout import colors

# //--------------------------------------------
##     ## ######## ##       ########  ######## ########
##     ## ##       ##       ##     ## ##       ##     ##
##     ## ##       ##       ##     ## ##       ##     ##
######### ######   ##       ########  ######   ########
##     ## ##       ##       ##        ##       ##   ##
##     ## ##       ##       ##        ##       ##    ##
##     ## ######## ######## ##        ######## ##     ##
# //--------------------------------------------

#-- Automatically close matplotlib plot after some time
def close_event():
    plt.close() #timer calls this function after 3 seconds and closes the window

# //--------------------------------------------
# //--------------------------------------------
#-- Printout training infos
def batchOutput(batch, logs):

    print("Finished batch: " + str(batch))
    print(logs)
# //--------------------------------------------
# //--------------------------------------------

#-- Write DNN input variables to a .txt file
def Write_Variables_To_TextFile(weight_dir, var_list):
    text_file = open(weight_dir + "ListVariables.txt", "w")
    for var in var_list:
        text_file.write(var)
        text_file.write("\n")
    text_file.close()
    # print("\n===> Saved list of variables in : " + weight_dir + "ListVariables.txt\n\n")
    print(colors.fg.lightgrey, '===> Saved list of variables in : ' + weight_dir + 'ListVariables.txt', colors.reset)

# //--------------------------------------------
# //--------------------------------------------
#-- Get execution time
class TimeHistory(tensorflow.keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.times = []

    def on_epoch_begin(self, batch, logs={}):
        self.epoch_time_start = time.time()

    def on_epoch_end(self, batch, logs={}):
        self.times.append(time.time() - self.epoch_time_start)

# //--------------------------------------------
# //--------------------------------------------
#-- Get name corresponding to the data-taking years which are considered in the DNN training
def Get_LumiName(lumi_years):

    # Set a unique name to each combination of year(s)
    if len(lumi_years) == 1:
        lumiName = lumi_years[0]
    elif len(lumi_years) == 2:
        if lumi_years[0] == "2016" and lumi_years[1] == "2017":
            lumiName = "201617"
        elif lumi_years[0] == "2016" and lumi_years[1] == "2018":
            lumiName = "201618"
        elif lumi_years[0] == "2017" and lumi_years[1] == "2018":
            lumiName = "201718"
    elif len(lumi_years) == 3:
        lumiName = "Run2"
    else:
        print(colors.bold, colors.bg.red, 'ERROR : wrong lumi_years values !', colors.reset)
        exit(1)

    return lumiName


# //--------------------------------------------
# //--------------------------------------------

#Sanity checks of input args
def SanityChecks_Parameters(processClasses_list, labels_list):

    if len(processClasses_list) == 0:
        print(colors.fg.red, 'ERROR : no process class defined...', colors.reset); exit(1)
    elif len(processClasses_list) is not len(labels_list):
        print(colors.fg.red, 'ERROR : sizes of lists processClasses_list and labels_list are different...', colors.reset); exit(1)
# //--------------------------------------------
# //--------------------------------------------

#Normalize input features
def normalize(val, m, dev):
    return (val-m)/dev
# //--------------------------------------------
# //--------------------------------------------

def normalize2(x_train, x_test):
    mu = np.mean(x_train, axis=0)
    std = np.std(x_train, axis=0)
    x_train_normalized = (x_train - mu) / std
    x_test_normalized = (x_test - mu) / std
    return x_train_normalized, x_test_normalized

# //--------------------------------------------
# //--------------------------------------------
