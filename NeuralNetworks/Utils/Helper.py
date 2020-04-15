# Nicolas Tonon (DESY)
# Python helper functions

import time   # time accounting
import ROOT
import numpy as np
from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
import tensorflow
import keras
import pandas as pd
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

    return

# //--------------------------------------------
# //--------------------------------------------
#-- Printout training infos
def batchOutput(batch, logs):

    print("Finished batch: " + str(batch))
    print(logs)

    return
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

    return

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

    return
# //--------------------------------------------
# //--------------------------------------------

#Normalize input features
def normalize(val, shift, scale):
    return (val-shift)/scale
# //--------------------------------------------
# //--------------------------------------------

#Alternative normalization #Separate train/test
'''
def normalize2(x_train, x_test):
    mu = np.mean(x_train, axis=0)
    std = np.std(x_train, axis=0)
    x_train_normalized = (x_train - mu) / std
    x_test_normalized = (x_test - mu) / std
    return x_train_normalized, x_test_normalized
'''
# //--------------------------------------------
# //--------------------------------------------

#Implementation from David's DeepPotato ; cf. Jonas' master thesis
#Using median and stddev from quantile is more robust against distributions with large tails
def get_normalization_iqr(np_array, q):
    """ Get shift and scale for events
    :param df: pandas DataFrame with events
    :param q: fraction of events that should be in the interval [-1, 1]
    :return: (shift, scale)
    """

    df = pd.DataFrame(data=np_array[0:,0:]) #Convert to panda DF

    q = (1 + q) / 2  # Tranform q so it works with quantile
    newDF = pd.DataFrame()
    for key in df.keys():
        newDF[key] = df[key].apply(lambda x: np.mean(x[np.nonzero(x)]) if hasattr(x, "__len__") else x)
    median = newDF.median()
    l = abs(newDF.quantile(1 - q) - median)
    r = abs(newDF.quantile(q) - median)
    for i, (il, ir) in enumerate(zip(l,r)):
        if il == ir:
            print(f"[WARNING] feature {df.keys()[i]} has no width")
            l[i] = 1.
            r[i] = 1.

    return median.values, np.maximum(l, r).values

# //--------------------------------------------
# //--------------------------------------------

#Printout the output of the first (=input) layer here for N events, e.g. to verify that the normalization layer works properly
def Printout_Outputs_FirstLayer(model, ilayer, xx):
    get_layer_output = keras.backend.function([model.layers[0].input], [model.layers[ilayer].output])
    layer_output = get_layer_output([xx])[0]
    print("\n", layer_output)

# //--------------------------------------------
# //--------------------------------------------
