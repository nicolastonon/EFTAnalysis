# Nicolas Tonon (DESY)
# Python helper functions

import time   # time accounting
import ROOT
import numpy as np
from numpy import random
from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
import tensorflow
import keras
import pandas as pd
import re
import os
from scipy.stats import ks_2samp, anderson_ksamp, chisquare
from matplotlib import pyplot as plt
from Utils.ColoredPrintout import colors

# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
 ######   ######## ##    ## ######## ########     ###    ##
##    ##  ##       ###   ## ##       ##     ##   ## ##   ##
##        ##       ####  ## ##       ##     ##  ##   ##  ##
##   #### ######   ## ## ## ######   ########  ##     ## ##
##    ##  ##       ##  #### ##       ##   ##   ######### ##
##    ##  ##       ##   ### ##       ##    ##  ##     ## ##
 ######   ######## ##    ## ######## ##     ## ##     ## ########
# //--------------------------------------------
# //--------------------------------------------
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
    """
    Get shift and scale for events
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
    # print('median', median)
    for i, (il, ir) in enumerate(zip(l,r)):
        if il == ir:
            print(f"[WARNING] feature {df.keys()[i]} has no width --> Set width = ", max(il, ir), ' (max. value)')
            # print('il', il); print('ir', ir)
            # l[i] = 1.; r[i] = 1.
            l[i] = il; r[i] = ir #CHANGED -- happens e.g. for discrete variables perfectly centered at 0. Better to return the max value, so that all values will effectively lie in [-1;+1]

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

#Shuffles coherently the rows of N arrays of same length
def unison_shuffled_copies(*arr):
    assert all(len(a) for a in arr)
    p = np.random.permutation(len(arr[0]))
    return (a[p] for a in arr)

#-- Alternative example :
# indices = np.arange(list_x_arrays_allClasses[i].shape[0])
# np.random.shuffle(indices) #Get shuffled indices
# list_x_arrays_allClasses[i] = list_x_arrays_allClasses[i][indices]

# //--------------------------------------------
# //--------------------------------------------

# See : https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.ks_2samp.html
# Computes the Kolmogorov-Smirnov statistic on 2 samples. This is a two-sided test for the null hypothesis that 2 independent samples are drawn from the same continuous distribution
# Returns (KS,pval). If the K-S statistic is small or the p-value is high, then we cannot reject the hypothesis that the distributions of the two samples are the same
def KS_test(values1, values2):

    # calculate the significance
    value, pvalue = ks_2samp(values1, values2)

    print('=== K-S test : KS statistics =', float('%.4g' % value), ' / p-value =', pvalue)
    if pvalue > 0.05:
    	print('Samples are likely drawn from the same distributions (fail to reject H0)')
    else:
    	print('Samples are likely drawn from different distributions (reject H0)')
    print('\n')

    return

# //--------------------------------------------
# //--------------------------------------------

#See : https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.anderson_ksamp.html
#The k-sample Anderson-Darling test is a modification of the one-sample Anderson-Darling test. It tests the null hypothesis that k-samples are drawn from the same population without having to specify the distribution function of that population. The critical values depend on the number of samples.
#stat=Normalized k-sample Anderson-Darling test statistic ; critical_values=The critical values for significance levels 25%, 10%, 5%, 2.5%, 1%, 0.5%, 0.1% ; significance_level=An approximate significance level at which the null hypothesis for the provided samples can be rejected. The value is floored / capped at 0.1% / 25%.
def Anderson_Darling_test(sample1, sample2):

    stat, critical_values, significance_level = anderson_ksamp(samples=[sample1,sample2])

    print('== A-D statistics =', stat)
    # print('critical_values =', critical_values)
    print('Approx. signif. level for rejecting hypothesis that both samples are drawn from same distrib. =', significance_level)

    return

# //--------------------------------------------
# //--------------------------------------------

def ChiSquare_test(obs, exp):

    chi2, pval = chisquare(obs, exp)
    print('== Chi-2 =', chi2)
    print('p-value =', pval)

    return

# //--------------------------------------------
# //--------------------------------------------

def my_training_batch_generator(features, labels, batch_size): # Create empty arrays to contain batch of features and labels

    batch_features = np.zeros((batch_size, features.shape[1]))
    batch_labels = np.zeros((batch_size, labels.shape[1]))

    while True:
        # choose random index in features
        for i in range(batch_size):
            index = random.choice(len(features),1)
            batch_features[i] = features[index]
            batch_labels[i] = labels[index]

        yield batch_features, batch_labels #'Yield' keyword returns a generator. It suspends the function, which resumes at next call ; allows to produce a series of objects over time, instead of computing/returning everything at once

# //--------------------------------------------
# //--------------------------------------------






# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
 ######  ########  ########  ######  #### ######## ####  ######
##    ## ##     ## ##       ##    ##  ##  ##        ##  ##    ##
##       ##     ## ##       ##        ##  ##        ##  ##
 ######  ########  ######   ##        ##  ######    ##  ##
      ## ##        ##       ##        ##  ##        ##  ##
##    ## ##        ##       ##    ##  ##  ##        ##  ##    ##
 ######  ##        ########  ######  #### ##       ####  ######
# //--------------------------------------------
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
def Initialization_And_SanityChecks(opts, lumi_years, processClasses_list, labels_list):

# //--------------------------------------------

#-- Sanity checks
    if len(processClasses_list) == 0:
        print(colors.fg.red, 'ERROR : no process class defined...', colors.reset); exit(1)
    elif len(processClasses_list) is not len(labels_list):
        print(colors.fg.red, 'ERROR : sizes of lists processClasses_list and labels_list are different...', colors.reset); exit(1)

    for procClass, label in zip(processClasses_list, labels_list):
        if "PrivMC" in label and len(procClass) > 1:
            print(colors.fg.red, 'ERROR : process classes containing a private EFT sample can only include that single sample', colors.reset); exit(1)

    onlyCentralSample=False #Check whether all training samples are central samples
    centralVSpureEFT=False #Check whether training samples are part central / part pure-EFT
    onlySMEFT=False #Check whether all training samples are SM/EFT private samples
    ncentralSamples=0; nPureEFTSamples=0; nSMEFTSamples=0
    for label in labels_list:
        if("PrivMC" in label and "_c" in label): nPureEFTSamples+=1 #E.g. 'PrivMC_tZq_ctz'
        elif("PrivMC" in label): nSMEFTSamples+=1 #E.g. 'PrivMC_tZq'
        else: ncentralSamples+=1 #E.g. 'tZq'

    if nSMEFTSamples == len(labels_list): onlySMEFT=True
    elif (nPureEFTSamples+ncentralSamples) == len(labels_list): centralVSpureEFT=True
    elif nPureEFTSamples == 0 and nSMEFTSamples==0: onlyCentralSample=True
    else: print(colors.fg.red, 'ERROR : sample naming conventions not recognized, or incorrect combination of samples', colors.reset); exit(1)

    if opts["parameterizedDNN"]==True and onlySMEFT==False:
        print(colors.bold, colors.fg.red, 'Parameterized DNN supported for SM/EFT samples only ! Setting parameterizedDNN to False', colors.reset);
        opts["parameterizedDNN"]=False;

    if opts["regress"]==True:
        if opts["target"] != "class": print(colors.fg.red, 'ERROR : target name not available for regression yet', colors.reset); exit(1)

# //--------------------------------------------

#-- Initialization

    #Year selection
    lumiName = Get_LumiName(lumi_years)

    # Set main output paths
    weightDir = "../weights/DNN/" + lumiName + '/'
    os.makedirs(weightDir, exist_ok=True)

    #Top directory containing all input ntuples
    ntuplesDir = "../input_ntuples/"

    #Model output name
    h5modelName = weightDir + 'model.h5'

    #Determine/store number of process classes
    opts["nofOutputNodes"] = len(processClasses_list) #1 output node per class
    if opts["regress"] is True or opts["nofOutputNodes"] == 2: #Special case : 2 classes -> binary classification -> 1 output node only
        opts["nofOutputNodes"] = 1
    if opts["parameterizedDNN"]:
        opts["nofOutputNodes"] = len(opts["listOperatorsParam"])+1 #1 output node for SM and each EFT operator

    if opts["parameterizedDNN"] == True: opts["maxEvents"] = opts["nEventsPerPoint"]
    else: opts["maxEvents"] = opts["maxEventsPerClass"]

    if opts["parameterizedDNN"] == True: opts["batchSize"] = opts["batchSizeClass"]
    else: opts["batchSize"] = opts["batchSizeEFT"]

    return opts, lumiName, weightDir, ntuplesDir, h5modelName, opts["batchSize"]

# //--------------------------------------------
# //--------------------------------------------

# //--------------------------------------------
# //--------------------------------------------





# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
######## ######## ########
##       ##          ##
##       ##          ##
######   ######      ##
##       ##          ##
##       ##          ##
######## ##          ##
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

def CheckName_EFTpoint_ID(old):

    new = re.sub('p(\d+)', r'.\1', old)
    new = re.sub('min(\d+)', r'-\1', new)

    # print(old, ' --> ', new)
    return new

# //--------------------------------------------
# //--------------------------------------------
