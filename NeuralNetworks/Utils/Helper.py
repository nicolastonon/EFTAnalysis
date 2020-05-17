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
import shutil
from datetime import datetime
from scipy.stats import ks_2samp, anderson_ksamp, chisquare
from matplotlib import pyplot as plt
from Utils.ColoredPrintout import colors
from tensorflow.keras.models import load_model

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

    plt.close('all') #timer calls this function after N seconds and closes all active figures

    return

# //--------------------------------------------
# //--------------------------------------------
#-- Printout training info
def batchOutput(batch, logs):

    print("Finished batch: " + str(batch))
    print(logs)

    return

# //--------------------------------------------
# //--------------------------------------------

#-- Write NN input variables to a .txt file
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

#Using median and stddev from quantile is more robust against distributions with large tails
#q is the fraction of events that should be in the interval [-1, 1]
def get_normalization_iqr(np_array, q):

    df = pd.DataFrame(data=np_array[0:,0:]) #Convert to panda DF for easier manipulation

    q = (1 + q) / 2  #Tranform q from percentage to quantile

    newDF = pd.DataFrame()
    for key in df.keys():
        newDF[key] = df[key].apply(lambda x: np.mean(x[np.nonzero(x)]) if hasattr(x, "__len__") else x) #Not sure what is done here
    median = newDF.median() #Get median
    l = abs(newDF.quantile(1 - q) - median) #Get lower boundary corresponding to quantile
    r = abs(newDF.quantile(q) - median) #Get upper boundary corresponding to quantile
    maximums = abs(newDF.max()) #Also get the maximum value found (if no proper boundary found, can still divide by max value)
    # print('median\n', median)
    # print('maximums\n', maximums)
    for i, (il, ir) in enumerate(zip(l,r)):
        # print('il', il); print('ir', ir)
        if il == ir:
            print(colors.dim, f"[WARNING] feature {df.keys()[i]} has no width --> Set width = ", maximums[i], ' (max. value)', colors.reset) #Happens e.g. for discrete variables perfectly centered at 0. Better to return the max value, so that all values will effectively lie in [-1;+1]
            l[i] = maximums[i]; r[i] = maximums[i]

    return median.values, np.maximum(l, r).values #Return median and quantile boundary for rescaling

# //--------------------------------------------
# //--------------------------------------------

#Printout the output of the first (=input) layer here for N events, e.g. to verify that the normalization layer works properly
def Printout_Outputs_Layer(model, ilayer, xx):
    print('--------------------------------------------')
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

    #Need 1D arrays
    if values1.ndim > 1: values1 = np.squeeze(values1)
    if values2.ndim > 1: values2 = np.squeeze(values2)

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

def Load_PreExisting_Model(h5modelName):

# //--------------------------------------------
#-- Can access weights and biases of any layer
# weights_layer, biases_layer = model.layers[0].get_weights(); print(weights_layer.shape); print(biases_layer.shape); print(weights_layer); print(biases_layer[0:2])
#-- Loads the latest checkpoint weights
# latest = tensorflow.train.latest_checkpoint(ckpt_dir)
# tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
# model = load_model(_h5modelName) # model has to be re-loaded
# model.load_weights(latest)
# //-------------------------------------------

    tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    model = load_model(h5modelName) # model has to be re-loaded


    return model

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

#-- Get name corresponding to the data-taking years which are considered in the NN training
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

#Perform sanitify check of user options and set internal options accordingly (updates the option dictionnary given in arg)
def Initialization_And_SanityChecks(opts, lumi_years, processClasses_list, labels_list):

# //--------------------------------------------
#-- Initialization

    opts["parameterizedNN"] = False
    opts["regress"] = False
    if opts["strategy"] in ["CARL", "CARL_multiclass", "ROLR", "RASCAL"]: opts["parameterizedNN"] = True
    if opts["strategy"] in ["regressor", "ROLR", "RASCAL"]: opts["regress"] = True

    #Year selection
    lumiName = Get_LumiName(lumi_years)

    # Set main output paths
    weightDir = "../weights/NN/" + lumiName + '/'

    #Model output name
    h5modelName = weightDir + 'model.h5'

    if opts["makeValPlotsOnly"] == False: #If training a new NN, remove previous output folder
        # os.remove(weightDir+"/*")
        shutil.rmtree(weightDir)
        os.makedirs(weightDir, exist_ok=True)
    else:
        print(colors.fg.orange, "Will only produce validation plots. Reading pre-existing NN model:", colors.reset, h5modelName)

    print(colors.dim, "Created clean output directory:", colors.reset, weightDir)

    #Top directory containing all input ntuples
    ntuplesDir = "../input_ntuples/"

    #Determine/store number of process classes, depending on strategy
    opts["nofOutputNodes"] = len(processClasses_list) #Multiclass classification --> 1 output node per process class
    if opts["strategy"] is "classifier" and len(processClasses_list) == 2: opts["nofOutputNodes"] = 1 #Binary classification --> single output node needed
    elif opts["strategy"] is "regressor": opts["nofOutputNodes"] = 1
    elif opts["strategy"] in ["CARL", "CARL_singlePoint"]: opts["nofOutputNodes"] = 1 #Binary classification
    elif opts["strategy"] is "CARL_multiclass":
        if len(opts["listOperatorsParam"]) == 1: opts["nofOutputNodes"] = 1 #Binary classification
        else: opts["nofOutputNodes"] = len(opts["listOperatorsParam"])+1 #1 output node for SM and each EFT operator
    elif opts["strategy"] is "ROLR": opts["nofOutputNodes"] = 1 #Regress on r
    elif opts["strategy"] is "RASCAL": opts["nofOutputNodes"] = 1 + len(opts["listOperatorsParam"]) #Regress on r and t ; t has 1 component per EFT operator

    if opts["parameterizedNN"] == True:
        opts["maxEvents"] = opts["nEventsPerPoint"]
        opts["batchSize"] = opts["batchSizeEFT"]
    else:
        opts["maxEvents"] = opts["maxEventsPerClass"]
        opts["batchSize"] = opts["batchSizeClass"]

# //--------------------------------------------
#-- Sanity checks

    if opts["strategy"] not in ["classifier", "regressor", "CARL", "CARL_multiclass", "CARL_singlePoint", "ROLR", "RASCAL"]:
        print(colors.fg.red, 'ERROR : strategy', opts["strategy"],'not recognized', colors.reset); exit(1)

    if opts["strategy"] is "CARL_singlePoint" and opts["maxEvents"] != -1:
        print(colors.dim, "For strategy CARL_singlePoint, all available events will be used. Setting maxEvents to -1 !", colors.reset)
        opts["maxEvents"] = -1

    if opts["refPoint"] is "SM" and opts["strategy"] is "CARL_singlePoint": print(colors.fg.red, 'ERROR : refPoint must be different than SM for strategy CARL_singlePoint (choose point to train against SM) !', colors.reset); exit(1)
    if opts["refPoint"] is not "SM":
        if opts["strategy"] is "CARL_multiclass": print(colors.fg.red, 'ERROR : refPoint', opts["refPoint"],' must be set to SM for strategy CARL_multiclass (will train SM against several EFT operators) !', colors.reset); exit(1)
        if not opts["refPoint"].startswith("rwgt_c"): print(colors.fg.red, 'ERROR : refPoint', opts["refPoint"],' must start with prefix [rwgt_]. Example: [rwgt_ctZ_-2.6_ctp_3.1] !', colors.reset); exit(1)
        operatorNames, operatorWCs, idx_SM = Parse_EFTpoint_IDs(opts["refPoint"]) #NB: returns 2D arrays
        for iop in range(len(operatorNames[0])):
            if str(operatorNames[0][iop]) not in opts["refPoint"] or not np.any(operatorWCs[0]): print(colors.fg.red, 'ERROR : refPoint', opts["refPoint"],' not understood !', colors.reset); exit(1)

    if opts["score_lossWeight"] <= 0:
        print(colors.dim, "Option score_lossWeight should not be <=0. Setting it to 1 !", colors.reset)
        opts["score_lossWeight"] = 1

    if opts["nHiddenLayers"]<0 or opts["nNeuronsPerLayer"]<0 or opts["dropoutRate"]<0: print(colors.fg.red, 'ERROR : Invalid negative values found in NN architecture options !', colors.reset); exit(1)

    if opts["activHiddenLayers"] is "lrelu" or opts["activInputLayer"] is "lrelu":  print(colors.fg.red, 'ERROR : leaky relu not properly implemented yet (see Model.py) !', colors.reset); exit(1)

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
    ncentralSamples=0; nPureEFTSamples=0; nSMEFTSamples=0;
    for label in labels_list:
        if("PrivMC" in label and "_c" in label): nPureEFTSamples+=1 #E.g. 'PrivMC_tZq_ctz'
        elif("PrivMC" in label): nSMEFTSamples+=1 #E.g. 'PrivMC_tZq'
        else: ncentralSamples+=1 #E.g. 'tZq'

    totalSamples = len(processClasses_list)
    if nSMEFTSamples == len(labels_list): onlySMEFT=True
    elif (nPureEFTSamples+ncentralSamples) == len(labels_list): centralVSpureEFT=True
    elif nPureEFTSamples == 0 and nSMEFTSamples==0: onlyCentralSample=True
    else: print(colors.fg.red, 'ERROR : sample naming conventions not recognized, or incorrect combination of samples', colors.reset); exit(1)

    if (opts["parameterizedNN"]==True or opts["strategy"] not in ["classifier", "regressor"]) and onlySMEFT==False: print(colors.bold, colors.fg.red, 'This NN strategy is supported for SM+EFT samples only !', colors.reset); exit(1)
    elif opts["strategy"] in ["classifier", "regressor"] and nSMEFTSamples > 0: print(colors.bold, colors.fg.red, 'This NN strategy is not supported for SM+EFT samples !', colors.reset); exit(1)
    if totalSamples < 2 and opts["strategy"] is "classifier": print(colors.bold, colors.fg.red, 'Classifier strategy requires at least 2 samples !', colors.reset); exit(1)
    if opts["nPointsPerOperator"] < 2: print(colors.bold, colors.fg.red, 'Parameter nPointsPerOperator must be >= 2 !', colors.reset); exit(1)

# //--------------------------------------------

    Write_Timestamp_toLogfile(weightDir, 0)
    Dump_NN_Options_toLogFile(opts, weightDir)

    return lumiName, weightDir, h5modelName, ntuplesDir, opts["batchSize"]

# //--------------------------------------------
# //--------------------------------------------

#Write information related to this NN training
def Dump_NN_Options_toLogFile(opts, weightDir):

    #-- Also append the names of the input/output nodes in the file "NN_info.txt" containing input features names, etc. (for later use in C++ code)
    text_file = open(weightDir + "NN_settings.txt", "a+") #Overwrite file

    text_file.write("\nOPTIONS\n")
    text_file.write("----------------- \n")

    for key in opts:
        text_file.write(str(key) + " --> " + str(opts[key]) + "\n")

    text_file.write("----------------- \n\n")
    text_file.close()

    return

# //--------------------------------------------
# //--------------------------------------------

#Write timestamp to NN logfile
#status: 0=start, 1=end
def Write_Timestamp_toLogfile(weightDir, status):

    if status == 0: mode = "w"
    else: mode = "a+"

    dateTimeObj = datetime.now()
    timestampStr = dateTimeObj.strftime("%d-%b-%Y (%H:%M:%S)")
    # print('Current Timestamp : ', timestampStr)

    text_file = open(weightDir + "NN_settings.txt", mode) #Overwrite file
    if status == 0: text_file.write("Start of NN training :" + str(timestampStr) + "\n")
    elif status == 1: text_file.write("End of NN training and evaluation :" + str(timestampStr) + "\n")
    text_file.close()

    return
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

def Remove_Unnecessary_EFTweights(array_EFTweights, array_EFTweightIDs):
    """
    Look for EFT weight names which do not follow the expected naming convention (of the kind 'rwgt_ctZ_5.2'), and removes coherently these elements in the arrays of EFT weights/names. Returns modified arrays.

    Parameters:
    array_EFTweights (ndarray of shape [n_events, n_points]) : reweights for all points, for all events
    array_EFTweightIDs (ndarray of shape [n_events, n_points]) : reweight names for all points, for all events. Example name: 'rwgt_ctZ_5_ctW_3'

    Returns:
    Same arrays, but without the columns (<-> reweights) which are not necessary for the morphing procedure (e.g. SM reweight which does not follow the same naming convention)
    """

    #Remove coherently the incorrect weights
    indices = np.char.find(array_EFTweightIDs[0,:].astype('U'), 'rwgt_c') #Get column indices of reweight names not following the expected naming convention, from first event; 0 <-> found ; -1 <-> not found
    indices = np.asarray(indices, dtype=bool) #Convert indices to booleans (True <-> does not follow naming convention)
    array_EFTweightIDs = np.delete(array_EFTweightIDs, indices==True, axis=1)
    array_EFTweights = np.delete(array_EFTweights, indices==True, axis=1)

    #Just removes first weight if it is 'rwgt_1' (baseline weight sometimes included by default by MG)
    # if array_EFTweightIDs[0][0] is "rwgt_1":
    #     array_EFTweightIDs = np.delete(array_EFTweightIDs, 0, axis=1)
    #     array_EFTweights = np.delete(array_EFTweights, 0, axis=1)

    # print(array_EFTweights.shape); print(array_EFTweightIDs.shape)
    return array_EFTweights, array_EFTweightIDs

# //--------------------------------------------
# //--------------------------------------------

  #####    ##   #####   ####  ######
  #    #  #  #  #    # #      #
  #    # #    # #    #  ####  #####
  #####  ###### #####       # #
  #      #    # #   #  #    # #
  #      #    # #    #  ####  ######

def Parse_EFTpoint_IDs(benchmarkIDs):
    """
    Parse an array of strings, each corresponding to the name of an EFT point

    Parameters:
    benchmarkIDs (ndarray of shape [n_points]) : array of strings, each corresponding to a separate EFT point. Example : 'rwgt_ctZ_5_ctW_3'

    Returns:
    operatorNames (ndarray of shape [n_points, n_operators]) : array of names of EFT operators, for all (n_points) EFT points
    operatorWCs (ndarray of shape [n_points, n_operators]) : array whose columns represent the WC value of each operator defining a given point, and rows represent different EFT points
    idx_SM : return index corresponding to SM point, if corresponding string is found
    """

    benchmarkIDs = np.atleast_1d(benchmarkIDs) #If a single point is given in arg, make sure it is treated as an array in the function (and not as a string)

    idx_SM = -1 #Store SM index
    operatorNames = []
    operatorWCs = []
    # for idx in range(len(benchmarkIDs)):
    for idx, ID in enumerate(benchmarkIDs):

        #For each EFT point, get the list of names/WCs for all operators
        list_operatorNames = []
        list_operatorWCs = []
        # ID = benchmarkIDs[idx]

        if not ID.startswith("rwgt_c"): #Every considered EFT operator is expected to start with this substring ; for others (e.g. 'rwgt_sm'), don't parse
            # list_operatorNames.append(ID) #Operator name
            # list_operatorWCs.append(float(0)) #Operator WC
            if ID=="rwgt_sm" or ID=="rwgt_SM": idx_SM = idx #SM point found
            continue

        ID = CheckName_EFTpoint_ID(ID) #Enforce proper naming convention

        prefix = 'rwgt_' #Naming convention, strip this substring
        if ID.startswith(prefix): ID = ID[len(prefix):]
        else: print('Error : naming convention in benchmark ID not recognized'); exit(1)

        list_keys = ID.split('_') #Split string into list of substrings (pairs [name,WC])
        for ikey in range(0, len(list_keys)-1, 2):
            list_operatorNames.append(list_keys[ikey]) #Operator name
            list_operatorWCs.append(float(list_keys[ikey+1])) #Operator WC

        #Append list for each EFT point
        operatorNames.append(list_operatorNames)
        operatorWCs.append(list_operatorWCs)

    #Convert list of lists into 2D array. Each list element must have an equal length
    operatorNames = np.array(operatorNames)
    operatorWCs = np.array(operatorWCs)

    return operatorNames, operatorWCs, idx_SM

# //--------------------------------------------
# //--------------------------------------------

#Translate reference points IDs (if different from SM) into the corresponding 2D array of WC values, properly ordered following same operator order as found in the sample
def Translate_EFTpointID_to_WCvalues(operatorNames_sample, refPointIDs):

    refPointIDs = np.atleast_1d(refPointIDs) #If a single point is given in arg, make sure it is treated as an array in the function (and not as a string)

    if refPointIDs[0] in ["SM", "sm"]:
        return np.zeros(len(operatorNames_sample)) #SM point <-> all WCs are null

    operatorNames_new, operatorWCs_new, _ = Parse_EFTpoint_IDs(refPointIDs)

    orderedList_WCvalues_allPoints = []
    for i_ID in range(len(refPointIDs)): #For each point ID
        orderedList_WCvalues = []
        for op_sample in operatorNames_sample: #For each operator found in sample
            for iop_ref in range(operatorNames_new.shape[1]): #For each operator in ID
                if str(operatorNames_new[i_ID][iop_ref]) == str(op_sample): orderedList_WCvalues.append(operatorWCs_new[i_ID][iop_ref]) #Append WC to list at correct position (according to ordering in sample)
                elif operatorNames_new[i_ID][iop_ref] not in operatorNames_sample: print(colors.bold, colors.fg.red, 'ERROR : refPoint seems not to be compatible with operators included in this sample... (check exact naming)', colors.reset); exit(1) #Operator specified in ID not found in sample

        orderedList_WCvalues_allPoints.append(orderedList_WCvalues)

    orderedList_WCvalues_allPoints = np.asarray(orderedList_WCvalues_allPoints) #List --> array

    return orderedList_WCvalues_allPoints

# //--------------------------------------------
# //--------------------------------------------
