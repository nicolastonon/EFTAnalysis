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
import math
import os
from os import path
import shutil
from datetime import datetime
from scipy.stats import ks_2samp, anderson_ksamp, chisquare
from matplotlib import pyplot as plt
from Utils.ColoredPrintout import colors
from tensorflow.keras.models import load_model
from Utils.LossOptimMetric import Get_Loss_Optim_Metrics
from Utils.InputFeatures import *


#-- Top directory containing all input ntuples
ntuplesDir = "../input_ntuples/" #LOCAL
# ntuplesDir = "/nfs/dust/cms/user/ntonon/CMSSW_10_2_20/src/potato_nicolas/potato-nicolas/nicolas/output/Analyzer3l-V10-AllSamples-d20201023-t151834/merged_ntuples/" #CMSSW


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

#-- Using median and stddev from quantile is more robust against distributions with large tails
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
            print(colors.dim, "[WARNING] feature", {df.keys()[i]}, "has no width --> Set width = ", maximums[i], ' (max. value)', colors.reset) #Happens e.g. for discrete variables perfectly centered at 0. Better to return the max value, so that all values will effectively lie in [-1;+1]
            l[i] = maximums[i]; r[i] = maximums[i]

    return median.values, np.maximum(l, r).values #Return median and quantile boundary for rescaling

# //--------------------------------------------
# //--------------------------------------------

#Printout the output of the first (=input) layer here for N events, e.g. to verify that the normalization layer works properly
def Printout_Outputs_Layer(model, ilayer, xx):
    print('\n', colors.bg.orange,'--------------------------------------------', colors.reset)
    print('LAYER ', ilayer,)
    get_layer_output = keras.backend.function([model.layers[0].input], [model.layers[ilayer].output])
    layer_output = get_layer_output([xx])[0]
    print("\n", layer_output)

# //--------------------------------------------
# //--------------------------------------------

#Printout the output of the first (=input) layer here for N events, e.g. to verify that the normalization layer works properly
def Printout_Weights_Layer(model):
    print('\n', colors.bg.orange,'--------------------------------------------', colors.reset)
    for ilayer, layer in enumerate(model.layers):
        print('\n\nLAYER ', ilayer, ' (', layer.name, ')')
        weights = layer.get_weights() # list of numpy arrays #First array=kernel (weights), second array=biases, then NB/PRelu/... ?
        for iarr in range(len(weights)):
            str = '\n-- Weights:'
            if 'batch_normalization' in layer.name:
                if iarr==0: str = '\n-- Gamma:'
                elif iarr==1: str = '\n-- Beta:'
                elif iarr==2: str = '\n-- Mean:'
                elif iarr==3: str = '\n-- Variance:'
            elif iarr==1: str = '\n-- Biases:'
            elif iarr>1: str = '\n-- (?):'
            print(str, weights[iarr])

# //--------------------------------------------
# //--------------------------------------------

#Shuffles coherently the rows of N arrays of same length
def unison_shuffled_copies(*arr):
    # print(len(arr))
    assert len(arr) > 1
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

''' #Obsolete?
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
'''

# //--------------------------------------------
# //--------------------------------------------

def Load_PreExisting_Model(h5modelName):
    '''
    #-- Can access weights and biases of any layer
    # weights_layer, biases_layer = model.layers[0].get_weights(); print(weights_layer.shape); print(biases_layer.shape); print(weights_layer); print(biases_layer[0:2])
    #-- Loads the latest checkpoint weights
    # latest = tensorflow.train.latest_checkpoint(ckpt_dir)
    # tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    # model = load_model(_h5modelName) # model has to be re-loaded
    # model.load_weights(latest)
    '''

    tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    model = load_model(h5modelName, compile=False) # model has to be re-loaded #compile=False <-> does not need to define any custom loss, since not needed for testing

    return model

# //--------------------------------------------
# //--------------------------------------------

def truncate(number, digits): #-> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper

def trunc_np(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)

# //--------------------------------------------
# //--------------------------------------------

def lineno():
    """Returns the current line number in our program."""
    import inspect
    return inspect.currentframe().f_back.f_lineno

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
def Initialization_And_SanityChecks(opts, lumi_years, processClasses_list, labels_list, list_features):

# //--------------------------------------------
#-- Initialization

    opts["trainAtManyEFTpoints"] = False
    opts["regress"] = False
    if opts["strategy"] in ["CARL", "CARL_multiclass", "ROLR", "RASCAL", "CASCAL"]: opts["trainAtManyEFTpoints"] = True #By construction these methods require to train the NN on a range of different hypotheses
    if opts["strategy"] in ["regressor", "ROLR", "RASCAL"]: opts["regress"] = True #Regressor strategies

    #-- NB: introduced difference between 'mixed-EFT' strategy (<-> requires training over many SMEFT hypotheses, e.g. CARL) and 'parameterized' strategies (for which we choose to add EFT WCs as additional inputs) #NB: parameterized NNs must be 'mixed-EFT', but mixed-EFT NNs don't need to be parameterized
    if opts["trainAtManyEFTpoints"] is False:
        if opts["parameterizedNN"] is True:
            opts["parameterizedNN"] = False
            print(colors.fg.orange, colors.bold, "\nWarning: you are not training on a mixture of many EFT points, hence option 'parameterizedNN' has been forced to False ! \n", colors.reset)

    #Year selection
    lumiName = Get_LumiName(lumi_years)

# //--------------------------------------------
#-- Define few necessary options

    if opts["testToy1D"] == True:
        if opts["strategy"] != "CARL": print(colors.fg.red, colors.bold, "\n\nERROR ! Debug option [testToy1D] is only valid for strategy CARL ! \n", colors.reset, h5modelName, '\n'); exit(1)
        print(colors.fg.orange, colors.bold, "\n\nWill train a NN on a dummy 1D toy example with only 1 input feature (dummy gaussian centered on theta parameter). See ref. article on paramNN.  HARDCODED FOR CTW OPERATOR ONLY !\n", colors.reset, h5modelName, '\n')
        opts["listOperatorsParam"] = ['ctw']
        opts["nPointsPerOperator"] = 3
        opts["minWC"] = -3
        opts["maxWC"] = 3

    #Determine/store number of process classes, depending on strategy
    opts["nofOutputNodes"] = len(processClasses_list) #Multiclass classification --> 1 output node per process class
    if opts["strategy"] is "classifier" and len(processClasses_list) == 2: opts["nofOutputNodes"] = 1 #Binary classification --> single output node needed
    elif opts["strategy"] is "regressor":
        opts["nofOutputNodes"] = 1
        if len(opts["targetVarIdx"])>1: opts["nofOutputNodes"] = len(opts["targetVarIdx"])
    elif opts["strategy"] in ["CARL", "CARL_singlePoint"]: opts["nofOutputNodes"] = 1 #Binary classification
    elif opts["strategy"] is "CARL_multiclass":
        if len(opts["listOperatorsParam"]) == 1: opts["nofOutputNodes"] = 1 #Binary classification
        else: opts["nofOutputNodes"] = len(opts["listOperatorsParam"])+1 #1 output node for SM and each EFT operator
        if opts["parameterizedNN"] is False: print(colors.fg.red, 'ERROR : strategy', opts["CARL_multiclass"],' requires option [parameterizedNN==True] (else need to rethink how to avoid cutting on input WCs to determine correct class, ...)', colors.reset); exit(1)
    elif opts["strategy"] is "ROLR": opts["nofOutputNodes"] = 1 #Regress on r
    elif opts["strategy"] in ["RASCAL","CASCAL"]: opts["nofOutputNodes"] = 1 + len(opts["listOperatorsParam"]) #Regress on r and t ; t has 1 component per EFT operator

    if opts["trainAtManyEFTpoints"] == True: #Use different data sampling options depending on case
        opts["maxEvents"] = opts["nEventsPerPoint"]
        opts["batchSize"] = opts["batchSizeEFT"]
    else:
        opts["maxEvents"] = opts["maxEventsPerClass"]
        opts["batchSize"] = opts["batchSizeClass"]

# //--------------------------------------------
#-- Sanity checks

    if opts["strategy"] not in ["classifier", "regressor", "CARL", "CARL_multiclass", "CARL_singlePoint", "ROLR", "RASCAL", "CASCAL"]:
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

    if "nNeuronsAllHiddenLayers" in opts and opts["nNeuronsAllHiddenLayers"]<0: print(colors.fg.red, 'ERROR : Invalid value for option [nNeuronsAllHiddenLayers] !', colors.reset); exit(1)
    if "nNeuronsAllHiddenLayers" not in opts and "nNeuronsPerHiddenLayer" not in opts: print(colors.fg.red, 'ERROR : number of neurons not set !', colors.reset); exit(1)
    if opts["nHiddenLayers"]<0 or opts["dropoutRate"]<0: print(colors.fg.red, 'ERROR : Invalid negative values found in NN architecture options !', colors.reset); exit(1)

    if "nNeuronsPerHiddenLayer" in opts and len(opts["nNeuronsPerHiddenLayer"]) != opts["nHiddenLayers"]: print(colors.fg.red, 'ERROR : list [nNeuronsPerHiddenLayer] must have a length == [nHiddenLayers] !', colors.reset); exit(1)

    if opts["activHiddenLayers"] is '': print(colors.fg.red, 'ERROR : Empty activation function for hidden layers !', colors.reset); exit(1)

    if opts["optimizer"] not in ['Adam','Nadam','Adadelta','AdaBound','RMSprop','SGD']: print(colors.fg.red, "ERROR: unknown optimizer algorithm", opts["optimizer"], colors.reset); exit(1)
    if opts["learnRate"] <= 0: print(colors.fg.red, "Wrong learning rate value ", opts["learnRate"], colors.reset); exit(1)

    if opts["regularizer"][0] is 'L1': opts["regularizer"][0] = 'l1'
    elif opts["regularizer"][0] is 'L2': opts["regularizer"][0] = 'l2'
    elif opts["regularizer"][0] is 'L1L2': opts["regularizer"][0] = 'l1l2'
    if opts["regularizer"][0] not in ['','l1','l2','l1l2'] or opts["regularizer"][1]<0: print(colors.fg.red, 'ERROR : wrong value for arg. regularizer !', colors.reset); exit(1)

    if len(processClasses_list) == 0:
        print(colors.fg.red, 'ERROR : no process class defined...', colors.reset); exit(1)
    elif len(processClasses_list) is not len(labels_list):
        print(colors.fg.red, 'ERROR : sizes of lists processClasses_list and labels_list are different...', colors.reset); exit(1)

    for procClass, label in zip(processClasses_list, labels_list):
        if "PrivMC" in label and "_c" not in label and len(procClass) > 1: print(colors.fg.red, 'ERROR: process classes containing a private EFT sample can only include that single sample. Maybe the process class labels must be changed (following conventions, cf. Helper.py). \nNB: if you are trying to set up a parameterized NN which considers the combination of >1 physics processes, just add them independently to the lists (merged automatically) !', colors.reset); exit(1)

    onlyCentralSample=False #Check whether all training samples are central samples
    centralVSpureEFT=False #Check whether training samples are part central / part pure-EFT
    onlySMEFT=False #Check whether all training samples are SM/EFT private samples
    ncentralSamples=0; nPureEFTSamples=0; nSMEFTSamples=0;
    for label in labels_list:
        if("PrivMC" in label and "_c" in label): nPureEFTSamples+=1 #E.g. 'PrivMC_tZq_ctz'
        elif("PrivMC" in label): nSMEFTSamples+=1 #E.g. 'PrivMC_tZq'
        else: ncentralSamples+=1 #E.g. 'tZq'

    totalSamples = len(processClasses_list)
    if nPureEFTSamples == 0 and nSMEFTSamples==0: onlyCentralSample=True
    elif nSMEFTSamples == len(labels_list): onlySMEFT=True
    elif (nPureEFTSamples+ncentralSamples) == len(labels_list): centralVSpureEFT=True
    else: print(colors.fg.red, 'ERROR : sample naming conventions not recognized, or incorrect combination of samples', colors.reset); exit(1)

    opts["samplesType"] = ""
    if onlyCentralSample: opts["samplesType"] = "onlyCentralSample"
    elif centralVSpureEFT: opts["samplesType"] = "centralVSpureEFT"
    elif onlySMEFT: opts["samplesType"] = "onlySMEFT"

    if (opts["trainAtManyEFTpoints"]==True or opts["strategy"] not in ["classifier", "regressor"]) and onlySMEFT==False: print(colors.bold, colors.fg.red, 'This NN strategy is supported for SM+EFT samples only !', colors.reset); exit(1) #Can only train on different EFT hypotheses if using SMEFT samples
    # elif opts["strategy"] in ["classifier", "regressor"] and nSMEFTSamples > 0: print(colors.bold, colors.fg.red, 'This NN strategy is not supported for SM+EFT samples !', colors.reset); exit(1)
    if totalSamples < 2 and opts["strategy"] is "classifier": print(colors.bold, colors.fg.red, 'Classifier strategy requires at least 2 samples !', colors.reset); exit(1)
    if opts["nPointsPerOperator"] < 2: print(colors.bold, colors.fg.red, 'Parameter nPointsPerOperator must be >= 2 !', colors.reset); exit(1)

    if "listMinMaxWC" in opts and len(opts["listMinMaxWC"]) != 2*len(opts["listOperatorsParam"]): print(colors.bold, colors.fg.red, 'ERROR : Length of list [listMinMaxWC] (', len(opts["listMinMaxWC"]), ') must be exactly twice that of [listOperatorsParam] (', len(opts["listOperatorsParam"]), ') !', colors.reset); exit(1)

    opts["loss"], _, opts["metrics"], _ = Get_Loss_Optim_Metrics(opts) #NB: these options are not used anywhere (will be obtained again in main function) ! Only read here so that they can be dumped into the logfile

    if opts["strategy"] is "regressor" and (not all([idx>=0 for idx in opts["targetVarIdx"]])): print(colors.bold, colors.fg.red, 'Wrong option [targetVarIdx] or [comparVarIdx] !', colors.reset); exit(1)
    elif opts["strategy"] is not "regressor": opts["targetVarIdx"] = []; opts["comparVarIdx"] = -1 #Set default
    elif opts["nofOutputNodes"] > 1 and opts["comparVarIdx"]>=0:
        print(colors.fg.red, '\nNB: there is more than 1 output node --> will *not* compare predictions to selected variable... !', colors.reset)
        opts["comparVarIdx"] = -1
    else:
        print('\n', colors.fg.orange, colors.underline, '-- Will use the following variable(s) as target(s) :', colors.reset, [list_features[opts["targetVarIdx"][idx]] for idx in opts["targetVarIdx"]])
        if(opts["comparVarIdx"] >= 0): print('\n\n', colors.fg.orange, colors.underline, '-- Will compare predictions to the following variable:', colors.reset, list_features[opts["comparVarIdx"]], '\n\n')

    #-- NN strategy (for easier interfacing with my analysis code)
    opts["NN_strategy"] = "MVA_SM" #Default, can use full MVA distribution directly in Combine
    if centralVSpureEFT is True or opts["strategy"] is "CARL_singlePoint" or (opts["trainAtManyEFTpoints"] is True and opts["parameterizedNN"] is False): opts["NN_strategy"] = "MVA_EFT" #Will need to consider each MVA bin separately, for individual EFT parametrization
    elif opts["parameterizedNN"] is True: opts["NN_strategy"] = "MVA_param" #Will need to produce MVA templates for each and every point considered for signal extraction

    #-- Protections on dataset splitting
    if opts["trainAtManyEFTpoints"] == True and (opts["nEventsTot_train"]!=-1. or opts["nEventsTot_val"]!=-1. or opts["nEventsTot_test"]!=-1): #Mixed-EFT strategy <-> don't consider these options
        print(colors.fg.red, 'Warning: dataset splitting options [nEventsTot_train/nEventsTot_val/nEventsTot_test] can not be used when training at many EFT points... Will use [splitTrainValTestData] proportions instead !', colors.reset)
        opts["nEventsTot_train"]=-1.; opts["nEventsTot_val"]=-1.; opts["nEventsTot_test"]=-1.;

    if opts["splitTrainValTestData"][0]>0. and (opts["nEventsTot_train"]!=-1. or opts["nEventsTot_val"]!=-1. or opts["nEventsTot_test"]!=-1): #splitTrainValTestData option superseeds these other options
        print(colors.fg.red, 'Warning: Overriding dataset splitting options [nEventsTot_train/nEventsTot_val/nEventsTot_test] with superseeding options [splitTrainValTestData] !', colors.reset)
        opts["nEventsTot_train"]=-1.; opts["nEventsTot_val"]=-1.; opts["nEventsTot_test"]=-1.;

    if len(opts["splitTrainValTestData"]) is not 3 or opts["splitTrainValTestData"][0]==1. or (opts["splitTrainValTestData"][0]>0. and opts["splitTrainValTestData"][0]+opts["splitTrainValTestData"][1]+opts["splitTrainValTestData"][2] != 1.): #Check option validity
        print(colors.fg.red, 'ERROR: Wrong option [splitTrainValTestData]', colors.reset); exit(1)
    if opts["nEventsTot_train"]!=-1. and opts["nEventsTot_val"]==-1. and opts["nEventsTot_test"]==-1.: #Check option validity
        print(colors.fg.red, 'ERROR: Wrong values for options [nEventsTot_train/nEventsTot_val/nEventsTot_test]', colors.reset); exit(1)

    if (opts["splitTrainValTestData"][1] != 0. and opts["splitTrainValTestData"][2] == 0.) or (opts["splitTrainValTestData"][2] != 0. and opts["splitTrainValTestData"][1] == 0.) or (opts["nEventsTot_val"] == -1. and opts["nEventsTot_test"] != -1.) or (opts["nEventsTot_val"] != -1. and opts["nEventsTot_test"] == -1.): #Convention: want nTest to be non-zero, not nVal (if not requiring independent val/test sets)
        print(colors.fg.orange, 'NB: you have either set testData=0 or valData=0; will hence use the same data for both val/test !', colors.reset)
        if opts["nEventsTot_train"] != -1. and opts["nEventsTot_test"] == -1.:
            opts["nEventsTot_test"] = opts["nEventsTot_val"]
            opts["nEventsTot_val"] = -1.
        elif opts["splitTrainValTestData"][0] != 0 and opts["splitTrainValTestData"][2] == 0.: #Convention: want nTest to be non-zero, not nVal
            opts["splitTrainValTestData"][2] = opts["splitTrainValTestData"][1]
            opts["splitTrainValTestData"][1] = 0.


# //--------------------------------------------
#-- Set/update list of input features

    #-- User can choose to use a case-specific list of input features (hard-coded in 'InputFeatures.py'); otherwise, use list of features defined in the main code
    if opts["useHardCodedListInputFeatures"] == True:
        list_features = [] #Re-initialize

        #-- Hardcode different use-cases here
        if opts["strategy"] is "classifier": #NN-SM
            list_features = features_SM

        elif (opts["trainAtManyEFTpoints"] == True and len(opts["listOperatorsParam"])==1) or opts["strategy"]=="CARL_singlePoint":
            print(len(opts["listOperatorsParam"]))
            if 'tZq' in labels_list[0]:
                if (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=='ctz') or 'ctz' in opts["refPoint"]: list_features = features_CARL_tZq_ctz
                elif (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="ctw") or 'ctw' in opts["refPoint"]: list_features = features_CARL_tZq_ctw
                elif (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="cpqm") or 'cpqm' in opts["refPoint"]: list_features = features_CARL_tZq_cpqm
                elif (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="cpq3") or 'cpq3' in opts["refPoint"]: list_features = features_CARL_tZq_cpq3
                elif (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="cpt") or 'cpt' in opts["refPoint"]: list_features = features_CARL_tZq_cpt
            elif 'ttZ' in labels_list[0]:
                if (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="ctz") or 'ctz' in opts["refPoint"]: list_features = features_CARL_ttZ_ctz
                elif (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="ctw") or 'ctw' in opts["refPoint"]: list_features = features_CARL_ttZ_ctw
                elif (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="cpqm") or 'cpqm' in opts["refPoint"]: list_features = features_CARL_ttZ_cpqm
                elif (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="cpq3") or 'cpq3' in opts["refPoint"]: list_features = features_CARL_ttZ_cpq3
                elif (opts["refPoint"]=="SM" and opts["listOperatorsParam"][0]=="cpt") or 'cpt' in opts["refPoint"]: list_features = features_CARL_ttZ_cpt

        elif (opts["trainAtManyEFTpoints"] == True and len(opts["listOperatorsParam"])==5):
            if 'tZq' in labels_list[0]: list_features = features_CARL_tZq_all
            elif 'ttZ' in labels_list[0]: list_features = features_CARL_ttZ_all

        if len(list_features)==0: print(colors.fg.red, 'ERROR : option [useHardCodedListInputFeatures=False], but can not find a dedicated list of input features for the particular use-case you are currently considering (cf. user-options). Set to True to use the list defined in the main code, or define the use-case in [InputFeatures.py and Helper.py] !', colors.reset); exit(1)

    # '''
    if opts["useLowLevelFeatures"]: #Also include (hardcoded) low-level input features
        list_features.append("lep1_pt")
        list_features.append("lep2_pt")
        list_features.append("lep3_pt");
        list_features.append("lep1_eta")
        list_features.append("lep2_eta")
        list_features.append("lep3_eta")
        list_features.append("lep1_phi")
        list_features.append("lep2_phi")
        list_features.append("lep3_phi")

        list_features.append("jet1_pt")
        list_features.append("jet2_pt")
        list_features.append("jet3_pt")
        list_features.append("jet1_eta")
        list_features.append("jet2_eta")
        list_features.append("jet3_eta")
        list_features.append("jet1_phi")
        list_features.append("jet2_phi")
        list_features.append("jet3_phi")
        list_features.append("jet1_DeepJet")
        list_features.append("jet2_DeepJet")
        list_features.append("jet3_DeepJet")

        # list_features.append("jet4_pt")
        # list_features.append("jet4_eta")
        # list_features.append("jet4_phi")
        # list_features.append("jet4_DeepJet")

        ## list_features.append("jet1_DeepCSV")
        ## list_features.append("jet2_DeepCSV")
        ## list_features.append("jet3_DeepCSV")
        ## list_features.append("jet4_DeepCSV")
    # '''

# //--------------------------------------------
#-- Define path-related variables

    weightDir = "../weightsMVA/NN/" + lumiName + '/' #Base output dir

    if opts ["storeInTestDirectory"] is True: #Store all output files in common test dir., and overwrite previous files
        weightDir+= 'tmp/'
        print(colors.dim, "Created clean test output directory:", colors.reset, weightDir)
    else: #Store all output files in training-specific subdir., following same path conventions as in main analysis code
        if opts["NN_strategy"] is "MVA_SM":
            weightDir+= 'SM/'
            if opts["nofOutputNodes"] == 1:
                if "tZq" in labels_list[0]: weightDir+= 'tZq/'
                elif "ttZ" in labels_list[0]: weightDir+= 'ttZ/'
                else: weightDir+= 'Other/'
            else: weightDir+= 'Multiclass/'
        else:
            weightDir+= 'EFT/'
            if "tZq" in labels_list[0]: weightDir+= 'tZq/'
            elif "ttZ" in labels_list[0]: weightDir+= 'ttZ/'
            else: weightDir+= 'Other/'

    if opts["storeInTestDirectory"] == False and opts["storePerOperatorSeparately"] == True and opts["strategy"] in ["CARL","CARL_singlePoint","CARL_multiclass","ROLR","RASCAL","CASCAL"]: #Store in dedicated operator-dependent output dir.
        if len(opts["listOperatorsParam"])==1: weightDir+= opts["listOperatorsParam"][0] + '/'
        elif len(opts["listOperatorsParam"])==5: weightDir+= 'all' + '/'

    #Model output name
    h5modelName = weightDir + 'model.h5'

    if opts["makeValPlotsOnly"] == False: #If training a new NN, remove previous output folder
        # os.remove(weightDir+"/*")
        if path.exists(weightDir): shutil.rmtree(weightDir, ignore_errors=True)
        os.makedirs(weightDir, exist_ok=True)
    else: print(colors.fg.orange, colors.bold, "\n\nWill only produce validation plots. Reading pre-existing NN model:\n", colors.reset, h5modelName, '\n')

# //--------------------------------------------

    #-- Overwrite previously existing file if any (then, other codes will open it in 'append' mode)
    if opts["makeValPlotsOnly"]==False:
        # shutil.move(weightDir + "NN_info.txt", weightDir + "NN_info_bck.txt") #First, copy previous NN_info file to bck path... can be useful if rerun code by mistake and previous file gets overwritten
        text_file = open(weightDir + "NN_info.txt", "w") #'w' to overwrite
        text_file.close()

    Write_Timestamp_toLogfile(weightDir, 0)
    Dump_NN_Options_toLogFile(opts, weightDir) #Write user-options to dedicated logfile

    return lumiName, weightDir, h5modelName, opts["batchSize"], list_features

# //--------------------------------------------
# //--------------------------------------------

#Write information related to this NN training
#NB: also append the names of the input/output nodes in separate output file "NN_info.txt" containing names of input features, etc. (for later use in C++ code)
def Dump_NN_Options_toLogFile(opts, weightDir):

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

def Get_ListPointsSampling_SingleOp(operator_scan, range_step):
    '''
    Set list of 'WC points' for which output plots will be produced.

    Return list of point names, and the corresponding WC values.
    '''

    list_points_sampling = ['SM']
    WCs = [0]

    nPoints = 1 + (range_step[1] - range_step[0]) / range_step[2]
    # print('nPoints', nPoints)
    for x in np.linspace(range_step[0], range_step[1], num=int(nPoints)):
        # print('x', x)
        # if x==0: continue
        point = 'rwgt_' + operator_scan + '_' + str(x)
        list_points_sampling.append(point)
        WCs.append(x)

    # print('list_points_sampling: ', list_points_sampling)
    return list_points_sampling, WCs

# //--------------------------------------------
# //--------------------------------------------

def Make_Animation_fromParamOutputPlots(standaloneValDir, list_labels, list_points_sampling, operator_scan, WCs, ROC=False):

    delay = 100 #in 1/100th of a second
    outname = standaloneValDir+'Overtrain_paramNN.gif' #Output name
    if ROC: outname = standaloneValDir+'ROCs_paramNN.gif'

    cmd = 'convert -delay ' + str(delay) + ' -loop 0 '

    for ipt in range(len(list_points_sampling)):
        if ipt==0: continue #Don't want SM point
        WC = str(WCs[ipt]).replace('.0','')
        # WC = str(WCs[ipt])
        plotname = standaloneValDir + 'Overtraining_NN_' + list_labels[0] + '_' + operator_scan + '_' + WC + '.png'
        if ROC: plotname = standaloneValDir + 'ROC' + '_' + operator_scan + '_' + WC + '.png'
        cmd+= plotname + ' '

    cmd+= outname
    # print(cmd)
    os.system(cmd)

    #Additional command to add the same gif in reverse, and overwrite => continuous loop without break
    cmd = "convert " + outname +  " -coalesce -duplicate 1,-2-1 -quiet -layers OptimizePlus -loop 0 " + outname
    # print(cmd)
    os.system(cmd)

    print(colors.fg.lightgrey,'\n---> Created animated GIF: ', colors.reset, outname)

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
    # new = new.lower()

    # print(old, ' --> ', new)
    return new

# //--------------------------------------------
# //--------------------------------------------

def Remove_Unnecessary_EFTweights(array_EFTweights, array_EFTweightIDs):
    """
    Look for EFT weight names which do not follow the expected naming convention (of the kind 'rwgt_ctZ_5.2'), and removes coherently these elements in the arrays of EFT weights/names. Returns modified arrays.

    NB: also remove SM benchmark weight, voluntarily (already stored in dedicated array anyway)

    Parameters:
    array_EFTweights (ndarray of shape [n_events, n_points]) : reweights for all points, for all events
    array_EFTweightIDs (ndarray of shape [n_events, n_points]) : reweight names for all points, for all events. Example name: 'rwgt_ctZ_5_ctW_3'

    Returns:
    Same arrays, but without the columns (<-> reweights) which are not necessary for the morphing procedure (e.g. SM reweight which does not follow the same naming convention)
    """

    #-- Remove coherently the incorrect weights
    indices = np.char.find(array_EFTweightIDs[0,:].astype('U'), 'rwgt_c') #Get column indices of reweight names not following the expected naming convention, from first event; 0 <-> found ; -1 <-> not found
    indices_TOP19001 = np.char.find(array_EFTweightIDs[0,:].astype('U'), 'EFTrwgt') #TOP19001 convention
    indices_TOP19001_SM = np.char.find(array_EFTweightIDs[0,:].astype('U'), 'EFTrwgt183_') #Special case: SM weight follows same conventions as other weights, but we want to remove it (it's not parsed as other benchmark weights)

    indices = indices * indices_TOP19001 #If any convention found, index = 0
    indices[indices_TOP19001_SM==0] = 1 #If TOP19001 SM weight found, we want to delete it (<-> set to non-zero)
    indices = np.asarray(indices, dtype=bool) #Convert indices to booleans (True <-> does not follow naming convention)

    array_EFTweightIDs = np.delete(array_EFTweightIDs, indices==True, axis=1) #Uncomment if reading EFTweightIDs for all events (instead of single event)
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
    Parse an array of strings, each corresponding to the name of an EFT point.

    Parameters:
    benchmarkIDs (ndarray of shape [n_points]) : array of strings, each corresponding to a separate EFT point. Example : 'rwgt_ctZ_5_ctW_3'

    Returns:
    operatorNames (ndarray of shape [n_points, n_operators]) : array of names of EFT operators, for all (n_points) EFT points
    operatorWCs (ndarray of shape [n_points, n_operators]) : array whose columns represent the WC value of each operator defining a given point, and rows represent different EFT points
    idx_SM : return index corresponding to SM point, if corresponding string is found
    """

    benchmarkIDs = np.atleast_1d(benchmarkIDs) #If a single point is given in arg, make sure it is treated as an array in the function (and not as a string)
    # print('benchmarkIDs.shape', benchmarkIDs.shape)

    prefix = 'rwgt_c' #Naming convention, strip this substring
    prefix_TOP19001 = 'EFTrwgt' #Naming convention of TOP19001

    idx_SM = -1 #Store SM index
    operatorNames = []
    operatorWCs = []
    # for idx in range(len(benchmarkIDs)):
    for idx, ID in enumerate(benchmarkIDs):

        #For each EFT point, get the list of names/WCs for all operators
        list_operatorNames = []
        list_operatorWCs = []
        # ID = benchmarkIDs[idx]

        if not ID.startswith(prefix) and (not ID.startswith(prefix_TOP19001) or ID.startswith("EFTrwgt183_")): #Every considered EFT operator is expected to start with this substring ; for others (e.g. 'rwgt_sm'), don't parse
            if ID=="rwgt_sm" or ID=="rwgt_SM" or ID.startswith("EFTrwgt183_"): idx_SM = idx #SM point found #Hard-coded naming conventions
            continue

        # print(ID)
        ID = CheckName_EFTpoint_ID(ID) #Enforce proper naming convention

        # if ID.startswith(prefix): list_keys = ID.split('_')[1:]
        if ID.startswith(prefix) or ID.startswith(prefix_TOP19001): list_keys = ID.split('_')[1:]
        else: print('Error : naming convention in benchmark ID not recognized'); exit(1)

        # print('list_keys', list_keys)
        for ikey in range(0, len(list_keys)-1, 2):
            list_operatorNames.append(list_keys[ikey].lower()) #Operator name #Force lowercase for cross-samples compatibility (using Madspin seems to make the reweight names lowercase)
            list_operatorWCs.append(float(list_keys[ikey+1])) #Operator WC

        #Append list for each EFT point
        operatorNames.append(list_operatorNames)
        operatorWCs.append(list_operatorWCs)

    #Convert list of lists into 2D array. Each list element must have an equal length
    operatorNames = np.array(operatorNames)
    operatorWCs = np.array(operatorWCs)

    if(len(operatorNames)==0): print(colors.fg.red, 'Parse_EFTpoint_IDs() -- WARNING: len(operatorNames)==0 ! benchmarkIDs: ', benchmarkIDs, colors.reset)
    return operatorNames, operatorWCs, idx_SM

# //--------------------------------------------
# //--------------------------------------------

def Translate_EFTpointID_to_WCvalues(operatorNames_sample, refPointIDs):
    '''
    Translate reference points IDs (if different from SM) into the corresponding 2D array of WC values, properly ordered following same operator order as found in the sample

    NB: only works for EFT points, not SM point (because SM point discarded in Parse_EFTpoint_IDs() )

    Return:
    orderedList_WCvalues_allPoints (array of shape [npoints,n_operators]) : array of WC values corresponding to all SMEFT points passed in arg
    '''

    refPointIDs = np.atleast_1d(refPointIDs) #If a single point is given in arg, make sure it is treated as an array in the function (and not as a string)

    if refPointIDs[0] in ["SM", "sm"]: return np.zeros(len(operatorNames_sample)) #SM point <-> all WCs are null

    operatorNames_new, operatorWCs_new, _ = Parse_EFTpoint_IDs(refPointIDs)

    orderedList_WCvalues_allPoints = []; n_skipped_points = 0
    for i_ID in range(len(refPointIDs)): #For each point ID
        if any(x in refPointIDs[i_ID] for x in ['sm','SM','EFTrwgt183']): #Don't consider SM point, cf. function description
            n_skipped_points+= 1 #Must update index accordingly
            continue
        orderedList_WCvalues = []
        for op_sample in operatorNames_sample: #For each operator found in sample
            found = False #Check whether operator in sample (needed for EFT parameterization) is found in current point's ID
            for iop_ref in range(operatorNames_new.shape[1]): #For each operator in ID
                # print('i_ID-n_skipped_points', i_ID-n_skipped_points, 'iop_ref', iop_ref)
                if str(operatorNames_new[i_ID-n_skipped_points][iop_ref]) == str(op_sample):
                    found = True
                    orderedList_WCvalues.append(operatorWCs_new[i_ID-n_skipped_points][iop_ref]) #Append WC to list at correct position (according to ordering in sample)
                elif operatorNames_new[i_ID-n_skipped_points][iop_ref] not in operatorNames_sample: print(colors.bold, colors.fg.red, 'ERROR : refPoint seems not to be compatible with operators included in this sample... (check exact naming)', colors.reset); exit(1) #Operator specified in ID not found in sample
            if found is False: orderedList_WCvalues.append(0) #If current point's ID does not include an operator found in sample (needed for parameterization), include it and set to 0

        orderedList_WCvalues_allPoints.append(orderedList_WCvalues)

    orderedList_WCvalues_allPoints = np.asarray(orderedList_WCvalues_allPoints) #List --> array

    return orderedList_WCvalues_allPoints

# //--------------------------------------------
# //--------------------------------------------

#-- If want to validate over a single operator ( e.g. 'rwgt_ctZ_3') but the DNN was training over more operators (e.g. 'rwgt_ctZ_3_ctW_0_cpQM_0_cpQ3_0_cpt_0', etc.), need to include missing operators into points names
def AddMissingOperatorsToValPointsNames(opts, list_points):

    #If a single point is given in arg, make sure it is treated as a list in the function (and not as a single string) #NB: can't use 'np.atleast_1d' here, because numpy does not treat strings properly (or should convert to object, etc.)
    if isinstance(list_points, str):
        if list_points is '': return list_points #Option 'evalPoint' may be voluntarily kept empty (see conventions)
        tmp = list_points; list_points = []; list_points.append(tmp)

    for ipt, point in enumerate(list_points): #For each point
        if point in ["SM", "sm"]: newname = "SM"
        else:
            newname = "rwgt"
            operatorNames, operatorWCs, _ = Parse_EFTpoint_IDs(point) #Get the operator names/values of point

            for opPoint in operatorNames[0,:]: #Sanity check
                if opPoint not in opts["listOperatorsParam"]: print(colors.bold, colors.fg.red, 'ERROR : validation point operator ', opPoint, ' was not used in training phase ; can not be used for evaluation', colors.reset); exit(1) #Operator specified in point's name was not used during training phase

            for opParam in opts["listOperatorsParam"]: #Add each operator used to parameterize the DNN to new name
                # print(opParam)
                newname+= "_"+opParam+"_"
                found = False
                for iOpPoint in range(operatorNames.shape[1]):
                    if opParam == operatorNames[0,iOpPoint]: #Operator found in point name -> add corresponding value next to operator name
                        found = True
                        newname+= str(operatorWCs[0,iOpPoint])
                        break
                if not found: newname+= "0" #Operator not found in point name -> add 0 next to operator name
        # print(newname)
        list_points[ipt] = newname

    return list_points

# //--------------------------------------------
# //--------------------------------------------

#-- Translate a point name (e.g. 'rwgt_ctZ_3_ctW_5_cpQM_0_cpQ3_0_cpt_0') to a proper plot legend name (e.g. 'ctZ=3,ctW=5')
def GetLegendNameEFTpoint(list_points):

    legendNames = []
    for ipt, point in enumerate(list_points):
        legname = ""
        if point in ["SM", "sm"]: legname = "SM"
        else:
            operatorNames, operatorWCs, _ = Parse_EFTpoint_IDs(point)
            # print('operatorNames', operatorNames)
            for iOpPoint in range(operatorNames.shape[1]):
                if operatorWCs[0,iOpPoint] == 0: continue #Don't show null operators
                if legname is not '': legname+= ','
                legname+= operatorNames[0,iOpPoint]+'='+str(operatorWCs[0,iOpPoint]).rstrip('0').rstrip('.')
        # print(legname)
        legendNames.append(legname)

    return legendNames

# //--------------------------------------------
# //--------------------------------------------

#-- Convert classifier response (s) <-> likelihood ratio (r)
def s_from_r(r):
    return np.clip(1.0 / (1.0 + r), 0.0, 1.0)

def r_from_s(s, epsilon=1.0e-6):
    return np.clip((1.0 - s + epsilon) / (s + epsilon), epsilon, None)

# //--------------------------------------------
# //--------------------------------------------

#-- Remove events which have w_EFT/w_SM > x from input arrays
def Remove_LargeEFTWeight_Events(w_EFT, threshold, remove_above_treshold):

    if remove_above_treshold is True: mask = np.where(np.abs(w_EFT) < np.mean(np.abs(w_EFT)) * threshold)
    else: mask = np.where(np.abs(w_EFT) > np.mean(np.abs(w_EFT)) * threshold)

    return mask

# //--------------------------------------------
# //--------------------------------------------

def Find_Components_LargeNofOperators(n_operators):
    '''
    See Utils/EFT.py --> Find_Components(): method to find all valid components entering the full amplitude becomes way too CPU-intensive for large nof operators (~ N!)
    In those cases, use this more hard-coded function: assume 'minPower_perOperator=0', 'maxPower_perOperator=2'
    '''

    #-- List to be used as template: list filled with 0, 1 entry for each operator (e.g.: 3 operators -> [0 0 0])
    template_list = []
    for iop in range(n_operators):
        template_list.append(0)

    components = [] #List of components to return
    components.append(template_list) #Append template_list (only zeros <-> corresponds to SM term)

    for iop in range(n_operators): #For each operator
        # print('iop', iop)
        tmp = template_list[:]; tmp[iop] = 2; components.append(tmp) #Append corresponding squared contribution (unique)
        tmp = template_list[:]; tmp[iop] = 1; components.append(tmp) #Append corresponding SM-EFT interference contribution (unique)

        for jop in range(iop+1,n_operators): #For all operator combinations not previously considered
            # print('jop', jop)
            tmp = template_list[:]; tmp[iop] = 1; tmp[jop] = 1; components.append(tmp) #Append corresponding EFT-EFT interference contribution (unique)

    components = sorted(set(map(tuple, components)), reverse=True) #Remove any duplicate, order

    # for i in components: print(i)
    # print(len(tmp))
    return components
