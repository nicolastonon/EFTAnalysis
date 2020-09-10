# Read ROOT files, shape/transform the data (x : input features, y : labels), compute event reweights, ...

import ROOT
from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
import numpy as np
from numpy import random
from root_numpy import root2array, tree2array, array2root, hist2array
import pandas as pd
import tensorflow
import keras
import re
from pathlib import Path
from sklearn.utils import class_weight, shuffle
from sklearn.feature_selection import RFE, SelectKBest, chi2
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, Normalizer, StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold
from tensorflow.keras import utils

from Utils.EFT import *
from Utils.ColoredPrintout import colors
from Utils.Helper import get_normalization_iqr, normalize, unison_shuffled_copies

#Filtering out manually some unimportant warnings #NB: sklean resets warning filters
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
warnings.filterwarnings("ignore", message="Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated")
warnings.filterwarnings("ignore", message="FutureWarning: in the future insert will treat boolean arrays and array-likes as boolean index instead of casting it to integer")

np.set_printoptions(threshold=np.inf) #If activated, will print full numpy arrays
# //--------------------------------------------
# //--------------------------------------------

# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
 ######   ######## ########    ########     ###    ########    ###
##    ##  ##          ##       ##     ##   ## ##      ##      ## ##
##        ##          ##       ##     ##  ##   ##     ##     ##   ##
##   #### ######      ##       ##     ## ##     ##    ##    ##     ##
##    ##  ##          ##       ##     ## #########    ##    #########
##    ##  ##          ##       ##     ## ##     ##    ##    ##     ##
 ######   ########    ##       ########  ##     ##    ##    ##     ##

########  #######  ########     ######## ########     ###    #### ##    ## #### ##    ##  ######
##       ##     ## ##     ##       ##    ##     ##   ## ##    ##  ###   ##  ##  ###   ## ##    ##
##       ##     ## ##     ##       ##    ##     ##  ##   ##   ##  ####  ##  ##  ####  ## ##
######   ##     ## ########        ##    ########  ##     ##  ##  ## ## ##  ##  ## ## ## ##   ####
##       ##     ## ##   ##         ##    ##   ##   #########  ##  ##  ####  ##  ##  #### ##    ##
##       ##     ## ##    ##        ##    ##    ##  ##     ##  ##  ##   ###  ##  ##   ### ##    ##
##        #######  ##     ##       ##    ##     ## ##     ## #### ##    ## #### ##    ##  ######
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

def Get_Data(opts, list_lumiYears, list_processClasses, list_labels, list_features, weightDir, ntuplesDir, lumiName, singleThetaName=""):
    '''
    Call sub-function to read/store/shape the data.
    '''

    #-- Get data from TFiles
    list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_SMweights_allClasses = Read_Data(opts, list_lumiYears, ntuplesDir, list_processClasses, list_labels, list_features)

    #-- For private EFT samples, get the per-event fit coefficients (for later extrapolation at any EFT point) #Central samples: empty arrays
    list_EFT_FitCoeffs_allClasses = Get_EFT_FitCoefficients_allEvents(opts, list_processClasses, list_labels, list_EFTweights_allClasses, list_EFTweightIDs_allClasses)

    #-- If the NN is parameterized on Wilson coeffs. (or training requires EFT reweighting), need to artificially extend the dataset to train on many different points in EFT phase space
    list_x_allClasses, list_weights_allClasses, list_thetas_allClasses, list_targetClass_allClasses, list_jointLR_allClasses, list_score_allClasses_allOperators = Extend_Augment_Dataset(opts, list_labels, list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses, list_SMweights_allClasses, singleThetaName)

    #-- Concatenate + reshape arrays, and modify them as needed
    x, list_weights_allClasses, thetas_allClasses, targetClass_allClasses, jointLR_allClasses, scores_allClasses_eachOperator, list_nentries_class = Shape_Data(opts, list_x_allClasses, list_weights_allClasses, list_thetas_allClasses, list_targetClass_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses, list_jointLR_allClasses, list_score_allClasses_allOperators, singleThetaName)

    #-- Define 'physical event weights' (for plotting, ...) and 'training weights' (rescaled arbitrarily to improve the training procedure)
    LearningWeights_allClasses, PhysicalWeights_allClasses = Get_Events_Weights(opts, list_processClasses, list_labels, list_weights_allClasses, targetClass_allClasses)

    #-- For parameterized NN, need to count theory parameters as additional inputs. Also, from there on, different classes correspond to SM vs EFT (not to different physics processes)
    list_labels, list_features = Update_Lists(opts, list_labels, list_features)

    #-- Define the targets 'y' (according to which the classification/regression is performed). Also keep track of the process class indices of all events ('y_process')
    y, y_process, x = Get_Targets(opts, list_features, list_processClasses, list_nentries_class, targetClass_allClasses, jointLR_allClasses, scores_allClasses_eachOperator, x)

    #-- Sanitize data
    x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses = Sanitize_Data(opts, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, singleThetaName)

    if singleThetaName is not "": return x, y, y_process, PhysicalWeights_allClasses, list_labels, list_features #Only for validation, don't need to split between train/test data

    #-- No validation dataset required --> shuffle and split the data between training / testing datasets only
    if opts["splitTrainValTestData"][1] == 0. and opts["nEventsTot_val"] < 0:
        x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test = Train_Test_Split(opts, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses)
        x_val = x_test; y_val = y_test; y_process_val = y_process_test; PhysicalWeights_val = PhysicalWeights_test; LearningWeights_val = LearningWeights_test #Take testing dataset as validation dataset
    else: x_train, x_val, x_test, y_train, y_val, y_test, y_process_train, y_process_val, y_process_test, PhysicalWeights_train, PhysicalWeights_val, PhysicalWeights_test, LearningWeights_train, LearningWeights_val, LearningWeights_test = Train_Val_Test_Split(opts, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses) #Split data into train / val / test datasets

    #FIXME
    # mask_largeEFTweights = Remove_LargeEFTWeight_Events(LearningWeights_train, 30)
    # x_train = x_train[mask_largeEFTweights]
    # y_train = y_train[mask_largeEFTweights]
    # y_process_train = y_process_train[mask_largeEFTweights]
    # PhysicalWeights_train = PhysicalWeights_train[mask_largeEFTweights]
    # LearningWeights_train = LearningWeights_train[mask_largeEFTweights]

    #-- Get rescaling parameters for each input feature, given to first NN layer to normalize features -- derived from train data alone
    transfType='quantile' #'quantile', 'range', 'gauss'. Defines the transformation applied to normalize input data.
    xTrainRescaled, shifts, scales = Transform_Inputs(weightDir, x_train, list_features, lumiName, opts["parameterizedNN"], transfType=transfType)

    print(colors.fg.lightblue, "\n===========")
    print("-- Total nof events : " + str(x.shape[0]))
    print("-- Will use " + str(x_train.shape[0]) + " training events !")
    print("-- Will use " + str(x_val.shape[0]) + " validation events !")
    print("-- Will use " + str(x_test.shape[0]) + " testing events !")
    print("===========\n", colors.reset)

    return x_train, x_val, x_test, y_train, y_val, y_test, y_process_train, y_process_val, y_process_test, PhysicalWeights_train, PhysicalWeights_val, PhysicalWeights_test, LearningWeights_train, LearningWeights_val, LearningWeights_test, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, shifts, scales, xTrainRescaled, list_labels, list_features

# //--------------------------------------------
# //--------------------------------------------

 #####  ######   ##   #####     #####    ##   #####   ##
 #    # #       #  #  #    #    #    #  #  #    #    #  #
 #    # #####  #    # #    #    #    # #    #   #   #    #
 #####  #      ###### #    #    #    # ######   #   ######
 #   #  #      #    # #    #    #    # #    #   #   #    #
 #    # ###### #    # #####     #####  #    #   #   #    #

def Read_Data(opts, list_lumiYears, ntuplesDir, list_processClasses, list_labels, list_features):
    '''
    Read the data from ROOT files and store it in np arrays.
    '''

    testdirpath = Path(ntuplesDir)
    if not testdirpath.is_dir():
        print('Ntuple dir. '+ntuplesDir+' not found ! Abort !')
        exit(1)

    cuts = opts["cuts"]

    list_x_allClasses = [] #List of x-arrays storing the values of input features (for all events, all considered years, all physics processes in a given class) -- 1 array per process class
    list_weights_allClasses = [] #Idem for event weights
    list_EFTweights_allClasses = [] #Idem for EFT reweights
    list_EFTweightIDs_allClasses = [] #Idem for EFT reweights IDs
    list_SMweights_allClasses = [] #Idem for SM reweight
    for procClass, label in zip(list_processClasses, list_labels): #Loop on classes of physics processes (e.g. 'Signal','Backgrounds') #NB: can not use predefined keyword 'class'

        print('\n', colors.fg.purple, colors.underline, '* Class :', label, colors.reset, '\n')

        list_x_proc = [] #List of x-arrays storing the values of input features for all events, for all considered years -- 1 array per physics process (sample)
        list_weights_proc = [] #Idem for central event weights
        list_EFTweights_proc = []; list_EFTweightIDs_proc = [] #Idem for EFT reweights
        list_SMweights_proc = [] #Store SM weights for easy access (for private samples only)
        for iproc, process in enumerate(procClass): #Loop on physics processes (samples)

            cuts_tmp = cuts
            isNPLsample = False
            if "npl" in process.lower() or "dy" in process.lower() or "ttbar" in process.lower():
                isNPLsample = True
                cuts_tmp = cuts.replace("_SR", "_SRFake"); cuts_tmp = cuts.replace("_SR", "_CRFake") #NPL sample <-> different category flags

            isPrivMCsample = False
            if ("PrivMC" in process and "PrivMC" not in label) or ("PrivMC" in label and "PrivMC" not in process):
                print('process', process); print('label', label)
                print('\n', colors.fg.red, 'Error : keyword \'PrivMC\' must be present both in process and class names, or not at all', colors.reset) #Avoid ambiguities

            elif "PrivMC" in process and "PrivMC" in label: #Check whether EFT reweights should be looked for
                isPrivMCsample = True
                # print('process', process, ', label', label, ' --> isPrivMCsample')

            isPureEFT = False
            if "PrivMC" in label and "_c" in label: isPureEFT = True #Naming convention associated with pure-EFT sample --> no SM point

            print(colors.fg.pink, '* Process :', colors.reset, process)

            for iyear in range(len(list_lumiYears)): #Concatenate data (input features, event weights) for all considered years

                filepath = ntuplesDir + list_lumiYears[iyear] + '/' + process + '.root'
                if not Path(filepath).is_file():
                    print('File '+filepath+' not found ! Abort !')
                    exit(1)

                # print(colors.fg.lightgrey, '* Opening file:', colors.reset, ' ', filepath)
                file = TFile.Open(filepath)
                # tree = file.Get('result')
                tree = file.Get(opts["TTree"])

                nevents = None
                if opts["parameterizedNN"] == False and opts["maxEventsPerClass"] > 0: nevents = opts["maxEventsPerClass"]
                wname_tmp = 'eventWeight' #By default (for my ntuples), read this variable for event weight
                if opts["TTree"] != 'result': wname_tmp = opts["eventWeightName"]
                print(colors.fg.lightgrey, 'Opened file:', colors.reset, filepath, '(Total nof entries:', tree2array(tree, branches=wname_tmp, selection=cuts_tmp).shape[0], 'entries)') #Dummy variable, just to read the nof entries
                if nevents is not None: print('(---> Will consider at most ' + str(nevents) + ' entries [<-> maxEventsPerClass])')
                print('\n\n')

                # list_x_proc.append(tree2array(tree, branches=list_features, selection=cuts_tmp)) #Store values of input features into array, append to list

                #-- root_numpy 'tree2array' function returns numpy structured array : 1D array whose length equals the nof events, and each element is a structure with multiple fields (1 per feature)
                #For manipulation, it is easier to convert structured arrays obtained in this way into regular numpy arrays (e.g. x will be 2D and have shape (n_events, n_features) )
                # x_tmp = tree2array(tree, branches=list_features, selection=cuts_tmp, start=0,stop=nevents)
                x_tmp = tree2array(tree, branches=list_features, selection=cuts_tmp); x_tmp = x_tmp[:nevents]
                x_tmp = np.column_stack([x_tmp[name] for name in x_tmp.dtype.names]) #1D --> 2D
                x_tmp = x_tmp.astype(np.float32) #Convert all to floats
                list_x_proc.append(x_tmp) #Append features to list

                #-- Store event weights into array, append to list
                # if opts["eventWeightName"] != '': list_weights_proc.append(tree2array(tree, branches=opts["eventWeightName"], selection=cuts_tmp, start=0,stop=nevents))
                # elif isPureEFT is True: list_weights_proc.append(tree2array(tree, branches="eventWeight", selection=cuts_tmp, start=0,stop=nevents)) #For pure-EFT samples, weights are non physical. Just use the baseline MG weight, don't multiply by lumi*xsec or divide by SWE... rescaled for training anyway (and validation weights are non-physical)
                # else: list_weights_proc.append(tree2array(tree, branches="eventWeight*eventMCFactor", selection=cuts_tmp, start=0,stop=nevents))
                if opts["eventWeightName"] != '': list_weights_proc.append(tree2array(tree, branches=opts["eventWeightName"], selection=cuts_tmp)[:nevents])
                elif isPureEFT is True: list_weights_proc.append(tree2array(tree, branches="eventWeight", selection=cuts_tmp)[:nevents]) #For pure-EFT samples, weights are non physical. Just use the baseline MG weight, don't multiply by lumi*xsec or divide by SWE... rescaled for training anyway (and validation weights are non-physical)
                else: list_weights_proc.append(tree2array(tree, branches="eventWeight*eventMCFactor", selection=cuts_tmp)[:nevents])

            if isPrivMCsample: #For private MC samples, get the EFT reweights (properly normalized) and their IDs
                EFTweights_proc_tmp, EFTweightIDs_proc_tmp, SMweights_proc_tmp = Read_Data_EFT_File(opts, list_lumiYears, list_weights_proc, ntuplesDir, process, cuts_tmp, isPureEFT, iproc, nevents)
                list_EFTweights_proc.append(EFTweights_proc_tmp); list_EFTweightIDs_proc.append(EFTweightIDs_proc_tmp); list_SMweights_proc.append(SMweights_proc_tmp)

        #-- Concatenate the different arrays (for all years, processes) corresponding to a single class of process, and append them to their lists --> 1 single array per process class
        list_x_allClasses.append(np.concatenate(list_x_proc))

        #-- Concatenate other arrays
        if "PrivMC" in label and isPureEFT is False: #For pure-EFT samples, don't care about EFT reweights (only use baseline weights)
            list_weights_allClasses.append(np.concatenate(list_EFTweights_proc,axis=0)[:,0]) #For private EFT samples, only the weights in 'list_EFTweights_proc' make sense (not those in 'list_weights_proc') ; but need single value per event --> use element corresponding to baseline (first column)
            list_EFTweights_allClasses.append(np.concatenate(list_EFTweights_proc,axis=0))
            list_EFTweightIDs_allClasses.append(np.concatenate(list_EFTweightIDs_proc,axis=0))
            list_SMweights_allClasses.append(np.concatenate(list_SMweights_proc,axis=0))
        else:
            weights_tmp = np.concatenate(list_weights_proc)
            nentries = weights_tmp.shape[0]
            list_weights_allClasses.append(weights_tmp)
            list_EFTweights_allClasses.append(np.full(shape=nentries, fill_value=-999)); list_EFTweightIDs_allClasses.append(np.full(shape=nentries, fill_value='-999', dtype=object)); list_SMweights_allClasses.append(np.full(shape=nentries, fill_value=-999)); #Append empty arrays (1 row per event) for proper ordering

    print('\n\n')

    return list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_SMweights_allClasses

# //--------------------------------------------
# //--------------------------------------------

def Read_Data_EFT_File(opts, list_lumiYears, list_weights_proc, ntuplesDir, process, cuts, isPureEFT, iproc, nevents):
    '''
    For private MC (EFT) samples, retrieve the EFT reweights and their IDs. Directly normalize properly the EFT weights.
    '''

    list_EFTweights_proc = []
    list_EFTweightIDs_proc = []
    list_SMweights_proc = []
    for iyear in range(len(list_lumiYears)):

        filepath = ntuplesDir + list_lumiYears[iyear] + '/' + process + '.root'
        if not Path(filepath).is_file():
            print('File '+filepath+' not found ! Abort !')
            exit(1)

        file = TFile.Open(filepath)
        # tree = file.Get('result')
        tree = file.Get(opts["TTree"])
        # print(colors.fg.lightgrey, '* Opened file:', colors.reset, filepath, '(', tree2array(tree, branches="eventWeight", selection=cuts).shape[0], 'entries )') #Dummy variable, just to read the nof entries

        #NB: don't want to store EFT reweights IDs for all events (always same names)
        #NB : assume that EFT reweights IDs are exactly the same, in the same order, for all considered years !
        #Only read a single event, store IDs once per private MC process (<-> single process class)
        # if iyear==0:
        #     array_EFTweightIDs_proc = tree2array(tree, start=0, stop=1, branches="mc_EFTweightIDs", selection=cuts)
        #     array_EFTweightIDs_proc = array_EFTweightIDs_proc[0] #Reshape array
        # elif np.array_equal(array_EFTweightIDs_proc, tree2array(tree, start=0, stop=1, branches="mc_EFTweightIDs", selection=cuts)) == False: #If IDs are already stored, simply verify that IDs for other years match
        #     print('\n', colors.fg.red, 'Error : EFT reweight IDs do not seem to match between the different years for sample:', colors.reset, process)
        #     exit(1)

        #Array of reweights from MG
        if iproc == 0: weightsProc = list_weights_proc[iyear] #If single process in current class (or first element) -> elements positions in list_weights_proc only depend on year
        else: weightsProc = list_weights_proc[iproc*len(list_lumiYears) + iyear] #Else, if multiple processes in current class, must also account for previous processes in list and update index position

        #Get the EFT reweights IDs
        # array_EFTweightIDs_proc = np.stack(tree2array(tree, branches="mc_EFTweightIDs", selection=cuts, start=0,stop=nevents)) #stack : array of arrays -> 2d array
        array_EFTweightIDs_proc = np.stack(tree2array(tree, branches="mc_EFTweightIDs", selection=cuts)[:nevents]) #stack : array of arrays -> 2d array

        #Get the EFT reweights, and normalization factor (because will multiply by baseline weight as a trick to apply the SFs to all weights ; must then divide by baseline weight 'weightMENominal')
        # array_EFTweights_proc = np.stack(tree2array(tree, branches="mc_EFTweights", selection=cuts, start=0,stop=nevents)) #stack : array of arrays -> 2d array
        # normWeights_proc = tree2array(tree, branches="weightMENominal", selection=cuts, start=0,stop=nevents) #Normalization factors
        array_EFTweights_proc = np.stack(tree2array(tree, branches="mc_EFTweights", selection=cuts)[:nevents]) #stack : array of arrays -> 2d array
        normWeights_proc = tree2array(tree, branches="weightMENominal", selection=cuts)[:nevents] #Normalization factors

        #Get the sums of weights (before any preselection) corresponding to each EFT point
        hist = file.Get("EFT_SumWeights") #Sums of weights for each EFT reweight is stored in histogram, read it
        array_EFT_SWE_proc = hist2array(hist) #Store into array
        array_EFT_SWE_proc = array_EFT_SWE_proc[0:array_EFTweights_proc.shape[1]] #Only need the SWE values for the considered reweight points

        #Store the MG reweight value for the SM point
        idx_SM = -1
        if not isPureEFT:
            for i in range(len(array_EFTweightIDs_proc[0])):
                if array_EFTweightIDs_proc[0][i] == "rwgt_sm" or array_EFTweightIDs_proc[0][i] == "rwgt_SM": idx_SM = i; break

        array_SMweights_proc = np.array([])
        if idx_SM != -1: #Get weights at SM point (and normalize them, cf. below)
            array_SMweights_proc = array_EFTweights_proc[:,i]; #print('idx_SM = ',idx_SM)
            array_SMweights_proc = array_SMweights_proc * weightsProc
            array_SMweights_proc = np.divide(array_SMweights_proc, normWeights_proc)
            array_SMweights_proc = np.divide(array_SMweights_proc, array_EFT_SWE_proc[idx_SM])
        elif not isPureEFT: print(colors.fg.red, 'ERROR: from naming convention, infer that this is a SM+EFT sample. However, the benchmark SM point from MG was not found (needed for proper normalization since the sample xsec is expected to correspond to SM. Abort ! (If this is not the desired behaviour, adapt the code...)', colors.reset); exit(1)

        #-- Normalize EFT weights: rwgt = rwgt_from_MG * (centralEventWeight/weightMENominal) / SWE_SM
        #-- 1st term is the MG reweight for a given point idx ; 2nd term corresponds to SF * L * xsec(SM) ; 3rd term is needed for proper normalization, it is the sum of weights before preselection at the SM point (consistent with using the xsec at SM point)
        array_EFTweights_proc = array_EFTweights_proc * weightsProc[:,None] #rwgt_from_MG * centralEventWeight
        array_EFTweights_proc = np.divide(array_EFTweights_proc, normWeights_proc[:,None]) # ... / weightMENominal
        SWE_SM = array_EFT_SWE_proc[idx_SM] # ... / SWE_SM
        if SWE_SM == 0: SWE_SM = 1 #For pure-EFT samples, SM=0 -> Set to 1 (weights are unphysical anyway)
        array_EFTweights_proc = np.divide(array_EFTweights_proc, SWE_SM)

        #Manually find and remove all weights with unproper naming conventions (for example 'rwgt_1' nominal weight is included by default by MG)
        array_EFTweights_proc, array_EFTweightIDs_proc = Remove_Unnecessary_EFTweights(array_EFTweights_proc, array_EFTweightIDs_proc)

        list_EFTweights_proc.append(array_EFTweights_proc) #Append array of EFT reweights (for given year) to list
        list_EFTweightIDs_proc.append(array_EFTweightIDs_proc) #Append array of EFT reweights IDs (for given year) to list
        list_SMweights_proc.append(array_SMweights_proc) #Append array of SM reweights (for given year) to list

    # return list_EFTweights_proc, list_EFTweightIDs_proc, list_SMweights_proc
    return np.concatenate(list_EFTweights_proc,axis=0), np.concatenate(list_EFTweightIDs_proc,axis=0), np.concatenate(list_SMweights_proc,axis=0)

# //--------------------------------------------
# //--------------------------------------------

 ###### ###### #####     ####   ####  ###### ###### ######  ####
 #      #        #      #    # #    # #      #      #      #
 #####  #####    #      #      #    # #####  #####  #####   ####
 #      #        #      #      #    # #      #      #           # ###
 #      #        #      #    # #    # #      #      #      #    # ###
 ###### #        #       ####   ####  ###### #      #       ####  ###

def Get_EFT_FitCoefficients_allEvents(opts, list_processClasses, list_labels, list_EFTweights_allClasses, list_EFTweightIDs_allClasses):
    '''
    Get the 'fit coefficients' (1 per component entering the squared matrix element of the EFT process) A satisfying : A.T=w, with w the benchmark weights and T the matrix of 'effective WCs' corresponding to the benchmark points.
    Once these fit coefficients are extracted for a given event, they can be used to extrapolate the event weight at any new EFT point.
    '''

    # for iclass in range(len(list_EFTweights_allClasses)): list_EFTweights_allClasses[iclass]=list_EFTweights_allClasses[iclass][:5]; list_EFTweightIDs_allClasses[iclass]=list_EFTweightIDs_allClasses[iclass][:5] #For debugging

    list_EFT_FitCoeffs_allClasses = []
    for iclass in range(len(list_processClasses)): #Loop on classes of physics processes

        # if "PrivMC" in list_processClasses[iclass][0] and "PrivMC" in list_labels[iclass]: #Check whether EFT reweights should be looked for
        if "PrivMC" in list_processClasses[iclass][0] and "PrivMC" in list_labels[iclass] and "_c" not in list_processClasses[iclass][0]: #Check whether EFT reweights should be looked for

            operatorNames, operatorWCs, _ = Parse_EFTpoint_IDs(list_EFTweightIDs_allClasses[iclass][0]) #Get the lists of operator names and WC values for this process #NB: assumes that they are identical for all events in this process
            n_components, components = Find_Components(operatorNames[0]) #Determine the components required to parameterize the event weight #NB: assumes that they are identical for all events in this process
            effWC_components = Get_EffectiveWC_eachComponent(n_components, components, operatorWCs) #Determine the 'effective WC' values associated with each component, for each benchmark point
            fit_coeffs = Get_FitCoefficients(effWC_components, list_EFTweights_allClasses[iclass], n_components, components) #Determine the fit coefficients of the events, based on the benchmark weights and 'effective WC' values

            #Sanity check: for parameterized strategies, all operators selected by user must be found in sample parametrization (<-> in MG reweight names)
            if opts["parameterizedNN"] is True:
                for op in opts["listOperatorsParam"]:
                    if op not in operatorNames[0,:]: print(colors.fg.red, 'ERROR: operator', op, 'selected in option [listOperatorsParam] not found in event weight parametrization. The operators found in sample are:', operatorNames[0,:],'. Abort !', colors.reset); exit(1)

            #-- Debugging
            # print('WCs',effWC_components.shape,  effWC_components[:5])
            # print('coeffs', fit_coeffs.shape, fit_coeffs[:5])
            # print('w0', np.dot(effWC_components[0],np.transpose(fit_coeffs[0])))
            # print('w1 ev1', np.dot(effWC_components[1],np.transpose(fit_coeffs[1])))
            # print('w2 ev2', np.dot(effWC_components[2],np.transpose(fit_coeffs[2])))
            # print('w2 ev2', np.dot(np.transpose(fit_coeffs[2]), effWC_components[2]))
            # print('bench', list_EFTweights_allClasses[iclass][:3])
            # exit(1)

            list_EFT_FitCoeffs_allClasses.append(fit_coeffs) #Append fit coeffs to list

        else: list_EFT_FitCoeffs_allClasses.append(np.full(shape=len(list_EFTweightIDs_allClasses[iclass]), fill_value=0)) #If process is not EFT, fill with dummy value for now (maintain ordering)

    return list_EFT_FitCoeffs_allClasses

# //--------------------------------------------
# //--------------------------------------------

  ####  #    #   ##   #####  ######
 #      #    #  #  #  #    # #
  ####  ###### #    # #    # #####
      # #    # ###### #####  #
 #    # #    # #    # #      #
  ####  #    # #    # #      ######

def Shape_Data(opts, list_x_arrays_allClasses, list_weights_allClasses, list_thetas_allClasses, list_targetClass_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses, list_jointLR_allClasses, list_score_allClasses_allOperators, singleThetaName=""):
    '''
    Properly shape the arrays and concatenate them (for all years, processes, etc.).

    NB: nominal weights are still returned as a list (want to retain process class info), and get concatenated in dedicated function Get_Events_Weights()
    '''

    #--- Get nof entries for each class
    list_nentries_class = []
    for iclass in range(len(list_x_arrays_allClasses)): list_nentries_class.append(len(list_x_arrays_allClasses[iclass]))

    maxEvents = opts["maxEvents"]

    #--- Max nof events for train/test phases (not applied for parameterized NN, since total nof events depends on other options)
    if maxEvents != -1 and opts["parameterizedNN"] == False:

        #Skim each process class (keep only maxEvents events) ; skim all relevant arrays coherently
        for iclass in range(len(list_x_arrays_allClasses)):

            if list_nentries_class[iclass] > maxEvents: #NB: could not make it work with lists... !

                #-- If only consider part of the process class data, shuffle events, so that the resulting dataset is representative of the event proportions of each process within the class (else, it could happen that e.g. only events from the first process get considered)
                #NB: only arrays potentially used for training (x, weight, theta, ...) must be shuffled simultaneously
                # if opts["strategy"] is "classifier" or opts["strategy"] is "regressor" or (np.all(list_EFTweights_allClasses[iclass] == -999)): list_tmp=[list_x_arrays_allClasses,list_weights_allClasses]
                # else: list_tmp=[list_x_arrays_allClasses,list_weights_allClasses,list_thetas_allClasses,list_targetClass_allClasses]
                if opts["strategy"] is "classifier" or opts["strategy"] is "regressor" or (np.all(list_EFTweights_allClasses[iclass] == -999)):
                    list_x_arrays_allClasses[iclass],list_weights_allClasses[iclass] = unison_shuffled_copies(list_x_arrays_allClasses[iclass],list_weights_allClasses[iclass]) # To verify
                    # list_tmp=[list_x_arrays_allClasses[iclass],list_weights_allClasses[iclass]]
                    list_x_arrays_allClasses[iclass] = list_x_arrays_allClasses[iclass][0:maxEvents]
                    list_weights_allClasses[iclass] = list_weights_allClasses[iclass][0:maxEvents]
                else:
                    list_x_arrays_allClasses[iclass],list_weights_allClasses[iclass],list_thetas_allClasses[iclass],list_targetClass_allClasses[iclass] = unison_shuffled_copies(list_x_arrays_allClasses[iclass],list_weights_allClasses[iclass],list_thetas_allClasses[iclass],list_targetClass_allClasses[iclass])
                    # list_tmp=[list_x_arrays_allClasses[iclass],list_weights_allClasses[iclass],list_thetas_allClasses[iclass],list_targetClass_allClasses[iclass]]
                    list_x_arrays_allClasses[iclass] = list_x_arrays_allClasses[iclass][0:maxEvents]
                    list_weights_allClasses[iclass] = list_weights_allClasses[iclass][0:maxEvents]
                    list_thetas_allClasses[iclass] = list_thetas_allClasses[iclass][0:maxEvents]
                    list_targetClass_allClasses[iclass] = list_targetClass_allClasses[iclass][0:maxEvents]

                list_nentries_class[iclass] = maxEvents

                # for il in range(len(list_tmp)): #NB: can't do 'for l in list_tmp', since then would create copies rather than modify list elements directly
                #     if len(list_tmp[il]) > 0:
                #         list_tmp[il] = list_tmp[il][0:maxEvents]
                #         print(len(list_tmp[il]))
                #         print(list_tmp[il][:5])

            '''
            #-- Only needed if want to propagate EFT arrays further through the code, and therefore need to shape them properly even when they are empty. Not needed for now
            for jclass in range(len(list_x_arrays_allClasses)):
                if iclass is not jclass:

                    #-- For central samples (dummy array of EFT weights is 1D instead of 2D), reshape to same nof columns (<-> nof reweight points) as in EFT samples, so that this array could be manipulated identically for all classes
                    if list_EFTweights_allClasses[iclass].ndim is 2 and list_EFTweights_allClasses[jclass].ndim is 1:
                        list_EFTweights_allClasses[jclass] = np.full(shape=(list_EFTweights_allClasses[iclass].shape[0],list_EFTweights_allClasses[iclass].shape[1]), fill_value=0.)
                        list_EFTweightIDs_allClasses[jclass] = np.full(shape=(list_EFTweightIDs_allClasses[iclass].shape[0],list_EFTweightIDs_allClasses[iclass].shape[1]), fill_value='-1', dtype=object)
                        list_EFT_FitCoeffs_allClasses[jclass] = np.full(shape=(list_EFT_FitCoeffs_allClasses[iclass].shape[0],list_EFT_FitCoeffs_allClasses[iclass].shape[1]), fill_value=0.)
                        list_thetas_allClasses[jclass] = np.full(shape=(list_thetas_allClasses[iclass].shape[0],list_EFT_FitCoeffs_allClasses[iclass].shape[1]), fill_value=0.)

                    #-- Also, different EFT samples may have different numbers of columns (<-> nof benchmark reweights) --> Force them to have same nof columns, for manipulation
                    elif list_EFTweights_allClasses[iclass].ndim is 2 and list_EFTweights_allClasses[jclass].ndim is 2 and list_EFTweights_allClasses[iclass].shape[1] > list_EFTweights_allClasses[jclass].shape[1]:
                        tmp = np.zeros((list_EFTweights_allClasses[jclass].shape[0],list_EFTweights_allClasses[iclass].shape[1]))
                        tmp[:,:-1] = list_EFTweights_allClasses[jclass]
                        list_EFTweights_allClasses[jclass] = tmp

                        tmp = np.full(shape=(list_EFTweightIDs_allClasses[jclass].shape[0],list_EFTweightIDs_allClasses[iclass].shape[1]), fill_value='-1', dtype=object)
                        tmp[:,:-1] = list_EFTweightIDs_allClasses[jclass]
                        list_EFTweightIDs_allClasses[jclass] = tmp

                        tmp = np.zeros(shape=(list_EFT_FitCoeffs_allClasses[jclass].shape[0],list_EFT_FitCoeffs_allClasses[iclass].shape[1]) )
                        tmp[:,:-1] = list_EFT_FitCoeffs_allClasses[jclass]
                        list_EFT_FitCoeffs_allClasses[jclass] = tmp

                        tmp = np.zeros(shape=(list_thetas_allClasses[jclass].shape[0],list_thetas_allClasses[iclass].shape[1]) )
                        tmp[:,:-1] = list_thetas_allClasses[jclass]
                        list_thetas_allClasses[jclass] = tmp
            '''

    #Default dummy arrays and lists
    jointLR_allClasses = np.array([]); scores_allClasses_eachOperator = [] #Scores are stored in a list, 1 element per component
    thetas_allClasses = np.zeros((1,1)); targetClass_allClasses = thetas_allClasses

    #-- Transform list of arrays --> single concatenated array
    x = np.concatenate(list_x_arrays_allClasses, 0)
    # EFTweights_allClasses = np.concatenate(list_EFTweights_allClasses, 0); EFTweightIDs_allClasses = np.concatenate(list_EFTweightIDs_allClasses, 0); EFT_FitCoeffs_allClasses = np.concatenate(list_EFT_FitCoeffs_allClasses, 0)

    if len(list_thetas_allClasses)is 0 and opts["makeValPlotsOnly"] is False and opts["parameterizedNN"] is True: print('Warning: len(list_thetas_allClasses)==0...')

    # if opts["strategy"] in ["CARL_singlePoint", "CARL_multiclass"]: #-- can not parameterize CARL_multiclass as I did, else NN relies ~ only on WC values
    if opts["strategy"] == "CARL_singlePoint": targetClass_allClasses = np.concatenate(list_targetClass_allClasses, axis=0)

    #-- parameterized NN: pass the values of the WCs as input features
    elif opts["parameterizedNN"]==True:
        thetas_allClasses = np.concatenate(list_thetas_allClasses, 0)
        targetClass_allClasses = np.concatenate(list_targetClass_allClasses, 0)

        if opts["strategy"] in ["ROLR", "RASCAL"]:
            jointLR_allClasses = np.concatenate(list_jointLR_allClasses, 0)
            if opts["strategy"] is "RASCAL": #Concatenate different classes (first dim) together ! But retain the ordering ot the operator components (second dim) and events (third dim)
                scores_allClasses_eachOperator = np.concatenate(np.array(list_score_allClasses_allOperators), 0)
                if len(scores_allClasses_eachOperator)==1: scores_allClasses_eachOperator = np.squeeze(scores_allClasses_eachOperator, 0) #Squeeze first dimension (<-> if there is a single proc class there was no concatenation --> remove this useless dimension)
                scores_allClasses_eachOperator = scores_allClasses_eachOperator.T #problems for single operator ?

        #-- Append WC values to input features, etc.
        # if singleThetaName is "SM" or singleThetaName is "sm": #If making validation plots for SM, still need to append WC values to x (all set to 0)
        #     tmp = np.zeros((len(x),len(opts["listOperatorsParam"])))
        #     x = np.append(x, tmp, axis=1)
        if opts["strategy"] is not "classifier": #Trick: in case we run StandaloneValidation code on classifier (not parameterized), want to evaluate on SMEFT sample but without including WCs as inputs !
            if singleThetaName is not "": #For standalone val in multiclass, need to evaluate on a single operator (<-> single class), but retain all operators used during training
                # print('singleThetaName', singleThetaName); print('evalPoint', opts['evalPoint'])

                tmp = thetas_allClasses

                #The lines below have 1 purpose: whatever the operators considered in 'singleThetaName' (<-> validation) or 'list_EFTweightIDs_allClasses' (<-> EFT parametrization in sample), we only want to keep the operators which were included in the training phase, for which we need to provide WC values
                #NB: really not sure about the correctness here...
                operatorNames_tmp, operatorWCs_tmp, _ = Parse_EFTpoint_IDs(list_EFTweightIDs_allClasses[0][0][0])
                indices_opToRemove = []
                for iOpSample, opSample in enumerate(operatorNames_tmp[0]): #NB: operatorNames_tmp is 2D
                    # print('opSample', opSample); print('opts[evalPoint]', opts['evalPoint']), print('opSample in opts[evalPoint]', str(opSample) in str(opts['evalPoint']))
                    if str(opts['evalPoint']) not in ['','SM','sm'] and opSample not in str(opts['evalPoint']): indices_opToRemove.append(iOpSample)
                    elif opts['evalPoint'] in ['','SM','sm'] and opSample not in opts['listOperatorsParam']: indices_opToRemove.append(iOpSample)
                # print(indices_opToRemove)
                if len(indices_opToRemove)>0: tmp = np.delete(thetas_allClasses, indices_opToRemove, axis=1)

                x = np.append(x, tmp, axis=1) #NB -- can remove input WC here

            else: #Otherwise, only consider selected operators
                #'thetas_allClasses' has as many columns as there are EFT operators generated in the sample (needed for extraction of fit coefficients from benchmark weights). But from there, only want to retain EFT operators which the NN will get trained on --> Only parameterize NN on such operators, not the others (not used)
                theta_tmp = thetas_allClasses[:, ~np.all(thetas_allClasses==0, axis=0)] #Only keep columns (operators) which were activated by the user #'~' is negation
                if opts["strategy"] is "CARL_multiclass": targetClass_allClasses = targetClass_allClasses[:, ~np.all(targetClass_allClasses==0, axis=0)] #Idem (only needed for multiclass, where class is encoded in multiple columns)
                targetClass_allClasses = np.squeeze(targetClass_allClasses) #If 2D with single column, squeeze into 1D array
                x = np.append(x, theta_tmp, axis=1) #NB -- can remove input WC here

    # print(x.shape)

    return x, list_weights_allClasses, thetas_allClasses, targetClass_allClasses, jointLR_allClasses, scores_allClasses_eachOperator, list_nentries_class

# //--------------------------------------------
# //--------------------------------------------

 #    # ###### #  ####  #    # #####  ####
 #    # #      # #    # #    #   #   #
 #    # #####  # #      ######   #    ####
 # ## # #      # #  ### #    #   #        #
 ##  ## #      # #    # #    #   #   #    #
 #    # ###### #  ####  #    #   #    ####

def Get_Events_Weights(opts, list_processClasses, list_labels, list_weights_allClasses, targetClass_allClasses):
    '''
    Compute and apply weights to training dataset to balance the training.
    There are 3 possibilities:
    1) Classification between physics processes --> rescale each process to same total training weight
    2) Classification between SM and single EFT point ('CARL_singlePoint' strategy) --> rescale each hypothesis to same total training weight (but merge all process classes together)
    3) Parametrized classifier of regressor --> set all training weights to 1 (because samples were already unweighted, to draw events according to their weights)

    NB: using absolute event weights only. Keras can deal with negative weights (inverse impact on the loss function, depending on which the neurons' weights get updated), but probably does not make sense: e.g. if a signal event was correctly classified, we don't want a negative weight to treat that as an error ! But as a consequence of only using abs(w), we are biasing the training phase space --> sub-optimal performance.
    '''

    parameterizedNN = opts["parameterizedNN"]

    #Duplicate list of weight arrays, storing absolute weights (can't handle negative weights in training)
    list_weights_allClasses_abs = []
    for weights_class in list_weights_allClasses: list_weights_allClasses_abs.append(np.absolute(weights_class))

    #Also create corresponding array of physical event weights, to get correct plots, etc.
    PhysicalWeights_allClasses = np.concatenate(list_weights_allClasses, 0)

    #Create array for learning weights (weights used for training)
    list_LearningWeights_allClasses = []

    if opts["balancedClasses"] is False: return PhysicalWeights_allClasses, PhysicalWeights_allClasses #Choose not to balance classes --> Learning weights are taken as abs(physical weights)

    if opts["strategy"] is "CARL_singlePoint": #For this (non-parameterized strategy), there is only 1 process/class (SMEFT sample) --> Check 'y_target' to determine whether an event is signal or bkg #Apply SFs at the end of function

        allweights = np.concatenate(list_weights_allClasses_abs)
        yield_SM = allweights[targetClass_allClasses==0].sum()
        yield_EFT = allweights[targetClass_allClasses==1].sum()
        SF_SM = 100. / yield_SM
        SF_EFT = 100. / yield_EFT

        print(colors.ital, '* SM:', colors.reset)
        print('Default yield = ', float('%.4g' % yield_SM))
        print('Rescaling factor = ', float('%.3g' % SF_SM))
        print('==> Rescaled yield :', float('%.2g' % (yield_SM*SF_SM)), '\n')
        print(colors.ital, '* EFT:', colors.reset)
        print('Default yield = ', float('%.4g' % yield_EFT))
        print('Rescaling factor = ', float('%.3g' % SF_EFT))
        print('==> Rescaled yield :', float('%.2g' % (yield_EFT*SF_EFT)), '\n')

        for iclass in range(len(list_processClasses)):
            list_LearningWeights_allClasses.append(list_weights_allClasses_abs[iclass])

    elif parameterizedNN == False: #Other non-parameterized strategies (e.g. classifier)

        #Compute 'yields' (from *absolute* weights) to reweight classes
        list_yields_abs_allClasses = []
        yield_abs_total = 0
        for iclass in range(len(list_processClasses)):
            list_yields_abs_allClasses.append(list_weights_allClasses_abs[iclass].sum())
            yield_abs_total+= list_weights_allClasses_abs[iclass].sum()

        #Compute scale factors to rescale each class to 'yield_abs_total'
        list_SFs_allClasses = []
        for iclass in range(len(list_processClasses)):

            list_SFs_allClasses.append(100. / list_yields_abs_allClasses[iclass]) #Compute SF for each process so that its total yield equals N (arbitrary)
            # list_SFs_allClasses.append(yield_abs_total / list_yields_abs_allClasses[iclass])

            print(colors.ital, '* Class', list_labels[iclass], '(' + str(len(list_weights_allClasses[iclass])) + ' entries):', colors.reset)
            print('Default yield = ', float('%.4g' % list_yields_abs_allClasses[iclass]), '[absolute weights]')
            print('Rescaling factor = ', float('%.3g' % list_SFs_allClasses[iclass]))
            print('==> Rescaled yield :', float('%.2g' % (list_yields_abs_allClasses[iclass]*list_SFs_allClasses[iclass])), '\n')

        for iclass in range(len(list_processClasses)):
            list_LearningWeights_allClasses.append(list_weights_allClasses_abs[iclass]*list_SFs_allClasses[iclass]) #Training weights = phys weights * norm_SF
            # print(list_weights_allClasses_abs[iclass][:10]) #Before
            # print(list_LearningWeights_allClasses[iclass][:10]) #After
            # print(list_LearningWeights_allClasses[iclass].sum()) #Recheck sum

    else: #Param. NN <-> unweighted samples <-> training weights = 1

        for iclass in range(len(list_processClasses)):
            list_LearningWeights_allClasses.append(np.ones(len(list_weights_allClasses_abs[iclass])))

    #-- Can artificially manipulate class weights here
    # list_LearningWeights_allClasses[0]*= 0.

    LearningWeights_allClasses = np.concatenate(list_LearningWeights_allClasses, 0)

    if opts["strategy"] is "CARL_singlePoint":
        LearningWeights_allClasses[targetClass_allClasses==0] = LearningWeights_allClasses[targetClass_allClasses==0] * SF_SM
        LearningWeights_allClasses[targetClass_allClasses==1] = LearningWeights_allClasses[targetClass_allClasses==1] * SF_EFT

    return LearningWeights_allClasses, PhysicalWeights_allClasses

# //--------------------------------------------
# //--------------------------------------------

 #    # #####  #####    ##   ##### ######     ####  #####  ##### #  ####  #    #  ####
 #    # #    # #    #  #  #    #   #         #    # #    #   #   # #    # ##   # #
 #    # #    # #    # #    #   #   #####     #    # #    #   #   # #    # # #  #  ####
 #    # #####  #    # ######   #   #         #    # #####    #   # #    # #  # #      #
 #    # #      #    # #    #   #   #         #    # #        #   # #    # #   ## #    #
  ####  #      #####  #    #   #   ######     ####  #        #   #  ####  #    #  ####

def Update_Lists(opts, list_labels, list_features):
    '''
    Update the lists of features and class labels depending on the chosen strategy. For parameterized NN, need to include the theory parameters as inputs. And for parameterized NN, relevant classes from there on are SM/EFT of SM/operators, not physics processes anymore.
    '''

    if opts["parameterizedNN"] == True:

        if opts["testToy1D"]: list_features = ['x']

        list_features = np.append(list_features, opts["listOperatorsParam"]) #Treat the theory parameters theta (WC values of each operator) as input features

        refName = "SM"
        if opts["refPoint"] is not "SM": refName = "refPoint"

        if opts["strategy"] in ["CARL", "ROLR", "RASCAL"]: list_labels = []; list_labels.append("EFT"); list_labels.append(refName) #1=EFT, 0=SM
        # if opts["strategy"] in ["CARL", "ROLR", "RASCAL"]: list_labels = []; list_labels.append(refName); list_labels.append("EFT") #SM vs EFT
        elif opts["strategy"] is "CARL_multiclass": list_labels = opts["listOperatorsParam"][:]; list_labels.insert(0, refName) #SM vs op1 vs op2 vs ... #NB: must specify '[:]' to create a copy, not a reference

    elif opts["strategy"] is "CARL_singlePoint":
        list_labels = []; list_labels.append("EFT"); list_labels.append("SM") #1=EFT, 0=SM

    elif opts["strategy"] is "regressor":
        if all([idx >= 0 for idx in opts["targetVarIdx"]]): #Regressor: use this/these variable(s) as target(s) --> remove from training features
            # list_labels = []; [list_labels.append(list_features[v]) for v in opts["targetVarIdx"]] #Make the 'labels' represent the output variables rather than different processes... needed to get correspondance with each output node
            for idx in sorted(opts["targetVarIdx"], reverse=True): #need to delete in reverse order so that you don't throw off the subsequent indexes.
                del list_features[idx] #Remove feature(s) used as target(s) from list_features (--> don't use for training)
                if opts["comparVarIdx"] > idx: opts["comparVarIdx"] = opts["comparVarIdx"] -1; #Removed 1 feature from list --> Need to update other indices accordingly

    return list_labels, list_features

# //--------------------------------------------
# //--------------------------------------------

 #####   ##   #####   ####  ###### #####  ####
   #    #  #  #    # #    # #        #   #
   #   #    # #    # #      #####    #    ####
   #   ###### #####  #  ### #        #        #
   #   #    # #   #  #    # #        #   #    #
   #   #    # #    #  ####  ######   #    ####

def Get_Targets(opts, list_features, list_processClasses, list_nentries_class, targetClass_allClasses, jointLR_allClasses, scores_allClasses_eachOperator, x):
    '''
    Create and return array 'y' <-> target for classification/regression.
    Also create and return array 'y_process' <-> will keep track of which process each event belongs to (since for regression, target will differ from 0,1).

    NB: exception: for SM vs EFT-only classifiers, SM is treated as the 'signal' (since the corresponding sample is usually defined first) insteas as "background" as in other SM vs EFT strategies.
    '''

    #-- Binary or multiclass classification between different processes
    if opts["strategy"] is "classifier": #Separate SM processes, or SM/pure-EFT --> Target corresponds to process class itself

        #Create array of labels (1 row per event, 1 column per class)
        if opts["nofOutputNodes"] == 1: #binary, single column => sig 1, bkg 0
            y = np.ones(list_nentries_class[0]) #'1' = signal
            if len(list_processClasses)>1:
                y_integer_bkg = np.zeros(list_nentries_class[1]) #'0' = bkg
                y = np.concatenate((y, y_integer_bkg), axis=0)

        else: #multiclass, 1 column per class (y=1 for corresponding process, 0 for others)
            list_y_integer_allClasses = []

            #For each process class, create a 1D array with shape (n_entries) filled with integer associated with the class. Repeat for each class, concatenate, and one-hot encode
            for iclass in range(len(list_nentries_class)): #Concatenate subsequent classes
                list_y_integer_allClasses.append(np.full((list_nentries_class[iclass]), iclass) )
            y_integer = np.concatenate(list_y_integer_allClasses, 0)
            y = utils.to_categorical(y_integer, num_classes=opts["nofOutputNodes"])

        y_process = y #Classification <-> target corresponds to process ID

    #-- Regress some quantity (only 0,1 supported for now, for debug) #Only 2 process classes supported here yet
    elif opts["strategy"] is "regressor":

        mode = 4 #0 <-> regress label 0,1 ; 1 <-> regress on 2 gaussians ; 2 <-> regress dummy value ; 3 <-> regress some input feature; 4 <-> regress on first variable in input list (and remove it for training)

        if all([idx>=0 for idx in opts["targetVarIdx"]]): mode = 4 #Regress on target variable(s) selected by user

        if opts["nofOutputNodes"] != 1 and mode < 4: print('ERROR ! Not supported yet... check Get_Targets()'); exit(1)

        print('\n')
        if mode == 0:

            print(colors.ital, 'Regressor mode = 0: Regress on 0 (bkg) and 1 (sig)...', colors.reset)

            #Regress 0, 1
            y = np.ones(list_nentries_class[0]) #'1' = signal
            if len(list_processClasses) > 1:
                y_integer_bkg = np.zeros(list_nentries_class[1]) #'0' = bkg
                y = np.concatenate((y, y_integer_bkg), axis=0)
            y_process = y #Classification <-> target corresponds to process ID

        elif mode == 1:

            print(colors.ital, 'Regressor mode = 1: Regress on 2 dummy gaussians for sig/bkg...', colors.reset)

            #Targets randomly sampled around 0,1 -- testing
            y = np.random.normal(loc=0.2, scale=0.1, size=list_nentries_class[0])
            y_process = np.ones(list_nentries_class[0]) #'1' = signal
            if len(list_processClasses) > 1:
                y_integer_bkg = np.random.normal(loc=5, scale=0.1, size=list_nentries_class[1])
                y = np.concatenate((y, y_integer_bkg), axis=0)
                y_process = np.concatenate((y_process, np.zeros(list_nentries_class[1])), axis=0) #'0' = bkg

        elif mode == 2:

            print(colors.ital, 'Regressor mode = 2: Regress on fixed dummy value...', colors.reset)

            #Regress fixed dummy value
            dummy = 0.5
            y = np.full(list_nentries_class[0], dummy)
            y_process = np.ones(list_nentries_class[0]) #'1' = signal
            if len(list_processClasses) > 1:
                y_integer_bkg = np.full(list_nentries_class[1], dummy)
                y = np.concatenate((y, y_integer_bkg), axis=0)
                y_process = np.concatenate((y_process, np.zeros(list_nentries_class[1])), axis=0) #'0' = bkg

        elif mode == 3:

            print(colors.ital, 'Regressor mode = 3: Regress on first input feature (used in training ! for debugging only)...', colors.reset)

            #Regress an input feature
            y = x[:list_nentries_class[0], 0] #First feature
            y_process = np.ones(list_nentries_class[0]) #'1' = signal
            if len(list_processClasses) > 1:
                y_integer_bkg = x[list_nentries_class[0]:list_nentries_class[1]+list_nentries_class[1], 0]
                y = np.concatenate((y, y_integer_bkg), axis=0)
                y_process = np.concatenate((y_process, np.zeros(list_nentries_class[1])), axis=0) #'0' = bkg

        elif mode == 4:

            print(colors.ital, 'Regressor mode = 4: Regress on user-selected variable', [list_features[i] for i in opts["targetVarIdx"]], '(not used in training)...', colors.reset)

            #Regress on first input feature, and remove it for training
            y = x[:list_nentries_class[0], opts["targetVarIdx"]] # Use selected variable(s) as target(s)
            y_process = np.ones(list_nentries_class[0]) #'1' = signal
            if len(list_processClasses) > 1:
                y_integer_bkg = x[list_nentries_class[0]:list_nentries_class[1]+list_nentries_class[1], opts["targetVarIdx"]]
                y = np.concatenate((y, y_integer_bkg), axis=0)
                y_process = np.concatenate((y_process, np.zeros(list_nentries_class[1])), axis=0) #'0' = bkg
            x = np.delete(x, opts["targetVarIdx"], 1)  # Delete selected variable(s)=column(s) from x

    #-- For NNs separating SM from EFT, already defined target in dedicated function (not based on 'process class' like for regular classification)
    else:

        if "CARL" in opts["strategy"]: #For CARL, binary label 0/1 indicates whether event was generated at reference point (SM) or not. For CARL_multiclass, there as (1+n_operators) labels: 0=refpoint, 1=operator1 activated, 2=operator2 activated, etc.; i.e. only 1 EFT operator can be non-zero at once for the CARL_multiclass approach.
        # if opts["strategy"] in ["CARL", "CARL_multiclass"]: #For CARL, binary label 0/1 indicates whether event was generated at reference point (SM) or not. For CARL_multiclass, there as (1+n_operators) labels: 0=refpoint, 1=operator1 activated, 2=operator2 activated, etc.; i.e. only 1 EFT operator can be non-zero at once for the CARL_multiclass approach.
            y = targetClass_allClasses; y_process = y #Info already stored when defining EFT points to train on

        elif opts["strategy"] == "ROLR": #Target is joint likelihood ratio r
            y = jointLR_allClasses; y_process = targetClass_allClasses

        elif opts["strategy"] == "RASCAL": #Targets are joint likelihood ratio r and score t (1 component per operator)
            y = np.column_stack((jointLR_allClasses, scores_allClasses_eachOperator)); y_process = targetClass_allClasses

    return y, y_process, x

# //--------------------------------------------
# //--------------------------------------------

def Sanitize_Data(opts, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, singleThetaName=""):
    '''
    Sanitize the data provided as input to the NN.
    '''

    #-- Sanity check (NaN, infinite)
    removeEvent = False
    lists_names= ['x', 'y', 'y_process', 'PhysicalWeights_allClasses', 'LearningWeights_allClasses']
    for il, l in enumerate([x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses]):
        if np.isnan(l).any() or not np.isfinite(l).all():
            print('\n', colors.fg.red, 'WARNING : found a NaN/inf value in array', lists_names[il],'(returned by Get_Data). Removing this event from all arrays...', colors.reset)
            removeEvent = True

    if removeEvent is True:
        #-- Define masks to remove any event (row) containing a NaN/inf value. 'any(axis=1)' <-> look for any row containing a NaN. 'reshape' <-> convert 1D arrays to 2D arrays for convenience
        mask_nan = ~np.isnan(x.reshape(len(x), -1)).any(axis=1) & ~np.isnan(y.reshape(len(y), -1)).any(axis=1) & ~np.isnan(y_process.reshape(len(y_process), -1)).any(axis=1) & ~np.isnan(PhysicalWeights_allClasses.reshape(len(PhysicalWeights_allClasses), -1)).any(axis=1) & ~np.isnan(LearningWeights_allClasses.reshape(len(LearningWeights_allClasses), -1)).any(axis=1)
        mask_inf = ~np.isinf(x.reshape(len(x), -1)).any(axis=1) & ~np.isinf(y.reshape(len(y), -1)).any(axis=1) & ~np.isinf(y_process.reshape(len(y_process), -1)).any(axis=1) & ~np.isinf(PhysicalWeights_allClasses.reshape(len(PhysicalWeights_allClasses), -1)).any(axis=1) & ~np.isinf(LearningWeights_allClasses.reshape(len(LearningWeights_allClasses), -1)).any(axis=1)
        mask = mask_nan & mask_inf #Combine both masks
        x = x[mask]; y = y[mask]; y_process = y_process[mask]; PhysicalWeights_allClasses = PhysicalWeights_allClasses[mask]; LearningWeights_allClasses = LearningWeights_allClasses[mask]
        # print(len(x))

    if opts["strategy"] in ["ROLR", "RASCAL"] and singleThetaName is "":
        if np.all(y_process==0):
            print(colors.fg.orange, "WARNING : all the class labels are set to 0. I notice that you're doing regression, so that may not be an issue. Still, for automated validation plots to work, I will set half the labels to 1!", colors.reset)
            max_idx = int(len(y_process)/2) #Half events
            y_process[::2] = 1 #x[start:stop:step] syntax -> change label of 1 every 2 elements
        elif np.all(y_process==1):
            print(colors.fg.orange, "WARNING : all the class labels are set to 1. I notice that you're doing regression, so that may not be an issue. Still, for automated validation plots to work, I will set half the labels to 0!", colors.reset)
            max_idx = int(len(y_process)/2) #x[start:stop:step] syntax -> change label of 1 every 2 elements
            y_process[::2] = 0

    #-- Check for presence of very large weights
    remove_largeWeight_events = False #NB: not correct because can only remove events from training sample #NB: not correct because the mean should be computed separately for each process class... ! #Remove events with too large weights (may bias the training) -- not for training against pure-EFT samples, where weights are very imbalanced by construction
    if opts["samplesType"] is not "centralVSpureEFT" and len(PhysicalWeights_allClasses[PhysicalWeights_allClasses > np.mean(PhysicalWeights_allClasses)*100]) > 0:
        print('Warning: very large event weights found (global mean = ',np.mean(PhysicalWeights_allClasses),') :')
        print(PhysicalWeights_allClasses[PhysicalWeights_allClasses > np.mean(PhysicalWeights_allClasses)*100])

        if remove_largeWeight_events:
            print('---> Removing these events from all arrays...')
            #-- Define masks to remove any event (row) containing a NaN/inf value. 'any(axis=1)' <-> look for any row containing a NaN. 'reshape' <-> convert 1D arrays to 2D arrays for convenience
            mask = np.where(PhysicalWeights_allClasses <= np.mean(PhysicalWeights_allClasses)*100)
            x = x[mask]; y = y[mask]; y_process = y_process[mask]; PhysicalWeights_allClasses = PhysicalWeights_allClasses[mask]; LearningWeights_allClasses = LearningWeights_allClasses[mask]
            # print(len(x))

    return x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses

# //--------------------------------------------
# //--------------------------------------------

  ####  #####  #      # #####
 #      #    # #      #   #
  ####  #    # #      #   #
      # #####  #      #   #
 #    # #      #      #   #
  ####  #      ###### #   #

#Split into train/test datasets
def Train_Test_Split(opts, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses):
    '''
    Split data into training/testing datasets.

    #http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
    #Default args : shuffle=True <-> shuffle events ; random_state=None <-> random seed ; could also use stratify=y so that the final splitting respects the class proportions of the array y, if desired (else: random)
    '''

    if opts["nEventsTot_train"] != -1 and opts["nEventsTot_test"] != -1: #Specify nof train/test events
        _trainsize=opts["nEventsTot_train"]; _testsize=opts["nEventsTot_test"]
    else: #Specify train/test relative proportions
        _trainsize=opts["splitTrainValTestData"][0]; _testsize=1-_trainsize

    if opts["makeValPlotsOnly"] is True: _trainsize = 0.10 #If not training a NN, use ~ all data for validation ('training data' is meaningless in that case)

    return train_test_split(x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, train_size=_trainsize, test_size=1-_trainsize, shuffle=True)

#Split into train/val/test datasets
#NB: don't simply call sklearn.train_test_split twice, because this makes the fractions complicated (x% of y%...)
def Train_Val_Test_Split(opts, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses):
    '''
    Split data into training/validation/testing datasets.
    '''

    if opts["nEventsTot_train"] != -1 and opts["nEventsTot_val"] != -1 and opts["nEventsTot_test"] != -1: #Specify nof train/test events
        _trainsize=opts["nEventsTot_train"]; _valsize=opts["nEventsTot_val"] ; _testsize=opts["nEventsTot_test"]
    else: #Specify train/val/test relative proportions
        _trainsize=opts["splitTrainValTestData"][0]; _valsize=opts["splitTrainValTestData"][1]; _testsize=opts["splitTrainValTestData"][2]

    if opts["makeValPlotsOnly"] is True: _trainsize = 0.01; _valsize = 0.01; _testsize = 1 - (_trainsize+_valsize) #If not training a NN, use ~ all data for testing ('training/val data' is meaningless in that case)

    x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses = unison_shuffled_copies(x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses) #Shuffle all arrays
    total_nentries = len(x) #Determine total nof entries in arrays

    if _trainsize <= 1: _trainsize = int(_trainsize*total_nentries) #If _trainsize is expressed as a fraction, translate it to a nof events
    if _valsize <= 1: _valsize = int(_valsize*total_nentries) #If _valsize is expressed as a fraction, translate it to a nof events
    if _testsize <= 1: _testsize = int(_testsize*total_nentries) #If _testsize is expressed as a fraction, translate it to a nof events

    #Use train/val/test nof events to split the full arrays #NB: fill validation arrays last, because user may choose not to use validation data
    x_train = x[:_trainsize]; y_train = y[:_trainsize]; y_process_train = y_process[:_trainsize];  PhysicalWeights_train = PhysicalWeights_allClasses[:_trainsize]; LearningWeights_train = LearningWeights_allClasses[:_trainsize]
    x_test = x[_trainsize:(_trainsize+_testsize)]; y_test = y[_trainsize:(_trainsize+_testsize)]; y_process_test = y_process[_trainsize:(_trainsize+_testsize)];  PhysicalWeights_test = PhysicalWeights_allClasses[_trainsize:(_trainsize+_testsize)]; LearningWeights_test = LearningWeights_allClasses[_trainsize:(_trainsize+_testsize)]
    x_val = x[(_trainsize+_testsize):(_trainsize+_testsize+_valsize)]; y_val = y[(_trainsize+_testsize):(_trainsize+_testsize+_valsize)]; y_process_val = y_process[(_trainsize+_testsize):(_trainsize+_testsize+_valsize)];  PhysicalWeights_val = PhysicalWeights_allClasses[(_trainsize+_testsize):(_trainsize+_testsize+_valsize)]; LearningWeights_val = LearningWeights_allClasses[(_trainsize+_testsize):(_trainsize+_testsize+_valsize)]

    return x_train, x_val, x_test, y_train, y_val, y_test, y_process_train, y_process_val, y_process_test, PhysicalWeights_train, PhysicalWeights_val, PhysicalWeights_test, LearningWeights_train, LearningWeights_val, LearningWeights_test

# //--------------------------------------------
# //--------------------------------------------

 ##### #####    ##   #    #  ####  ######  ####  #####  #    #
   #   #    #  #  #  ##   # #      #      #    # #    # ##  ##
   #   #    # #    # # #  #  ####  #####  #    # #    # # ## #
   #   #####  ###### #  # #      # #      #    # #####  #    #
   #   #   #  #    # #   ## #    # #      #    # #   #  #    #
   #   #    # #    # #    #  ####  #       ####  #    # #    #

 # #    # #####  #    # #####  ####
 # ##   # #    # #    #   #   #
 # # #  # #    # #    #   #    ####
 # #  # # #####  #    #   #        #
 # #   ## #      #    #   #   #    #
 # #    # #       ####    #    ####

def Transform_Inputs(weightDir, x_train, list_features, lumiName, parameterizedNN, transfType='quantile'):
    '''
    Get normalization parameters from training data. Give these parameters to NN input layers to directly normalize all inputs.

    transfType: 'quantile', 'range', 'gauss'. Defines the transformation applied to normalize input data.
    '''

    nmax = 100000 #None <-> use all events; else: don't compute means shifts and scales on more than nmax events (slow)
    np.set_printoptions(precision=3)
    xTrainRescaled = None

    # print('Before transformation :', x_train[0:5,:])

    if transfType not in ['quantile', 'range', 'gauss']:
        print('\n', colors.fg.red, 'Warning : transfType', colors.reset, transfType, colors.fg.red, 'not known ! Using default [quantile]', colors.reset)
        transfType = 'quantile'

    #--- QUANTILE RESCALING
    # a = median ; b = scale
    if transfType == 'quantile': #frac = Fraction of event within [-1;+1] -- e.g. 0.68 or 0.95
        # frac=0.99
        frac=0.95
        shift_, scale_ = get_normalization_iqr(x_train[:nmax], frac)
        # if parameterizedNN==False: xTrainRescaled = normalize(x_train, shift_, scale_) #Apply transformation on train events -- for control
        xTrainRescaled = normalize(x_train, shift_, scale_)

    #--- RANGE SCALING
    # a = min ; b = scale
    elif transfType == 'range':
        scaler = MinMaxScaler(feature_range=(-0.5, 0.5)).fit(x_train[:nmax]) #Compute macro parameters
        shift_ = scaler.min_
        scale_ = scaler.scale_
        # if parameterizedNN==False: xTrainRescaled = scaler.transform(x_train) #Apply transformation on train events -- for control
        xTrainRescaled = scaler.transform(x_train)

    #--- RESCALE TO UNIT GAUSSIAN -- https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
    # a = mean ; b = stddev
    elif transfType == 'gauss':
        scaler = StandardScaler().fit(x_train[:nmax]) #Get params
        shift_ = scaler.mean_
        scale_ = scaler.scale_
        if parameterizedNN==False: xTrainRescaled = scaler.transform(x_train) #Apply transformation on train events -- for control
        xTrainRescaled = scaler.transform(x_train)

    text_file = open(weightDir + "NN_info.txt", "a+") #'w' to overwrite

    #Dump shift_ and scale_ params into txtfile
    for ivar in range(len(list_features)):
        # print('Variable', list_features[ivar])
        text_file.write(list_features[ivar]); text_file.write(' ')
        text_file.write(str(shift_[ivar])); text_file.write(' ')
        text_file.write(str(scale_[ivar])); text_file.write('\n')

    text_file.close()
    # print(colors.fg.lightgrey, '\n===> Saved NN infos (input/output nodes names, rescaling values, etc.) in : ', weightDir + "NN_infos.txt \n", colors.reset)

    # print('shift_', shift_); print('scale_', scale_)
    # print('After transformation :', xTrainRescaled[0:5,:])

    return xTrainRescaled, shift_, scale_

# //--------------------------------------------
# //--------------------------------------------
