# Read ROOT files, shape/transform the data (x : input features, y : labels), compute event reweights, ...
#NB : if want to add a validation test, could simply split the test sample again with train_test_split...

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

#Call sub-function to read/store/shape the data
def Get_Data(opts, list_lumiYears, list_processClasses, list_labels, list_features, weightDir, ntuplesDir, lumiName):

    #-- Get data from TFiles
    list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses = Read_Data(list_lumiYears, ntuplesDir, list_processClasses, list_labels, list_features, opts["cuts"])

    #-- For private EFT samples, get the per-event fit coefficients (for later extrapolation at any EFT point) #Central samples: empty arrays
    list_EFT_FitCoeffs_allClasses, list_indexSM_allClasses = Get_EFT_FitCoefficients(list_processClasses, list_labels, list_EFTweights_allClasses, list_EFTweightIDs_allClasses)

    #-- If the DNN is parameterized on Wilson coeffs., need to artificially extend the dataset to train on many different points in EFT phase space
    list_x_allClasses, list_weights_allClasses, list_thetas_allClasses, list_targetClass_allClasses = Extend_Dataset(opts, list_labels, list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses, list_indexSM_allClasses)

    #-- Concatenate and reshape arrays
    x, list_weights_allClasses, thetas_allClasses, targetClass_allClasses, list_nentries_class = Shape_Data(opts, list_x_allClasses, list_weights_allClasses, list_thetas_allClasses, list_targetClass_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses)

    #-- Define the targets 'y' (according to which the classification/regression is performed). Also keep track of the process class indices of all events ('y_process')
    y, y_process = Get_Targets(opts, list_processClasses, list_nentries_class, targetClass_allClasses)

    #-- Define 'physical event weights' (for plotting, ...) and 'training weights' (rescaled arbitrarily to improve the training procedure)
    LearningWeights_allClasses, PhysicalWeights_allClasses = Get_Events_Weights(list_processClasses, list_labels, list_weights_allClasses, opts["parameterizedDNN"])

    #-- Before we randomize the events, store the input values of the very first events (which belong to first process) --> Can use these known events for later validation/comparison
    x_control_firstNEvents = x[0:10,:]

    #-- Define training & testing subsamples (takes care of splitting + shuffling)
    #http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
    #Default args : shuffle=True <-> shuffle events ; random_state=None <-> random seed ; could also use stratify=y so that the final splitting respects the class proportions of the array y, if desired (else : random)
    #NB : to also split into val dataset, could call the function several times with same rnd seed (1:train, 2:test, 3:val)
    if (opts["nEventsTot_train"] is not -1) and (opts["nEventsTot_test"] is not -1): #Specify nof train/test events
        _trainsize=opts["nEventsTot_train"]; _testsize=opts["nEventsTot_test"]
    else: #Specify train/test relative proportions
        _trainsize=opts["splitTrainEventFrac"]; _testsize=1-opts["splitTrainEventFrac"]

    x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test = train_test_split(x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, train_size=_trainsize, test_size=_testsize, shuffle=True)
    # x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, EFTweights_train, EFTweights_test, EFTweightIDs_train, EFTweightIDs_test, EFT_FitCoeffs_train, EFT_FitCoeffs_test = train_test_split(x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, EFTweights_allClasses, EFTweightIDs_allClasses, EFT_FitCoeffs_allClasses, train_size=_trainsize, test_size=_testsize, shuffle=True)

    #-- Get rescaling parameters for each input feature, given to first DNN layer to normalize features -- derived from train data alone
    xTrainRescaled, shifts, scales = Transform_Inputs(weightDir, x_train, list_features, lumiName, opts["parameterizedDNN"])

    print(colors.fg.lightblue, "===========", colors.reset)
    print(colors.fg.lightblue, "-- Will use " + str(x_train.shape[0]) + " training events !", colors.reset)
    print(colors.fg.lightblue, "-- Will use " + str(x_test.shape[0]) + " testing events !", colors.reset)
    print(colors.fg.lightblue, "===========\n", colors.reset)

    return x_train, x_test, y_train, y_test, y_process_train, y_process_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, x, y, y_process, PhysicalWeights_allClasses, LearningWeights_allClasses, shifts, scales, x_control_firstNEvents, xTrainRescaled
# //--------------------------------------------
# //--------------------------------------------









# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
########  ########    ###    ########     ########     ###    ########    ###
##     ## ##         ## ##   ##     ##    ##     ##   ## ##      ##      ## ##
##     ## ##        ##   ##  ##     ##    ##     ##  ##   ##     ##     ##   ##
########  ######   ##     ## ##     ##    ##     ## ##     ##    ##    ##     ##
##   ##   ##       ######### ##     ##    ##     ## #########    ##    #########
##    ##  ##       ##     ## ##     ##    ##     ## ##     ##    ##    ##     ##
##     ## ######## ##     ## ########     ########  ##     ##    ##    ##     ##
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Read the data from ROOT files and store it in np arrays
def Read_Data(list_lumiYears, ntuplesDir, list_processClasses, list_labels, list_features, cuts):

    testdirpath = Path(ntuplesDir)
    if not testdirpath.is_dir():
        print('Ntuple dir. '+ntuplesDir+' not found ! Abort !')
        exit(1)

    list_x_allClasses = [] #List of x-arrays storing the values of input features (for all events, all considered years, all physics processes in a given class) -- 1 array per process class
    list_weights_allClasses = [] #Idem for event weights
    list_EFTweights_allClasses = [] #Idem for EFT reweights
    list_EFTweightIDs_allClasses = [] #Idem for EFT reweights IDs
    for procClass, label in zip(list_processClasses, list_labels): #Loop on classes of physics processes (e.g. 'Signal','Backgrounds') #NB: can not use predefined keyword 'class'

        print('\n', colors.fg.purple, colors.underline, '* Class :', label, colors.reset)

        list_x_proc = [] #List of x-arrays storing the values of input features for all events, for all considered years -- 1 array per physics process (sample)
        list_weights_proc = [] #Idem for central event weights
        list_EFTweights_proc = []; list_EFTweightIDs_proc = [] #Idem for EFT reweights
        for process in procClass: #Loop on physics processes (samples)

            isPrivMCsample = False
            if ("PrivMC" in process and "PrivMC" not in label) or ("PrivMC" in label and "PrivMC" not in process): #Avoid ambiguities
                print('\n', colors.fg.red, 'Error : keyword \'PrivMC\' must be present both in process and class names, or not at all', colors.reset)

            elif "PrivMC" in process and "PrivMC" in label: #Check whether EFT reweights should be looked for
                isPrivMCsample = True
                # print('process', process, ', label', label, ' --> isPrivMCsample')

            print('\n', colors.fg.pink, '* Process :', colors.reset, process)

            for iyear in range(len(list_lumiYears)): #Concatenate data (input features, event weights) for all considered years

                filepath = ntuplesDir + list_lumiYears[iyear] + '/' + process + '.root'
                if not Path(filepath).is_file():
                    print('File '+filepath+' not found ! Abort !')
                    exit(1)

                # print(colors.fg.lightgrey, '* Opening file:', colors.reset, ' ', filepath)
                file = TFile.Open(filepath)
                tree = file.Get('result')
                print(colors.fg.lightgrey, '* Opened file:', colors.reset, filepath, '(', tree2array(tree, branches="eventWeight", selection=cuts).shape[0], 'entries )') #Dummy variable, just to read the nof entries

                list_x_proc.append(tree2array(tree, branches=list_features, selection=cuts)) #Store values of input features into array, append to list

                #-- Store event weights into array, append to list
                list_weights_proc.append(tree2array(tree, branches="eventWeight*eventMCFactor", selection=cuts))
                # list_weights_proc.append(tree2array(tree, branches="eventWeight", selection=cuts))

            if isPrivMCsample: #For private MC samples, get the EFT reweights (properly normalized) and their IDs
                list_EFTweights_proc, list_EFTweightIDs_proc = Read_Data_EFT_File(list_lumiYears, ntuplesDir, process, cuts)

        #Concatenate the different arrays (for all years, processes) corresponding to a single class of process, and append them to their lists --> 1 array per process class
        list_x_allClasses.append(np.concatenate(list_x_proc))
        # list_weights_allClasses.append(np.concatenate(list_weights_proc))

        if "PrivMC" in label:
            list_weights_allClasses.append(np.concatenate(list_EFTweights_proc)[:,0]) #For private EFT samples, only the weights in 'list_EFTweights_proc' make sense (not those in 'list_weights_proc') ; but need single value per event --> use element corresponding to baseline (first column)
            list_EFTweights_allClasses.append(np.concatenate(list_EFTweights_proc))
            list_EFTweightIDs_allClasses.append(np.concatenate(list_EFTweightIDs_proc))

        else:
            weights_tmp = np.concatenate(list_weights_proc)
            nentries = weights_tmp.shape[0]
            list_weights_allClasses.append(weights_tmp)
            list_EFTweights_allClasses.append(np.full(shape=nentries, fill_value=-1)); list_EFTweightIDs_allClasses.append(np.full(shape=nentries, fill_value='-1', dtype=object)) #Append empty arrays (1 row per event) for proper ordering

    print('\n\n')

    return list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses


# //--------------------------------------------
# //--------------------------------------------


#For private MC (EFT) samples, retrieve the EFT reweights and their IDs. Directly normalize properly the EFT weights
def Read_Data_EFT_File(list_lumiYears, ntuplesDir, process, cuts):

    list_EFTweights_proc = []
    list_EFTweightIDs_proc = []
    for iyear in range(len(list_lumiYears)):

        filepath = ntuplesDir + list_lumiYears[iyear] + '/' + process + '.root'
        if not Path(filepath).is_file():
            print('File '+filepath+' not found ! Abort !')
            exit(1)

        file = TFile.Open(filepath)
        tree = file.Get('result')
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

        #Individually for each year, get the EFT reweights IDs
        array_EFTweightIDs_proc = np.stack(tree2array(tree, branches="mc_EFTweightIDs", selection=cuts)) #stack : array of arrays -> 2d array

        #Individually for each year, get the EFT reweights and directly divide their values by the corresponding SWEs
        array_EFTweights_proc = np.stack(tree2array(tree, branches="mc_EFTweights", selection=cuts)) #stack : array of arrays -> 2d array
        hist = file.Get("EFT_SumWeights") #Sums of weights for each EFT reweight is stored in histogram
        array_EFT_SWE_proc = hist2array(hist)
        array_EFT_SWE_proc = array_EFT_SWE_proc[0:array_EFTweights_proc.shape[1]] #Only need the SWE values for the considered reweight points
        array_EFTweights_proc = np.divide(array_EFTweights_proc, array_EFT_SWE_proc)
        normWeight_proc = tree2array(tree, start=0, stop=1, branches="eventWeight*eventMCFactor/weightMENominal", selection=cuts)[0] #Normalization factor
        array_EFTweights_proc = np.multiply(array_EFTweights_proc, normWeight_proc)

        #Manually find and remove all weights with unproper naming conventions (for example 'rwgt_1' nominal weight is included by default by MG)
        array_EFTweights_proc, array_EFTweightIDs_proc = Remove_Unnecessary_EFTweights(array_EFTweights_proc, array_EFTweightIDs_proc)

    list_EFTweights_proc.append(array_EFTweights_proc) #Append array of EFT reweights (for given year) to list
    list_EFTweightIDs_proc.append(array_EFTweightIDs_proc) #Append array of EFT reweights IDs (for given year) to list

    return list_EFTweights_proc, list_EFTweightIDs_proc










# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
######## ######## ########     ######   #######  ######## ######## ########
##       ##          ##       ##    ## ##     ## ##       ##       ##
##       ##          ##       ##       ##     ## ##       ##       ##
######   ######      ##       ##       ##     ## ######   ######   ######
##       ##          ##       ##       ##     ## ##       ##       ##
##       ##          ##       ##    ## ##     ## ##       ##       ##       ###
######## ##          ##        ######   #######  ######## ##       ##       ###
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Get the 'fit coefficients' (1 per component entering the squared matrix element of the EFT process) A satisfying : A.T=w, with w the benchmark weights and T the matrix of 'effective WCs' corresponding to the benchmark points.
#Once these fit coefficients are extracted for a given event, they can be used to extrapolate the event weight at any new EFT point
def Get_EFT_FitCoefficients(list_processClasses, list_labels, list_EFTweights_allClasses, list_EFTweightIDs_allClasses):

    # for iclass in range(len(list_EFTweights_allClasses)): list_EFTweights_allClasses[iclass]=list_EFTweights_allClasses[iclass][:5]; list_EFTweightIDs_allClasses[iclass]=list_EFTweightIDs_allClasses[iclass][:5] #For debugging

    list_indexSM_allClasses = [] #Store index of SM point for all classes
    list_EFT_FitCoeffs_allClasses = []
    for iclass in range(len(list_processClasses)): #Loop on classes of physics processes

        if "PrivMC" in list_processClasses[iclass][0] and "PrivMC" in list_labels[iclass]: #Check whether EFT reweights should be looked for

            operatorNames, operatorWCs, idx_SM = Parse_EFTpoint_IDs(list_EFTweightIDs_allClasses[iclass][0]) #Get the lists of operator names and WC values for this process #NB: assumes that they are identical for all events in this process
            n_components, components = Find_Components(operatorNames[0]) #Determine the components required to parameterize the event weight #NB: assumes that they are identical for all events in this process
            effWC_components = Get_EffectiveWC_eachComponents(n_components, components, operatorWCs) #Determine the 'effective WC' values associated with each component, for each benchmark point
            fit_coeffs = Get_FitCoefficients(effWC_components, benchmark_weights=list_EFTweights_allClasses[iclass]) #Determine the fit coefficients of the events, based on the benchmark weights and 'effective WC' values

            #-- Debug
            # print('WCs',effWC_components.shape,  effWC_components[:5])
            # print('coeffs', fit_coeffs.shape, fit_coeffs[:5])
            # print('w0', np.dot(effWC_components[0],np.transpose(fit_coeffs[0])))
            # print('w1 ev1', np.dot(effWC_components[1],np.transpose(fit_coeffs[1])))
            # print('w2 ev2', np.dot(effWC_components[2],np.transpose(fit_coeffs[2])))
            # print('w2 ev2', np.dot(np.transpose(fit_coeffs[2]), effWC_components[2]))
            # print('bench', list_EFTweights_allClasses[iclass][:3])
            # exit(1)

            list_EFT_FitCoeffs_allClasses.append(fit_coeffs) #Append fit coeffs to list
            list_indexSM_allClasses.append(idx_SM)

        else: list_EFT_FitCoeffs_allClasses.append(np.full(shape=len(list_EFTweightIDs_allClasses[iclass]), fill_value=0)) #If process is not EFT, fill with dummy value for now (maintain ordering)

    return list_EFT_FitCoeffs_allClasses, list_indexSM_allClasses











# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
 ######  ##     ##    ###    ########  ########    ########     ###    ########    ###
##    ## ##     ##   ## ##   ##     ## ##          ##     ##   ## ##      ##      ## ##
##       ##     ##  ##   ##  ##     ## ##          ##     ##  ##   ##     ##     ##   ##
 ######  ######### ##     ## ########  ######      ##     ## ##     ##    ##    ##     ##
      ## ##     ## ######### ##        ##          ##     ## #########    ##    #########
##    ## ##     ## ##     ## ##        ##          ##     ## ##     ##    ##    ##     ##
 ######  ##     ## ##     ## ##        ########    ########  ##     ##    ##    ##     ##
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Properly shape the arrays and concatenate them (for all years, processes, etc.)
#NB: nominal weights get concatenated in dedicated function
def Shape_Data(opts, list_x_allClasses, list_weights_allClasses, list_thetas_allClasses, list_targetClass_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses):

    #-- root_numpy 'tree2array' function returns numpy structured array : 1D array whose length equals the nof events, and each element is a structure with multiple fields (1 per feature)
    #For manipulation, it is easier to convert structured arrays obtained in this way into regular numpy arrays (e.g. x will be 2D and have shape (n_events, n_features) )
    #NB: EFT weights/IDs arrays already have proper 2D shapes (due to np.stack ?)
    list_x_arrays_allClasses = []
    for iclass in range(len(list_x_allClasses)):
        list_x_arrays_allClasses.append( list_x_allClasses[iclass].view(np.float32).reshape( (len(list_x_allClasses[iclass]), -1) ) ) #np.view <-> different view of same data ; here, enforces proper data type. Reshape is used to 'unroll' 1d elements into 2d ('-1' can be used when the new dimension can be guessed by numpy from the input data)

    #--- Get nof entries for each class
    list_nentries_class = []
    for iclass in range(len(list_x_arrays_allClasses)): list_nentries_class.append(len(list_x_arrays_allClasses[iclass]))

    maxEvents = opts["maxEvents"]

    #--- Max nof events for train/test phases
    if maxEvents is not -1 and opts["parameterizedDNN"] is False:
    # if maxEvents is not -1:

        #Skim each process class (keep only maxEvents events) ; skim all relevant arrays coherently
        for iclass in range(len(list_x_arrays_allClasses)):

            if list_nentries_class[iclass] > maxEvents:

                #-- If only consider part of the process class data, shuffle events, so that the resulting dataset is representative of the event proportions of each process within the class (else, it could happen that e.g. only events from the first process get considered)
                #NB: only arrays potentially used for training (x, weight, theta) must be shuffled simultaneously
                if len(list_EFTweights_allClasses[iclass]) <= 1: arr=[list_x_arrays_allClasses,list_weights_allClasses]
                else: arr=[list_x_arrays_allClasses,list_weights_allClasses,list_thetas_allClasses,list_targetClass_allClasses]
                unison_shuffled_copies(arr)

                list_nentries_class[iclass] = maxEvents
                list_tmp = [list_x_arrays_allClasses, list_weights_allClasses, list_thetas_allClasses,list_targetClass_allClasses]
                for l in list_tmp:
                    if len(l)>0: l[iclass] = l[iclass][0:maxEvents]

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

    #-- Transform list of arrays --> single concatenated array
    x = np.concatenate(list_x_arrays_allClasses, 0)
    # EFTweights_allClasses = np.concatenate(list_EFTweights_allClasses, 0)
    # EFTweightIDs_allClasses = np.concatenate(list_EFTweightIDs_allClasses, 0)
    # EFT_FitCoeffs_allClasses = np.concatenate(list_EFT_FitCoeffs_allClasses, 0)

    #-- Parameterized DNN: pass the values of the WCs as input features
    if opts["parameterizedDNN"]==True and len(list_thetas_allClasses)>0:
        thetas_allClasses = np.concatenate(list_thetas_allClasses, 0)
        targetClass_allClasses = np.concatenate(list_targetClass_allClasses, 0)

        #Theta has as many columns as there are EFT operators generated in the sample (needed for extraction of fit coefficients from benchmark weights). But from there, only want to retain EFT operators which the DNN will get trained on --> Only parameterize DNN on such operators, not the others (not used)
        theta_tmp = thetas_allClasses[:, ~np.all(thetas_allClasses==0, axis=0)] #Only keep columns (operators) which were activated by the user #'~' is negation
        targetClass_allClasses = targetClass_allClasses[:, ~np.all(targetClass_allClasses==0, axis=0)]
        x = np.append(x, theta_tmp, axis=1)
    else:
        thetas_allClasses = np.zeros((1,1)); targetClass_allClasses = thetas_allClasses #Dummy arrays

    # print(x.shape)

    return x, list_weights_allClasses, thetas_allClasses, targetClass_allClasses, list_nentries_class
    # return x, list_weights_allClasses, EFTweights_allClasses, EFTweightIDs_allClasses, EFT_FitCoeffs_allClasses, thetas_allClasses, list_nentries_class









# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
########    ###    ########   ######   ######## ########  ######
   ##      ## ##   ##     ## ##    ##  ##          ##    ##    ##
   ##     ##   ##  ##     ## ##        ##          ##    ##
   ##    ##     ## ########  ##   #### ######      ##     ######
   ##    ######### ##   ##   ##    ##  ##          ##          ##
   ##    ##     ## ##    ##  ##    ##  ##          ##    ##    ##
   ##    ##     ## ##     ##  ######   ########    ##     ######
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Create and return array 'y' <-> target for classification/regression
#Create and return array 'y_process' <-> will keep track of which process each event belongs to (since for regression, target will differ from 0,1)
def Get_Targets(opts, list_processClasses, list_nentries_class, targetClass_allClasses):

#-- CLASSIFICATION
#NB: execute these commands also for regression, in order to get array 'y_process'

    if opts["parameterizedDNN"] == False: #Separate SM processes, or SM/pure-EFT --> Target corresponds to process class itself

        #Create array of labels (1 row per event, 1 column per class)
        if opts["nofOutputNodes"] == 1: #binary, single column => sig 1, bkg 0
            y_integer_sig = np.ones(list_nentries_class[0]) #'1' = signal
            y = y_integer_sig
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

    else: #Separate SM/EFT at any point in EFT phase space --> Target corresponds to EFT operator which is activated in a given event (0 <-> SM)
        y = targetClass_allClasses #Info already stored when defining EFT points to train on

    y_process = y #For classification, target corresponds to process

#-- REGRESSION

    if opts["regress"]==True:

        if opts["nofOutputNodes"] == 1: #Target = 0,1
            y_integer_sig = np.ones(list_nentries_class[0]) #'1' = signal
            y = y_integer_sig

            #Targets randomly sampled around 0,1 -- testing
            # y_integer_sig = np.random.normal(loc=1.0, scale=0.05, size=list_nentries_class[0])
            # y_integer_bkg = np.random.normal(loc=0.0, scale=0.05, size=list_nentries_class[1])

            if len(list_processClasses)>1:
                y_integer_bkg = np.zeros(list_nentries_class[1]) #'0' = bkg
                y = np.concatenate((y, y_integer_bkg), axis=0)

        else: print('ERROR ! Not supported yet...'); exit(1)

    return y, y_process









# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
##      ## ######## ####  ######   ##     ## ########  ######
##  ##  ## ##        ##  ##    ##  ##     ##    ##    ##    ##
##  ##  ## ##        ##  ##        ##     ##    ##    ##
##  ##  ## ######    ##  ##   #### #########    ##     ######
##  ##  ## ##        ##  ##    ##  ##     ##    ##          ##
##  ##  ## ##        ##  ##    ##  ##     ##    ##    ##    ##
 ###  ###  ######## ####  ######   ##     ##    ##     ######
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Compute and apply weights to training dataset to balance the training
def Get_Events_Weights(list_processClasses, list_labels, list_weights_allClasses, parameterizedDNN):

    #Dupplicate list of weight arrays, storing absolute weights (can't handle negative weights in training)
    list_weights_allClasses_abs = []
    for weights_class in list_weights_allClasses:
        list_weights_allClasses_abs.append(np.absolute(weights_class))

    #Compute 'yields' (from *absolute* weights) to reweight classes
    list_yields_abs_allClasses = []
    yield_abs_total = 0
    for i in range(len(list_processClasses)):
        list_yields_abs_allClasses.append(list_weights_allClasses_abs[i].sum())
        yield_abs_total+= list_weights_allClasses_abs[i].sum()

    #Compute scale factors to rescale each class to 'yield_abs_total'
    list_SFs_allClasses = []
    for i in range(len(list_processClasses)):
        if parameterizedDNN == False:
            list_SFs_allClasses.append(100. / list_yields_abs_allClasses[i]) #Compute SF for each process so that its total yield equals N (arbitrary)
            # list_SFs_allClasses.append(yield_abs_total / list_yields_abs_allClasses[i])
        else: list_SFs_allClasses.append(1)

        print('* Class', list_labels[i], ' :')
        print('-- Default yield = ', float('%.4g' % list_yields_abs_allClasses[i]))
        print('-- Rescaling factor = ', float('%.3g' % list_SFs_allClasses[i]))
        print('===> Rescaled yield :', float('%.2g' % (list_yields_abs_allClasses[i]*list_SFs_allClasses[i])), '\n')
        # print('Class', list_labels[i], ' / Yield = ', float('%.4g' % list_yields_abs_allClasses[i]), ' / Scale factor = ', float('%.3g' % list_SFs_allClasses[i]), '===> Rescaled yield :', float('%.2g' % (list_yields_abs_allClasses[i]*list_SFs_allClasses[i])) )
    print('\n')

    #Get array of reweighted 'training' weights, i.e. used for training only and which are not physical
    list_LearningWeights_allClasses = []
    for i in range(len(list_processClasses)):
        list_LearningWeights_allClasses.append(list_weights_allClasses_abs[i]*list_SFs_allClasses[i])
    LearningWeights_allClasses = np.concatenate(list_LearningWeights_allClasses, 0)

    #Also create corresponding array of physical event weights, to get correct plots, etc.
    PhysicalWeights_allClasses = np.concatenate(list_weights_allClasses, 0)

    return LearningWeights_allClasses, PhysicalWeights_allClasses









# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
######## ########     ###    ##    ##  ######  ########  #######  ########  ##     ##    #### ##    ## ########  ##     ## ########  ######
   ##    ##     ##   ## ##   ###   ## ##    ## ##       ##     ## ##     ## ###   ###     ##  ###   ## ##     ## ##     ##    ##    ##    ##
   ##    ##     ##  ##   ##  ####  ## ##       ##       ##     ## ##     ## #### ####     ##  ####  ## ##     ## ##     ##    ##    ##
   ##    ########  ##     ## ## ## ##  ######  ######   ##     ## ########  ## ### ##     ##  ## ## ## ########  ##     ##    ##     ######
   ##    ##   ##   ######### ##  ####       ## ##       ##     ## ##   ##   ##     ##     ##  ##  #### ##        ##     ##    ##          ##
   ##    ##    ##  ##     ## ##   ### ##    ## ##       ##     ## ##    ##  ##     ##     ##  ##   ### ##        ##     ##    ##    ##    ##
   ##    ##     ## ##     ## ##    ##  ######  ##        #######  ##     ## ##     ##    #### ##    ## ##         #######     ##     ######
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#-- Get normalization parameters from training data. Give these parameters to DNN input layers to directly normalize all inputs
def Transform_Inputs(weightDir, x_train, list_features, lumiName, parameterizedDNN, transfType='quantile'):

    nmax = 100000 #Don't compute means shifts and scales on more than nmax events (slow)
    np.set_printoptions(precision=3)
    xTrainRescaled = None

    # print('Before transformation :', x_train[0:5,:])

    if transfType not in ['quantile', 'range', 'gauss']:
        print('\n', colors.fg.red, 'Warning : transfType', colors.reset, transfType, colors.fg.red, 'not known ! Using default [quantile]', colors.reset)
        transfType = 'quantile'

    #--- QUANTILE RESCALING
    # a = median ; b = scale
    if transfType == 'quantile':
        frac=0.95 #Fraction of event within [-1;+1] -- 0.68 or 0.95
        shift_, scale_ = get_normalization_iqr(x_train[:nmax], frac)
        # if parameterizedDNN==False: xTrainRescaled = normalize(x_train, shift_, scale_) #Apply transformation on train events -- for control
        xTrainRescaled = normalize(x_train, shift_, scale_)

    #--- RANGE SCALING
    # a = min ; b = scale
    elif transfType == 'range':
        scaler = MinMaxScaler(feature_range=(-1, 1)).fit(x_train[:nmax]) #Conpute macro parameters
        shift_ = scaler.min_
        scale_ = scaler.scale_
        # if parameterizedDNN==False: xTrainRescaled = scaler.transform(x_train) #Apply transformation on train events -- for control
        xTrainRescaled = scaler.transform(x_train)

    #--- RESCALE TO UNIT GAUSSIAN -- https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
    # a = mean ; b = stddev
    elif transfType == 'gauss':
        scaler = StandardScaler().fit(x_train[:nmax]) #Get params
        shift_ = scaler.mean_
        scale_ = scaler.scale_ # = np.sqrt(var_)
        if parameterizedDNN==False: xTrainRescaled = scaler.transform(x_train) #Apply transformation on train events -- for control
        xTrainRescaled = scaler.transform(x_train)

    text_file = open(weightDir + "DNN_infos.txt", "w")

    #Dump shift_ and scale_ params into txtfile
    for ivar in range(len(list_features)):
        text_file.write(list_features[ivar]); text_file.write(' ')
        text_file.write(str(shift_[ivar])); text_file.write(' ')
        text_file.write(str(scale_[ivar])); text_file.write('\n')

    text_file.close()
    print(colors.fg.lightgrey, '\n===> Saved DNN infos (input/output nodes names, rescaling values, etc.) in : ', weightDir + "DNN_infos.txt \n", colors.reset)

    # print('shift_', shift_); print('scale_', scale_)
    # print('After transformation :', xTrainRescaled[0:5,:])

    return xTrainRescaled, shift_, scale_
