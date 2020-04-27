# Author : Nicolas Tonon, DESY
# Helper functions related to the EFT treatment in the DNN (parameterization of weight .vs. WCs, extrapolation, computation of score/likelihood, etc.

DEBUG_ = 0 #0: no debug printouts; 1: minimal printouts; 2: maximal printouts (arrays, etc.)

# //--------------------------------------------

'''

#TODO :
-

# NOTES :
- Only applies to private EFT samples
- Invert matrices for all events at once ?
- Interpret power=0 <-> SM and power < max <-> interference
- properly restrict size of arrays to minimum
- how to deal with rwgt_1, rwgt_sm ?

# LIMITATIONS / ASSUMPTIONS :
- Can not sum several processes into a process class including an EFT process. Each EFT process must constitute a separate class.
- Assume that a given EFT sample has the exact same list of reweight points (in same order) for all considered years (could change that in the future by estimating fit coefficients directly when reading the sample, separately for each year)
- If single operator identified, assume sample is pure-EFT only
- Can only rescale array of event if all events have exact same EFT parameterization

# NAMING CONVENTIONS :
- 'WC' = Wilson coefficient <-> the scalar value scaling a given EFT operator in the Lagrangian
- 'Parameter'/'Operator' = EFT operator
- 'Benchmark point' = EFT point (<-> given WC values for given EFT operators) at which the event weight has been evaluated by MG. The ensemble of the benchmark points is also refered to as the 'basis'
- 'Benchmark weight' = corresponding weight value (evaluated by MG with the reweighting module)
- 'Benchmark ID' = string encoding of the name of the benchmark point
- 'Components' = individual terms contributing to the squared matrix element of the process. For example, if considering a single operator, there are 3 components to account for (corresponding to the SM squared term, the interference term, and the pure-EFT term)
- ...

'''

# //--------------------------------------------

# from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
# from root_numpy import root2array, tree2array, array2root
import numpy as np
import pandas as pd
from pathlib import Path
import itertools
from Utils.Helper import CheckName_EFTpoint_ID
from Utils.ColoredPrintout import colors

import warnings
warnings.filterwarnings("ignore", message="in the future insert will treat boolean arrays and array-likes as boolean index instead of casting it to integer")

# //--------------------------------------------
# //--------------------------------------------

 #     #
 #     # ###### #      #####  ###### #####
 #     # #      #      #    # #      #    #
 ####### #####  #      #    # #####  #    #
 #     # #      #      #####  #      #####
 #     # #      #      #      #      #   #
 #     # ###### ###### #      ###### #    #

# //--------------------------------------------
# //--------------------------------------------

def Remove_Unnecessary_EFTweights(array_EFTweights, array_EFTweightIDs):
    """
    Look for EFT weight names which do not follow the expected naming convention (of the kind 'rwgt_ctZ_5.2'), and removes coherently these elements in the arrays of EFT weights/names. Returns modified arrays.
    """

    indices = np.char.find(array_EFTweightIDs[0,:].astype('U'), 'rwgt_c') #Get column indices of reweight names not following the expected naming convention
    indices = np.asarray(indices, dtype=bool) #Convert indices to booleans (True=does not follow naming convention)

    #Remove coherently the incorrect weights
    array_EFTweightIDs = np.delete(array_EFTweightIDs, indices==True, axis=1)
    array_EFTweights = np.delete(array_EFTweights, indices==True, axis=1)

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
    operatorWCs (ndarray of shape [n_points, n_operators]) : array of corresponding WC values, for all (n_points) EFT points
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

    if DEBUG_:
        print('\n-- Parse_EFTpoint_IDs'); print(operatorNames.shape); print(operatorWCs.shape)
        if DEBUG_>1: print(operatorNames); print(operatorWCs)
    return operatorNames, operatorWCs, idx_SM

# //--------------------------------------------
# //--------------------------------------------

  ####   ####  #    # #####   ####  #    # ###### #    # #####  ####
 #    # #    # ##  ## #    # #    # ##   # #      ##   #   #   #
 #      #    # # ## # #    # #    # # #  # #####  # #  #   #    ####
 #      #    # #    # #####  #    # #  # # #      #  # #   #        #
 #    # #    # #    # #      #    # #   ## #      #   ##   #   #    #
  ####   ####  #    # #       ####  #    # ###### #    #   #    ####

#Adapted from madminer code (https://github.com/diana-hep/madminer/blob/080e7fefc481bd6aae16ba094af0fd6bfc301dff/madminer/utils/morphing.py#L116)
#NB : shall interpret row full of 0 <-> pure SM ; any row with sum of powers < 2 <-> interference with SM
def Find_Components(operatorNames):
    """
    Determine the number of components necessary to compute the squared matrix element, and at what power each EFT operator enters each component.
    The components are determined based on the list of considered EFT operators alone, assuming the maximal power at which each operator may contribute to the squared matrix element (2 by default).

    Parameters:
    operatorNames (ndarray of shape [n_operators]) : names of the EFT operators

    Returns:
    n_components (int) : number of components
    components (ndarray of shape [n_components, n_operators]) : array whose elements represent the power at which an operator 'enters' (scales?) a component
    """

    #Want to determine the components based on a single EFT point. If many points provided, assume they all follow the same convention, and only consider the first point
    if operatorNames.ndim > 1:
        print(colors.ital, 'Warning : operatorNames has dimension', operatorNames.ndim, ', suggesting that lists of operator names have been provided for several EFT points. It is expected that all point follow the same naming convention and include the same operators. Will consider only the first row of this array.', colors.reset)
        operatorNames = operatorNames[0, :] #Restrict to first EFT point

    minPower_perOperator = 0
    maxPower_perOperator = 2
    if len(operatorNames) == 1: #Assumption : 1 operator <-> pure EFT
        minPower_perOperator = 2
        maxPower_perOperator = 2

    #Set the min/max powers at which each operator can enter a component
    parameter_max_power = []
    for iter in range(len(operatorNames)):
        parameter_max_power.append(maxPower_perOperator)

    components = []
    for i in range(len(operatorNames)):

        #Get all arrays corresponding to all possible combinations of operator/power
        powers_each_component = [range(max_power+1) for max_power in parameter_max_power] #'+1' because index starts at 0
        for powers in itertools.product(*powers_each_component): #List all possible combinations of range-arrays in 'powers_each_component' (<-> get all possible contributions for each operator satisfying the overall constraints) #itertools.product() computes the cartesian product of input iterables. It is equivalent to nested for-loops.  For example, product(A, B) returns the same as ((x,y) for x in A for y in B).

            powers = np.array(powers, dtype=np.int) #Get corresponding array for given operator

            if np.sum(powers) < minPower_perOperator or np.sum(powers) > maxPower_perOperator: continue #Check the total EFT power

            if not any((powers == x).all() for x in components): #Append combination if not yet included
                components.append(np.copy(powers))

    components = np.array(components, dtype=np.int) #Convert list to array
    n_components = len(components)

    if DEBUG_:
        print('\n-- Find_Components'); print('n_components', n_components); print(components.shape)
        if DEBUG_>1: print(components)
    return n_components, components

# //--------------------------------------------
# //--------------------------------------------

 ###### ###### ######        #    #  ####
 #      #      #             #    # #    #
 #####  #####  #####         #    # #
 #      #      #      ###    # ## # #
 #      #      #      ###    ##  ## #    #
 ###### #      #      ###    #    #  ####

def Get_EffectiveWC_eachComponents(n_components, components, operatorWCs):
    """
    Determine the 'effective WC' (or 'scaling factor') for each component, i.e. the product of the corresponding WC values

    Parameters:
    n_components (int) : number of components
    components (ndarray of shape [n_components, n_operators]) : array whose elements represent the power at which an operator 'enters' (scales?) a component
    operatorWCs (ndarray of shape [n_points, n_operators]) : array of corresponding WC values, for all (n_points) EFT points

    Returns:
    effWC_components (ndarray of shape [n_points, n_components]) : 'effective WC' of given component, for all (n_points) EFT points
    """

    effWC_components = np.zeros((len(operatorWCs), n_components)) #Create empty array of shape (n_points, n_components)

    for y in range(len(operatorWCs)): #For each EFT point
        for x in range(n_components): #For each component
            factor = 1.0
            for p in range(len(np.transpose(operatorWCs))): #For each operator (transpose so that 1st dim corresponds to operators)
                # if DEBUG_>1: print('x', x, 'y', y, 'p', p)
                # print(operatorWCs[p]); print(components[x, p]); print(factor)
                factor *= float(operatorWCs[y, p] ** components[x, p]) #Contribution of this operator to this component
                effWC_components[y, x] = factor

    if DEBUG_:
        print('\n-- Get_EffectiveWC_eachComponents'); print(effWC_components.shape)
        if DEBUG_>1: print(effWC_components)
    return effWC_components

# //--------------------------------------------
# //--------------------------------------------

 ###### # #####     ####   ####  ###### ###### ######
 #      #   #      #    # #    # #      #      #
 #####  #   #      #      #    # #####  #####  #####
 #      #   #      #      #    # #      #      #      ###
 #      #   #      #    # #    # #      #      #      ###
 #      #   #       ####   ####  ###### #      #      ###

def Get_FitCoefficients(effWC_components, benchmark_weights):
    """
    Determine the 'fit coefficients', i.e. the structure constants associated with the components
    Corresponds to matrix A, as in : T . A = w <-> A = w . T^-1 (where w are the benchmark weights, and T the 'effective WCs' for all components and all benchmark points)

    Parameters:
    effWC_components (ndarray of shape [n_points, n_components]) : array of 'effective WC' for all components, for all 'n_points' considered EFT points ; retain only the necessary 'n_components' first benchmark weights
    benchmark_weights (ndarray of shape [n_events, n_points]) : benchmark weights for each event ; retain only the necessary 'n_components' first benchmark weights

    Returns:
    fit_coeffs (ndarray of shape [n_events, n_components]) : fit coefficients associated to the components, for all (n_points) EFT points
    """

    #Only retain minimal required nof benchmark weights / components
    if effWC_components.shape[0] > effWC_components.shape[1]:
        effWC_components = effWC_components[:effWC_components.shape[1], :]
    if benchmark_weights.shape[1] > effWC_components.shape[1]:
        benchmark_weights = benchmark_weights[:, :effWC_components.shape[1]]

    #Need square matrix, i.e. as many benchmark points as there are components
    #  SVD ? (https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.svd.html)
    # if effWC_components.shape[0] is not effWC_components.shape[1]:
    #     print('Error ! Need square T matrix for now...'); exit(1)
    assert(effWC_components.shape[0] == effWC_components.shape[1]), 'Error ! Need square T matrix for now...'
    if DEBUG_: print(benchmark_weights.shape); print(effWC_components.shape); print(np.linalg.inv(effWC_components).shape)

    fit_coeffs = np.dot(np.linalg.inv(effWC_components), np.transpose(benchmark_weights) ) #Matrix inversion ; A = T^-1 . w
    fit_coeffs = np.transpose(fit_coeffs) #Transpose to get desired shape (n_events, 0
    if DEBUG_:
        print('\n-- Get_FitCoefficients'); print(fit_coeffs.shape)
        if DEBUG_>1:
            # print(fit_coeffs)
            print('fit_coeffs', fit_coeffs[0])
            print('benchmark_weights', benchmark_weights[0])
            print('WC', effWC_components[0])
            print('inv WC', np.linalg.inv(effWC_components)[0])
            print('benchweight0 : coeff*WC0=', np.dot(effWC_components[0],fit_coeffs[0]) )
            print('benchweight1 : coeff*WC1=', np.dot(effWC_components[1],fit_coeffs[0]) )
            print('benchweight2 : coeff*WC2=', np.dot(effWC_components[2],fit_coeffs[0]) )
    return fit_coeffs

# //--------------------------------------------
# //--------------------------------------------

 ###### #    # ##### #####    ##   #####   ####  #
 #       #  #    #   #    #  #  #  #    # #    # #
 #####    ##     #   #    # #    # #    # #    # #
 #        ##     #   #####  ###### #####  #    # #      ###
 #       #  #    #   #   #  #    # #      #    # #      ###
 ###### #    #   #   #    # #    # #       ####  ###### ###

def Extrapolate_EFTweight(effWC_components, fit_coeffs):
    """
    Compute (extrapolate) new event weights, for all (n_events) considered events and all (n_points) new EFT points

    Parameters:
    effWC_components (ndarray of shape [n_points, n_components]) : 'effective WC' of given component, for all (n_points) EFT points
    fit_coeffs (ndarray of shape [n_events, n_components]) : fit coefficients associated with the components, for all (n_points) EFT points

    Returns:
    extrapolWeights (ndarray of shape [n_events, n_points]) : new weights extrapolated at the considered EFT points, for each considered event
    """

    if fit_coeffs.ndim==2 and effWC_components.ndim==2 and fit_coeffs.shape[1] is not effWC_components.shape[1]:
        print(fit_coeffs.shape[1])
        print(effWC_components.shape[1])
        print(colors.fg.red, 'Error : Arrays fit_coeffs and effWC_components have different numbers of columns !', colors.reset)
        return

    # print(effWC_components.shape); print(fit_coeffs.shape)
    extrapolWeights = np.dot(effWC_components, np.transpose(fit_coeffs)) #w' = t . A^T
    extrapolWeights = np.transpose(extrapolWeights) #Transpose to get desired shape (n_events, n_points)
    # print(extrapolWeights.shape)

    if DEBUG_:
        print('\n-- Extrapolate_EFTweight'); print(extrapolWeights.shape)
        if DEBUG_>1: print(extrapolWeights)
    return extrapolWeights

# //--------------------------------------------
# //--------------------------------------------

def Extrapolate_EFTxsec(effWC_components, fit_coeffs):
    """
    Idem as Extrapolate_EFTweight(), but sum all events to get cross sections.

    Returns:
    extrapolXsecs (ndarray of shape [n_points]) : new cross sections extrapolated at the considered EFT points
    """

    extrapolWeights = Extrapolate_EFTweight(effWC_components, fit_coeffs) #Get the extrapolated weights
    extrapolXsecs = np.sum(extrapolWeights, axis=0)

    if DEBUG_:
        print('\n-- Extrapolate_EFTxsec'); print(extrapolXsecs.shape)
        if DEBUG_>1: print(extrapolXsecs)
    return extrapolXsecs

# //--------------------------------------------
# //--------------------------------------------

 ###### #    # ##### ###### #    # #####
 #       #  #    #   #      ##   # #    #
 #####    ##     #   #####  # #  # #    #
 #        ##     #   #      #  # # #    #
 #       #  #    #   #      #   ## #    #
 ###### #    #   #   ###### #    # #####

#NB: only extending arrays used in training (x, weight, theta). Could also resize accordingly other arrays (EFTweights, EFTweightIDs, EFT_FitCoeffs) to keep track of this per-event info, but not needed for now
def Extend_Dataset(parameterizedDNN, listOperatorsParam, labels_list, list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses, list_indexSM_allClasses, nPointsPerOp=10, minWC=-5, maxWC=5):
    """
    Extend the original data so that the DNN can be parameterized on WCs. Events need to be reweighted at many different EFT points for the DNN to learn how to interpolate between them.
    First, define the points 'theta' on which the DNN will get trained. For now, only training DNN at points where a single operator is non-zero at once. Then, for the same input features x, will duplicate events many times at different points, with different reweights.

    Parameters:
    Lists of features, weights, etc. obtained from previous functions
    nPointsPerOp, minWC, maxWC : define the interval [min, max, step] in which points will be drawn uniformly, separately for each operator ; these are the EFT points on which the DNN will get trained

    Returns:
    Same lists, after extension
    """

    if parameterizedDNN == False: #Return empty lists
        list_thetas_allClasses = []; list_targetClass_allClasses = []
        return list_x_allClasses, list_weights_allClasses, list_thetas_allClasses, list_targetClass_allClasses

    if DEBUG_>1:
        for iclass in range(len(list_x_allClasses)):
            list_x_allClasses[iclass]=list_x_allClasses[iclass][:5]; list_weights_allClasses[iclass]=list_weights_allClasses[iclass][:5]; list_EFTweights_allClasses[iclass]=list_EFTweights_allClasses[iclass][:5]; list_EFTweightIDs_allClasses[iclass]=list_EFTweightIDs_allClasses[iclass][:5]; list_EFT_FitCoeffs_allClasses[iclass]=list_EFT_FitCoeffs_allClasses[iclass][:5]

    for iclass in range(len(labels_list)):
        operatorNames, _, _ = Parse_EFTpoint_IDs(list_EFTweightIDs_allClasses[iclass][0,0]) #Asumme that benchmark points are identical for all events in the sample
        operatorNames = operatorNames[0] #Single-element list -> array
        for op1 in listOperatorsParam:
            if "PrivMC" in labels_list[iclass] and op1 not in operatorNames and op1.lower() not in operatorNames: print(colors.fg.red, 'Error : parameterized operator ', op1,' not found in class', labels_list[iclass], colors.reset, ' (Operators found : ', operatorNames, ')'); exit(1)

# //--------------------------------------------

    #-- Define the hypotheses theta on which the DNN will get trained. Draw values uniformly in given interval for each operator, translate into array. Also define the target to train on (<-> the operator which is activated at a given point)
    thetas, targetClass = Get_ThetaParameters(listOperatorsParam, operatorNames, nPointsPerOp, minWC, maxWC)
    # print(thetas)

    #-- Determine the components necessary to parameterize the event weights
    n_components, components = Find_Components(operatorNames) #Determine the components required to parameterize the event weight #NB: assumes that they are identical for all events in this process
    # print(components)

    list_x_allClasses_tmp = []
    list_weights_allClasses_tmp = []
    list_thetas_allClasses_tmp = []; list_targetClass_allClasses_tmp = []
    #Loop on process classes (which should each contain a single process for now)
    for iclass in range(len(labels_list)):

        #-- Get the 'effective WC' values scaling each fit component
        effWC_components_thetas_class = Get_EffectiveWC_eachComponents(n_components, components, thetas) #Determine the 'effective WC' values associated with each component, for each benchmark point
        # print(effWC_components_thetas_class)

        #-- Extrapolate the event weights at the new points theta, for all events
        newWeights = Extrapolate_EFTweight(effWC_components_thetas_class, list_EFT_FitCoeffs_allClasses[iclass])

        #-- Now that we have all the necessary ingredients (inputs x, original event weights, new extrapolated weights, WCs of the new points), can concatenate arrays to obtain the full training/testing data used for training and validation
        # //--------------------------------------------
        # weights_refPoint_class = np.tile(list_weights_allClasses[iclass], len(thetas)) #Half of the total events are generated at a ref. point, and the other half at all the EFT points
        weights_refPoint_class = np.tile(list_EFTweights_allClasses[iclass][:,list_indexSM_allClasses[iclass]], len(thetas)) #Half of the total events are generated at a ref. point, and the other half at all the EFT points #Choose SM as reference point (since clasifying against EFT operators)

        weights_allThetas_class = newWeights.flatten() #Shape (n,m) -> (n*m)

        x_class = np.tile(list_x_allClasses[iclass], len(thetas)) #Each point shares the same input features x. Only the event weight, and the values of the WCs change

        thetas_class = np.tile(thetas, (len(list_x_allClasses[iclass]), 1) ) #Each original event gets repeated at all the different points theta

        targetClass_class = np.tile(targetClass, (len(list_x_allClasses[iclass]), 1) ) #Keep track of which operator gets activated (for labelling)

        #Create array of same shape as targetClass_class, but with only first column (reference point) filled with ones --> Events generated at reference point (half of the total events fed to DNN) will later have label '0' (one-hot encoding of first column)
        targetClassRef_class = np.zeros((targetClass_class.shape)); targetClassRef_class[:,0] = 1
        # //--------------------------------------------

        #-- Concatenate same data twice : once with SM weights, once with reweights corresponding to EFT points to train on
        list_x_allClasses_tmp.append( np.concatenate((x_class, x_class)) )
        list_weights_allClasses_tmp.append( np.concatenate((weights_refPoint_class, weights_allThetas_class)) )
        list_thetas_allClasses_tmp.append(np.concatenate((thetas_class, thetas_class)) )
        list_targetClass_allClasses_tmp.append(np.concatenate((targetClassRef_class, targetClass_class)) )
        # print(list_x_allClasses_tmp[iclass].shape); print(list_weights_allClasses_tmp[iclass].shape); print(list_thetas_allClasses_tmp[iclass].shape); print(list_targetClass_allClasses_tmp[iclass].shape)

    return list_x_allClasses_tmp, list_weights_allClasses_tmp, list_thetas_allClasses_tmp, list_targetClass_allClasses_tmp

# //--------------------------------------------
# //--------------------------------------------

  ####  ###### #####    ##### #    # ###### #####   ##
 #    # #        #        #   #    # #        #    #  #
 #      #####    #        #   ###### #####    #   #    #
 #  ### #        #        #   #    # #        #   ######
 #    # #        #        #   #    # #        #   #    #
  ####  ######   #        #   #    # ######   #   #    #

def Get_ThetaParameters(listOperatorsParam, operatorNames, n_pointsPerOperator, minVal, maxVal):

    list_thetas_allOperators = []
    list_targetClass = []
    for i_op_ToParameterize in range(len(listOperatorsParam)):
        # print(i_op_ToParameterize)

        thetas = np.zeros((n_pointsPerOperator, len(operatorNames))) #Shape (n_operators_inSample, n_pointsPerOperator)
        targetClass = np.zeros((n_pointsPerOperator, len(operatorNames)+1)) #Keep track of which operator is activated in each element of theta ; necessary to later specify target for DNN classification/regression #Shape (n_operators_inSample, n_pointsPerOperator) #Additional dim represents SM class (not filled here)

        for i_opInSample in range(len(operatorNames)):
            # print(i_opInSample)
            # print(operatorNames[i_opInSample])
            # print(listOperatorsParam[i_op_ToParameterize])
            # print(listOperatorsParam[i_op_ToParameterize].lower())

            if listOperatorsParam[i_op_ToParameterize] == operatorNames[i_opInSample] or listOperatorsParam[i_op_ToParameterize].lower() == operatorNames[i_opInSample]:

                iter = 0
                for x in np.linspace(minVal, maxVal, num=n_pointsPerOperator+1):

                    if x == 0: continue #skip SM point
                    # print('x', x)

                    thetas[iter, i_opInSample] = x
                    targetClass[iter, i_opInSample+1] = 1
                    iter+= 1

        list_thetas_allOperators.append(thetas)
        list_targetClass.append(targetClass)

    thetas_allOperators = np.concatenate(list_thetas_allOperators)
    targetClass = np.concatenate(list_targetClass)

    if DEBUG_:
        print('\n-- Get_ThetaParameters'); print(thetas_allOperators)
        if DEBUG_>1: print(thetas_allOperators)
    return thetas_allOperators, targetClass

# //--------------------------------------------
# //--------------------------------------------

# //--------------------------------------------
# //--------------------------------------------

# //--------------------------------------------
# //--------------------------------------------

# //--------------------------------------------
# //--------------------------------------------

















# //--------------------------------------------
