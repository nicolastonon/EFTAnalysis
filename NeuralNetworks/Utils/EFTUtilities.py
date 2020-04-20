# Author : Nicolas Tonon, DESY
# Helper functions related to the EFT treatment in the DNN (parameterization of weight .vs. WCs, extrapolation, computation of score/likelihood, etc.

DEBUG_ = False #True <-> debug printouts

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

import warnings
warnings.filterwarnings("ignore", message="FutureWarning: in the future insert will treat boolean arrays and array-likes as boolean index instead of casting it to integer")


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
    Parse an array of strings, each representing

    Parameters:
    benchmarkIDs (ndarray of shape [n_points]) : array of strings, each corresponding to a separate EFT point. Example : 'rwgt_ctZ_5_ctW_3'

    Returns:
    operatorNames (ndarray of shape [n_points, n_operators]) : array of names of EFT operators, for all (n_points) EFT points
    operatorWCs (ndarray of shape [n_points, n_operators]) : array of corresponding WC values, for all (n_points) EFT points
    """

    operatorNames = []
    operatorWCs = []

    prefix = 'rwgt_'
    for ID in benchmarkIDs:

        #For each EFT point, get the list of names/WCs for all operators
        list_operatorNames = []
        list_operatorWCs = []

        if not ID.startswith("rwgt_c"): #Every considered EFT operator is expected to start with this substring ; for others (e.g. 'rwgt_sm'), don't parse
            list_operatorNames.append('xxx') #Operator name
            list_operatorWCs.append(float(0)) #Operator WC
            continue

        ID = CheckName_EFTpoint_ID(ID) #Enforce proper naming convention

        if ID.startswith(prefix): ID = ID[len(prefix):] #Naming convention, strip this substring
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

    if DEBUG_: print(operatorNames.shape); print(operatorWCs.shape); print(operatorNames); print(operatorWCs)
    return operatorNames, operatorWCs

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
    Determine the number of components necessary to compute the squared matrix element, and at what power each EFT operator enters each component

    Parameters:
    operatorNames (ndarray of shape [n_operators]) : names of the EFT operators

    Returns:
    n_components (int) : number of components
    components (ndarray of shape [n_components, n_operators]) : array whose elements represent the power at which an operator 'enters' (scales?) a component
    """

    minPower_perOperator = 0
    maxPower_perOperator = 2
    if len(operatorNames) == 1: #Assumption : 1 operator <-> pure EFT
        minPower_perOperator = 2
        maxPower_perOperator = 2

    #Set the min/max powers at which each operator can enter a component
    parameter_max_power = []
    for i in range(len(operatorNames)):
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
                # print("Adding component %s", powers)
            # else:
                # print("Not adding component %s again", powers)

    components = np.array(components, dtype=np.int) #Convert list to array
    n_components = len(components)

    if DEBUG_: print('n_components', n_components); print(components.shape); print(components)
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
            for p in range(len(np.transpose(operatorWCs))): #For each operator (transpose so that 1rst dim represent operators)
                if DEBUG_: print('x', x, 'y', y, 'p', p)
                factor *= float(operatorWCs[y, p] ** components[x, p]) #Contribution of this operator to this component
                # print(operatorWCs[p]); print(components[x, p]); print(factor)
                effWC_components[y, x] = factor

    if DEBUG_: print(effWC_components.shape); print(effWC_components)
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
    effWC_components (ndarray of shape [n_points, n_components]) : 'effective WC' of given component, for all (n_points) EFT points
    benchmark_weights (ndarray of shape [n_components, n_operators]) : array whose elements represent the power at which an operator 'enters' (scales?) a component

    Returns:
    fit_coeffs (ndarray of shape [n_points, n_components]) : fit coefficients associated with the components, for all (n_points) EFT points
    """

    #Remove bad strings

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

    fit_coeffs = np.dot(benchmark_weights, np.linalg.inv(effWC_components)) #Matrix inversion

    if DEBUG_: print(fit_coeffs.shape); print(fit_coeffs)
    return fit_coeffs

# //--------------------------------------------
# //--------------------------------------------

 ###### #    # ##### #####    ##   #####   ####  #
 #       #  #    #   #    #  #  #  #    # #    # #
 #####    ##     #   #    # #    # #    # #    # #
 #        ##     #   #####  ###### #####  #    # #      ###
 #       #  #    #   #   #  #    # #      #    # #      ###
 ###### #    #   #   #    # #    # #       ####  ###### ###

def Get_Extrapolated_EFTweight(effWC_components, fit_coeffs):
    """
    Compute (extrapolate) new event weights, for all (n_events) considered events and all (n_points) new EFT points

    Parameters:
    effWC_components (ndarray of shape [n_points, n_components]) : 'effective WC' of given component, for all (n_points) EFT points
    fit_coeffs (ndarray of shape [n_events, n_components]) : fit coefficients associated with the components, for all (n_points) EFT points

    Returns:
    extrapolWeights (ndarray of shape [n_events, n_points]) : new weights extrapolated at the considered EFT points, for each considered event
    """

    # print(effWC_components.shape); print(fit_coeffs.shape)
    extrapolWeights = np.dot(fit_coeffs, np.transpose(effWC_components)) #FIXME -- check

    # print(extrapolWeights)
    return extrapolWeights

# //--------------------------------------------
# //--------------------------------------------


# //--------------------------------------------
# //--------------------------------------------


# //--------------------------------------------
# //--------------------------------------------


# //--------------------------------------------
# //--------------------------------------------


# //--------------------------------------------
# //--------------------------------------------

















# //--------------------------------------------
