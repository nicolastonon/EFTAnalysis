# Author : Nicolas Tonon, DESY
# Helper functions related to the EFT treatment in the DNN (parameterization of weight .vs. WCs, extrapolation, computation of score/likelihood, etc.

'''

#TODO :
- func : get baseline points and reweights for all events ?
- func : args = baseline weights, baseline points ; matrix inversion to get coeff matrix (for all events)
- func : args = EFT point, coeff matrix ; compute and return extrapolated weight

# NOTES :
- Only applies to private EFT samples
- Invert matrices for all events at once ?

# LIMITATIONS / ASSUMPTIONS :
- Can not sum several processes into a process class including an EFT process. Each EFT process must constitute a separate class.
- Assume that a given EFT sample has the exact same list of reweight points (in same order) for all considered years (could change that in the future by estimating fit coefficients directly when reading the sample, separately for each year)
- If single operator identified, assume sample is pure-EFT only

# NAMING CONVENTIONS :
- 'WC' = Wilson coefficient <-> the scalar value scaling a given EFT operator in the Lagrangian
- 'Parameter'/'Operator' = EFT operator
- 'Benchmark point' = EFT point (<-> given WC values for given EFT operators) at which the event weight has been evaluated by MG. The ensemble of the benchmark points is also refered to as the 'basis'
- 'Benchmark weight' = corresponding weight value (evaluated by MG with the reweighting module)
- 'Benchmark ID' = string encoding of the name of the benchmark point
- 'Components' = individual terms contributing to the squared matrix element of the process. For example, if considering a single operator, there are 3 components to account for (corresponding to the SM squared term, the interference term, and the pure-EFT term)
- ...

'''

# from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
# from root_numpy import root2array, tree2array, array2root
import numpy as np
import pandas as pd
from pathlib import Path
import itertools

# //--------------------------------------------

# //--------------------------------------------
# //--------------------------------------------

def Parse_EFTpoint_IDs(benchmarkIDs):
    """
    Parse an array of strings, each representing

    Parameters:
    benchmarkIDs (ndarray) : array of strings, each corresponding to a separate EFT point. Example : 'rwgt_ctZ_5_ctW_3'

    Returns:
    operatorNames (ndarray of shape [n_points, n_operators]) : array of names of EFT operators, for all (n_points) EFT points
    operatorWCs (ndarray of shape [n_points, n_operators]) : array of corresponding WC values, for all (n_points) EFT points
    """

    operatorNames = []
    operatorWCs = []

    prefix = 'rwgt_'
    for ID in benchmarkIDs:
        if ID.startswith(prefix): #Naming convention, strip this substring
            ID = ID[len(prefix):]
        else:
            print('Error : naming convention in benchmark ID not recognized'); exit(1)

        list_keys = ID.split('_') #Split string into list of substrings (pairs [name,WC])

        #For each EFT point, get the list of names/WCs for all operators
        list_operatorNames = []
        list_operatorWCs = []
        for ikey in range(0, len(list_keys)-1, 2):
            list_operatorNames.append(list_keys[ikey]) #Operator name
            list_operatorWCs.append(float(list_keys[ikey+1])) #Operator WC

        #Append list for each EFT point
        operatorNames.append(list_operatorNames)
        operatorWCs.append(list_operatorWCs)

    #Convert list of lists into 2D array. Each list element must have an equal length
    operatorNames = np.array(operatorNames)
    operatorWCs = np.array(operatorWCs)

    # print(operatorNames.shape); print(operatorWCs.shape); print(operatorNames); print(operatorWCs)
    return operatorNames, operatorWCs

# //--------------------------------------------
# //--------------------------------------------

#Adapted from madminer code (https://github.com/diana-hep/madminer/blob/080e7fefc481bd6aae16ba094af0fd6bfc301dff/madminer/utils/morphing.py#L116)
#NB : row full of 0 <-> pure SM ; any row with sum of powers < 2 <-> interference with SM
def Find_Components(operatorNames):
    """
    Determine the number of components necessary to compute the squared matrix element, and at what power each EFT operator enters each component

    Parameters:
    operatorNames (ndarray of shape [n_operators]) : names of the EFT operators

    Returns:
    n_components (int) : number of components
    components (ndarray of shape [n_components, n_operators]) : array whose elements represent the power at which an operator 'enters' (scales?) a component
    """

    maxPower_perOperator = 2
    if len(operatorNames) == 1: #Assumption : 1 operator <-> pure EFT
        maxPower_perOperator = 1

    #Set the maximal power at which each operator can enter a component
    parameter_max_power = []
    for i in range(len(operatorNames)):
        parameter_max_power.append(maxPower_perOperator)

    components = []
    for i in range(len(operatorNames)):

        #Get all arrays corresponding to all possible combinations of operator/power
        powers_each_component = [range(max_power + 1) for max_power in parameter_max_power]
        for powers in itertools.product(*powers_each_component):
            powers = np.array(powers, dtype=np.int)

            if np.sum(powers) > maxPower_perOperator: #Must not exceed maximal total EFT power
                continue

            if not any((powers == x).all() for x in components): #If combination not yet included, add it
                components.append(np.copy(powers))
                # print("Adding component %s", powers)
            # else:
                # print("Not adding component %s again", powers)

    components = np.array(components, dtype=np.int) #Convert list to array
    n_components = len(components)

    # print('n_components', n_components); print(components.shape); print(components)
    return n_components, components

# //--------------------------------------------
# //--------------------------------------------

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
            for p in range(len(operatorWCs)): #For each operator
                factor *= float(this_basis[y, p] ** components[x, p]) #Contribution of this operator to this component
                # print(operatorWCs[p]); print(components[x, p]); print(factor)
                effWC_components[y, x] = factor

    # print(effWC_components.shape); print(effWC_components)
    return effWC_components

# //--------------------------------------------
# //--------------------------------------------

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

    #Need square matrix, i.e. as many benchmark points as there are components
    #  SVD ? (https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.svd.html)
    # if effWC_components.shape[0] is not effWC_components.shape[1]:
    #     print('Error ! Need square T matrix for now...'); exit(1)
    assert(effWC_components.shape[0] == effWC_components.shape[1]), 'Error ! Need square T matrix for now...'

    # print(benchmark_weights.shape); print(effWC_components.shape); print(np.linalg.inv(effWC_components).shape)
    fit_coeffs = np.dot(benchmark_weights, np.linalg.inv(effWC_components)) #Matrix inversion

    # print(fit_coeffs.shape); print(fit_coeffs)
    return fit_coeffs

# //--------------------------------------------
# //--------------------------------------------

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
