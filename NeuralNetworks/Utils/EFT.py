# Helper functions related to the EFT treatment in the NN (parameterization of weight .vs. WCs, extrapolation, computation of score/likelihood, etc.

DEBUG_ = 0 #0: no debug printouts; 1: minimal printouts; 2: maximal printouts (arrays, etc.)

# //--------------------------------------------

'''
# NOTES :
-

# LIMITATIONS / ASSUMPTIONS :
- Can not sum several processes into a process class including an EFT process. Each EFT process must constitute a separate class.
- Assume that a given EFT sample has the exact same list of reweight points (in same order) for all considered years (could change that in the future by estimating fit coefficients directly when reading the sample, separately for each year)
- If single operator identified, assume sample is pure-EFT only
- Can only rescale array of event if all events have exact same EFT parameterization
- Negative weights are ignored

# NAMING CONVENTIONS :
NB: Squared amplitude M^2 = a0 + a1.c1 + a2.c2, for a single EFT operator. The 3 terms correspond to SM, interference, and pure-EFT.
- 'WC' = Wilson coefficient <-> the scalar value scaling a given EFT operator in the Lagrangian
- 'Parameter'/'Operator' = EFT operator
- 'EFT point'/'point'/'theta' = point in EFT phase space, corresponding to given WC values of considered operators
- 'Benchmark point' = EFT point at which the event weight has been evaluated by MG. The ensemble of the benchmark points is also refered to as the 'basis'
- 'Benchmark weight' = corresponding weight value (evaluated by MG with the reweighting module)
- 'Benchmark ID' = string encoding of the name of the benchmark point
- 'Components' = individual terms contributing to the squared matrix element of the process. For example, if considering a single operator, there are 3 components to account for (corresponding to the SM squared term, the interference term, and the pure-EFT term)
- 'Effective WC'/'component weight' = factors 'c_i' for all components. Example: for the interfence term between ctZ=3 and ctW=5, c=3*5=15. They are not the WC values of the operators, but depend directly on them, hence labelled "effective WC"
- 'Fit coefficients'/'coefficients' = factors 'a_i' scaling the components ; these are the coefficients determined per-event from the benckmark points, which are then used to parameterize the event weight
- 'JLR' = joint likelihood ratio, denoted r. Along with score t, these variables correspond to the augmented data extracted from the generator (see reference summary article: https://arxiv.org/abs/1805.00020)
- ...

'''

# //--------------------------------------------

# from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
# from root_numpy import root2array, tree2array, array2root
import numpy as np
import pandas as pd
from pathlib import Path
import itertools
from Utils.Helper import *
from Utils.ColoredPrintout import colors

import warnings
warnings.filterwarnings("ignore", message="in the future insert will treat boolean arrays and array-likes as boolean index instead of casting it to integer")

# //--------------------------------------------
# //--------------------------------------------

  ####   ####  #    # #####   ####  #    # ###### #    # #####  ####
 #    # #    # ##  ## #    # #    # ##   # #      ##   #   #   #
 #      #    # # ## # #    # #    # # #  # #####  # #  #   #    ####
 #      #    # #    # #####  #    # #  # # #      #  # #   #        #
 #    # #    # #    # #      #    # #   ## #      #   ##   #   #    #
  ####   ####  #    # #       ####  #    # ###### #    #   #    ####

#Adapted from madminer code (https://github.com/diana-hep/madminer/blob/080e7fefc481bd6aae16ba094af0fd6bfc301dff/madminer/utils/morphing.py#L116)
#NB: row full of zeros (<-> not impacted by any operator) corresponds to SM (will correspond to a fix coefficient value of 0^0=1)
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

def Get_EffectiveWC_eachComponent(n_components, components, operatorWCs):
    """
    Determine the 'effective WC' (see notes at b.o.f.) for each component, i.e. the product of the corresponding WC values

    Parameters:
    n_components (int) : number of components
    components (ndarray of shape [n_components, n_operators]) : array whose elements represent the power at which an operator 'enters' (scales?) a given component
    operatorWCs (ndarray of shape [n_points, n_operators]) : array whose columns represent the WC value of each operator defining a given point, and rows represent different EFT points

    Returns:
    effWC_components (ndarray of shape [n_points, n_components]) : 'effective WC' of given component, for all (n_points) EFT points
    """

    operatorWCs = np.atleast_2d(operatorWCs) #Need 2D array

    effWC_components = np.zeros((len(operatorWCs), n_components)) #Create empty array of shape (n_points, n_components)

    for y in range(len(operatorWCs)): #For each EFT point
        for x in range(n_components): #For each component
            effWC_components[y,x] = np.prod(operatorWCs[y,:] ** components[x,:]) #Product of contributions of all operators to this component (e.g. if consider 2 operators, element corresponding to their interference will have value: WC1*WC2; element corresponding to squared operator1 will have value: WC1*WC1 ; etc.)

            # factor = 1.0
            # for p in range(len(np.transpose(operatorWCs))): #For each operator (transpose so that 1st dim corresponds to operators)
            #     factor *= float(operatorWCs[y, p] ** components[x, p]) #Contribution of this operator to this component
                # if DEBUG_>1: print('x', x, 'y', y, 'p', p)
                # print(operatorWCs[y, p]); print(components[x, p]); print(factor)

    if DEBUG_:
        print('\n-- Get_EffectiveWC_eachComponent'); print(effWC_components.shape)
        if DEBUG_>1: print(effWC_components)
    return effWC_components

# //--------------------------------------------
# //--------------------------------------------

  ####  #####    ##   #####  # ###### #    # #####  ####
 #    # #    #  #  #  #    # # #      ##   #   #   #
 #      #    # #    # #    # # #####  # #  #   #    ####
 #  ### #####  ###### #    # # #      #  # #   #        #
 #    # #   #  #    # #    # # #      #   ##   #   #    #
  ####  #    # #    # #####  # ###### #    #   #    ####

#Adapted from MadMiner code: https://github.com/diana-hep/madminer/blob/43d5b59c6857b56cfbd2d06f683ccf32a0aa5440/madminer/utils/morphing.py#L462-L530
def Get_Gradients_EffectiveWC_eachComponent(n_components, components, operatorWCs):
    """
    Determine the *gradients* of the 'effective WC' (see notes at b.o.f.) for each component, for all points. One gradient w.r.t. each considered operator.

    Parameters:
    n_components (int) : number of components
    components (ndarray of shape [n_components, n_operators]) : array whose elements represent the power at which an operator 'enters' (scales?) a given component
    operatorWCs (ndarray of shape [n_points, n_operators]) : array whose columns represent the WC value of each operator defining a given point, and rows represent different EFT points

    Returns:
    gradients_effWC_components (ndarray of shape [n_operators, n_points, n_components]) : 'effective WC' of given component, for all (n_points) EFT points
    """

    gradients_effWC_components = np.zeros((components.shape[1], len(operatorWCs), n_components)) #Create empty array of shape (n_operators, n_points, n_components)

    for z in range(components.shape[1]): #For each operator (<-> loop on array elements. Will differentiate M^2 w.r.t. this operator)

        gradient_effWC_components_operator = np.zeros((len(operatorWCs), n_components)) #Create empty sub-array of shape (n_points, n_components)
        for y in range(len(operatorWCs)): #For each EFT point
            for x in range(n_components): #For each component (entering the reweighting parameterization)

                factor = 1.0
                for p in range(components.shape[1]): #For each operator (loop on operators, get the product of their individual impacts on the considered component)

                    if p == z and components[x,p] >= 1: #If current component (x,p) contains a contribution from the current operator (>=1), *and this is the operator w.r.t. which we are computing the gradient* (p==z) --> (x^n)' = n*(x^n-1)
                        factor *= float(components[x,p] * operatorWCs[y,p] ** (components[x,p]-1))
                    elif p == z: #Else if current component (x,p) *does not* contain a contribution from the current operator, *and this is the operator w.r.t. which we are computing the gradient* (p==z) --> (x^n)' = 0
                        factor = 0.0; break
                    else: #Else if current component (x,p) contains a contribution from the current operator (>=1), and this is *not* the operator w.r.t. which we are computing the gradient (p!=z) --> Multiplies total factor by constant x^n (NB: but if the current element does not depend at all on operator w.r.t. which we are differentiating, the total factor will eventually be null)
                        factor *= float(operatorWCs[y,p] ** components[x,p])

                    # if z==0 and y==0:
                    #     print('z', z, 'y', y, 'x', x, 'p', p)
                    #     print(operatorWCs[y,p]); print(components[x,p]); print(factor)
                    #     if DEBUG_>1: print('z', z, 'y', y, 'x', x, 'p', p)

                gradient_effWC_components_operator[y,x] = factor #Fill sub-array

        gradients_effWC_components[z, ...] = gradient_effWC_components_operator #Fill total array with sub-arrays

    # print('\n-- Get_Gradients_EffectiveWC_eachComponent');
    # print(gradients_effWC_components.shape)
    # print('z=0,p=0', operatorWCs[0,:], '\n', components[:,:], '\n', gradients_effWC_components[0,0,:])
    # print('z=0,p=1', operatorWCs[1,:], '\n', components[:,:], '\n', gradients_effWC_components[0,1,:])

    # if DEBUG_:
    #     print('\n-- Get_Gradients_EffectiveWC_eachComponent'); print(gradients_effWC_components.shape)
    #     if DEBUG_>1: print(gradients_effWC_components)
    return gradients_effWC_components

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
            print('fit_coeffs0', fit_coeffs[0])
            print('benchmark_weights0', benchmark_weights[0])
            print('effWC0', effWC_components[0])
            print('inv effWC0', np.linalg.inv(effWC_components)[0])
            print('benchweight0 : coeff*WC0=', np.dot(effWC_components[0],fit_coeffs[0]) )
            print('benchweight1 : coeff*WC1=', np.dot(effWC_components[1],fit_coeffs[0]) )
            print('benchweight2 : coeff*WC2=', np.dot(effWC_components[2],fit_coeffs[0]) )
    return fit_coeffs

# //--------------------------------------------
# //--------------------------------------------

 ###### #    # ##### #####    ##   #####   ####  #        ##   ##### ######
 #       #  #    #   #    #  #  #  #    # #    # #       #  #    #   #
 #####    ##     #   #    # #    # #    # #    # #      #    #   #   #####
 #        ##     #   #####  ###### #####  #    # #      ######   #   #
 #       #  #    #   #   #  #    # #      #    # #      #    #   #   #
 ###### #    #   #   #    # #    # #       ####  ###### #    #   #   ######

def Extrapolate_EFTweights(effWC_components, fit_coeffs):
    """
    Compute ('extrapolate') new event weights, for all (n_events) considered events and all (n_points) new EFT points, by multiplying the 'fit coefficients' of the events with the 'effective WC values' defining the new points

    Parameters:
    effWC_components (ndarray of shape [n_points, n_components]) : 'effective WC' of given component, for all (n_points) EFT points
    fit_coeffs (ndarray of shape [n_events, n_components]) : fit coefficients associated with the components, for all (n_points) EFT points

    Returns:
    extrapolatedWeights (ndarray of shape [n_events, n_points]) : new weights extrapolated at the considered EFT points, for each considered event
    """

    if fit_coeffs.ndim==2 and effWC_components.ndim==2 and fit_coeffs.shape[1] != effWC_components.shape[1]:
        print(fit_coeffs.shape[1]); print(effWC_components.shape[1])
        print(colors.fg.red, 'Error : Arrays fit_coeffs and effWC_components have different numbers of columns !', colors.reset)
        return

    # print(effWC_components.shape); print(fit_coeffs.shape)
    extrapolatedWeights = np.dot(effWC_components, np.transpose(fit_coeffs)) #w' = t . A^T
    extrapolatedWeights = np.transpose(extrapolatedWeights) #Transpose to get desired shape (n_events, n_points)
    # print(extrapolatedWeights.shape)

    if DEBUG_:
        print('\n-- Extrapolate_EFTweights'); print(extrapolatedWeights.shape)
        if DEBUG_>1: print(extrapolatedWeights)
    return extrapolatedWeights

# //--------------------------------------------
# //--------------------------------------------

def Extrapolate_EFTxsecs(effWC_components, fit_coeffs):
    """
    Same as Extrapolate_EFTweights(), but sum all events to get cross sections.
    NB: if extrapolated weights were already computed (e.g. because they need to be stored), call Extrapolate_EFTxsecs_fromWeights() instead

    Returns:
    extrapolatedXsecs (ndarray of shape [n_points]) : new cross sections extrapolated at the considered EFT points
    """

    extrapolatedWeights = Extrapolate_EFTweights(effWC_components, fit_coeffs) #Get the extrapolated weights
    extrapolatedXsecs = np.sum(extrapolatedWeights, axis=0)

    if DEBUG_:
        print('\n-- Extrapolate_EFTxsecs'); print(extrapolatedXsecs.shape)
        if DEBUG_>1: print(extrapolatedXsecs)
    return extrapolatedXsecs

# //--------------------------------------------
# //--------------------------------------------

def Extrapolate_EFTxsecs_fromWeights(extrapolatedWeights):
    """
    Sum extrapolated weights, at each point, to obtain extrapolated cross sections at all points

    Parameters:
    extrapolatedWeights (ndarray of shape [n_events, n_points]) : new weights extrapolated at the considered EFT points, for each considered event

    Returns:
    extrapolatedXsecs (ndarray of shape [n_points]) : new cross sections extrapolated at the considered EFT points
    """

    extrapolatedXsecs = np.sum(extrapolatedWeights, axis=0)

    if DEBUG_:
        print('\n-- Extrapolate_EFTxsecs_fromWeights'); print(extrapolatedXsecs.shape)
        if DEBUG_>1: print(extrapolatedXsecs)
    return extrapolatedXsecs

# //--------------------------------------------
# //--------------------------------------------

def Extrapolate_SM_Weights(n_components, components, fit_coeffs, listOperatorsParam):
    """
    Extrapolate event weights at SM point.
    """

    WCs_SM = np.zeros(len(listOperatorsParam)) #SM <-> All operators are set to 0
    effWC_components_SM = Get_EffectiveWC_eachComponent(n_components, components, WCs_SM) #Get corresponding matrix
    weight_SM = Extrapolate_EFTweights(effWC_components_SM, fit_coeffs) #Get corresponding event weights
    weight_SM = np.squeeze(weight_SM) #2D -> 1D


    weight_SM = Extrapolate_EFTweights(effWC_components_SM, fit_coeffs)

    return weight_SM

# //--------------------------------------------
# //--------------------------------------------

def Extrapolate_SM_Xsec(opts, weights_SM):
    """
    Extrapolate cross section at SM point.
    """

    effWC_components_SM = np.zeros(len(opts["listOperatorsParam"])) #SM <-> All operators are set to 0

    xsec_SM = extrapolatedXsecs = np.sum(weights_SM, axis=0)

    return xsec_SM

# //--------------------------------------------
# //--------------------------------------------

  ####  ###### #####    ##### #    # ###### #####   ##
 #    # #        #        #   #    # #        #    #  #
 #      #####    #        #   ###### #####    #   #    #
 #  ### #        #        #   #    # #        #   ######
 #    # #        #        #   #    # #        #   #    #
  ####  ######   #        #   #    # ######   #   #    #

#FIXME : different thetas for strategies. different y_targetClass
def Get_ThetaParameters(opts, operatorNames):
    """
    Sets the grid of EFT points on which the NN will get trained (<-> at which training events will be extrapolated).
    Returned arrays are defined in terms of operators found in the sample (<-> on which event weights are parameterized), not only the operators selected by the user. This is because these arrays will have to be used to parameterize weights, hence need to follow sample's convention.

    Parameters:
    listOperatorsParam : subset of operators on which to train the NN
    operatorNames : all operators found in the sample
    nPointsPerOperator, minWC, maxWC : define the interval [min, max, step] in which points will be drawn uniformly, separately for each operator ; these are the EFT points on which the NN will get trained

    Returns:
    thetas_allOperators (ndarray of shape [n_points, n_operators]) : WC values of all operators, for each point included in NN training
    targetClass (ndarray of shape [n_operators]) : translates each point into a target class (for event labelling). '1' means that the corresponding operator is non-zero, '0' otherwise. 1 column for each operator included in sample (remove un-necessary operators at later step)
    """

    if opts["strategy"] == "CARL_singlePoint": return np.array([]), np.array([])

    #Read options
    listOperatorsParam = opts["listOperatorsParam"]
    nPointsPerOperator = opts["nPointsPerOperator"]
    minWC = opts["minWC"]
    maxWC = opts["maxWC"]

    list_thetas_allOperators = []
    list_targetClass = []

#--- Sample thetas uniformly between [min;max], *independently for each operator*. No more than 1 operator activated (non-zero) at once, to qllow for unambiguous labelling
    if opts["strategy"] == "CARL_multiclass":

        for i_op_ToParameterize in range(len(listOperatorsParam)): #Loop on operators selected by user
            # print(i_op_ToParameterize)

            thetas = np.zeros((nPointsPerOperator, len(operatorNames))) #Shape (n_operators_inSample, nPointsPerOperator)
            targetClass = np.zeros((nPointsPerOperator, opts["nofOutputNodes"])) #CARL (SM vs EFT) <-> all EFT points have class ID = 1 ; CARL_multiclass (SM vs op1 vs op2 vs ...) <-> one-hot encoded multilabels. 0=SM (not filled here), 1=operator1 activated, 2=operator2 activated, etc. #NB: in multiclass, for ease, assign 1 column per operator found in sample for now. Later, will remove zero-only columns (<-> only keep operators chosen by user)

            for i_opInSample in range(len(operatorNames)): #1 column per operator found in sample for now (not only operators selected by user)
                # print(i_opInSample)
                # print(operatorNames[i_opInSample])
                # print(listOperatorsParam[i_op_ToParameterize])
                # print(listOperatorsParam[i_op_ToParameterize].lower())

                if listOperatorsParam[i_op_ToParameterize] == operatorNames[i_opInSample] or listOperatorsParam[i_op_ToParameterize].lower() == operatorNames[i_opInSample]:

                    iter = 0
                    for x in np.linspace(minWC, maxWC, num=nPointsPerOperator+1):

                        if x == 0: continue #skip SM point
                        # print('x', x)

                        thetas[iter, i_opInSample] = x
                        if opts["nofOutputNodes"] == 1: targetClass[iter] = 1 #Single column <-> any EFT point has label 1 (only SM point has label 0)
                        else: targetClass[iter, i_opInSample+1] = 1 #Multi-column
                        iter+= 1

            list_thetas_allOperators.append(thetas)
            list_targetClass.append(targetClass)

        thetas_allOperators = np.concatenate(list_thetas_allOperators)
        targetClasses_allOperators = np.concatenate(list_targetClass)

#-- Sample thetas randomly betweem [min;max] in multi-dimensional EFT phase space
    else:

        rng = np.random.default_rng() #Init random generator

        thetas_allOperators = np.random.uniform(low=minWC, high=maxWC, size=(nPointsPerOperator*len(listOperatorsParam), len(operatorNames)) )
        targetClasses_allOperators = np.ones(len(thetas_allOperators)) #Binary label 0=reference point (SM, not filled here), 1=any EFT point

        for i_opInSample in range(len(operatorNames)): #1 column per operator found in sample for now (not only operators selected by user)
            if operatorNames[i_opInSample] not in listOperatorsParam: thetas_allOperators[:,i_opInSample] = 0 #Set columns corresponding to operators present in sample but not selected by user to 0

# //--------------------------------------------

    if DEBUG_:
        print('\n-- Get_ThetaParameters'); print(thetas_allOperators)
        if DEBUG_>1: print(thetas_allOperators)
    return thetas_allOperators, targetClasses_allOperators

# //--------------------------------------------
# //--------------------------------------------

   ##   #    #  ####  #    # ###### #    # ##### ###### #####
  #  #  #    # #    # ##  ## #      ##   #   #   #      #    #
 #    # #    # #      # ## # #####  # #  #   #   #####  #    #
 ###### #    # #  ### #    # #      #  # #   #   #      #    #
 #    # #    # #    # #    # #      #   ##   #   #      #    #
 #    #  ####   ####  #    # ###### #    #   #   ###### #####

 #####    ##   #####   ##
 #    #  #  #    #    #  #
 #    # #    #   #   #    #
 #    # ######   #   ######
 #    # #    #   #   #    #
 #####  #    #   #   #    #

#These 2 functions compute joint quantities with which the data can be augmented to follow the 'gold mining' approach (train NN to regress on true likelihood ratio by using optimally the reweighting informating from generator)
#More info can be found in reference summary article: https://arxiv.org/abs/1805.00020

def Compute_JointLR(weights_refPoint, xsec_refPoint, weights_thetas, xsecs_thetas):
    """
    Compute the 'joint likelihood ratio' (JLR) quantity between 2 EFT hypotheses, defined as (M^2_0/xsec_0) / (M^2_1/xsec_1). The reference hypothesis theta1 at denominator is fixed, and not trained upon. Reminder: taking advantage of the LR property: LR(hA,hB) = LR(hA,h0) / LR(hB,h0)

    Parameters:
    weights_refPoint (ndarray of shape [n_events]) : weights extrapolated at reference point, for all events
    xsec_refPoint (float) : xsec value at the reference point
    weights_thetas (ndarray of shape [n_events, n_points]) : weights extrapolated at reference point, for all events and all EFT points
    xsecs_thetas (ndarray of shape [n_points]) : xsec values at all EFT points

    Returns:
    JLRs (ndarray of shape [n_events, n_points]) : JLR values computed for all events, for all EFT points
    """

    dw_thetas = np.divide(weights_thetas, xsecs_thetas)
    dw_ref = np.divide(weights_refPoint, xsec_refPoint)

    JLRs = np.divide(dw_thetas, dw_ref[:,np.newaxis])

    return JLRs

# //--------------------------------------------
# //--------------------------------------------

def Compute_Score_Component(weights_thetas, xsecs_thetas, gradWeights_thetas, gradXsecs_thetas):
    """
    Compute the 'score' quantity for a given EFT hypothesis, defined as grad_theta(log(JLR)) = (grad_theta(w|theta0) / w_theta0) - (grad_theta(xsec|theta0) / xsec_theta0)
    The score is a vectorial quantity, with as many elements as there are EFT operators defining 'theta0' (can differentiate theta0 w.r.t. each operator). This function returns the score component for 1 specific operator at once.

    Parameters:
    weights_thetas (ndarray of shape [n_events, n_points]) : weights extrapolated at all EFT points, for all events and all EFT points
    xsecs_thetas (ndarray of shape [n_points]) : xsec values extrapolated at all EFT points
    gradWeights_thetas (ndarray of shape [n_events, n_points]) : weight gradients extrapolated at all EFT points, for all events
    gradXsecs_thetas (ndarray of shape [n_points]) : xsec gradients extrapolated at all EFT points

    Returns:
    score_component (ndarray of shape [n_events, n_points]) : score component corresponding to differentiation w.r.t. considered operator, for all events and all EFT points
    """

    dw_thetas = np.divide(gradWeights_thetas, weights_thetas) #A = grad(w)/w
    dxs_thetas = np.divide(gradXsecs_thetas, xsecs_thetas) #B = grad(xsec)/xsec

    score_component = np.subtract(dw_thetas, dxs_thetas) #A - B

    return score_component

# //--------------------------------------------
# //--------------------------------------------

 #####  ###### #####  #    #  ####
 #    # #      #    # #    # #    #
 #    # #####  #####  #    # #
 #    # #      #    # #    # #  ###
 #    # #      #    # #    # #    #
 #####  ###### #####   ####   ####

def CrossCheck_FewEvents(opts, n_components, components, operatorNames, fit_coeffs, EFT_weights, EFT_weightIDs, list_weights, list_SMweights):
    """
    Debug printouts for first events. Check fit coefficients, compare event weights/xsecs from MG benchmark weights and from extrapolation procedure, etc.
    """

    print(colors.fg.orange, '=========== DEBUG =========\n\n', colors.reset)

    print('\nComponents --> \n', components)
    print('\nFit coeffs [:5] --> \n', fit_coeffs[:5])

    print("-----------")
    print('\nSM MG reweights [:5] --> \n', list_SMweights[:5])
    sm_weights = Extrapolate_SM_Weights(n_components, components, fit_coeffs, opts["listOperatorsParam"])
    print('SM extrapol. reweights [:5] --> \n', sm_weights[:5])

    print("-----------")
    sm_xsec = Extrapolate_SM_Xsec(opts, sm_weights)
    print('SM xsec from extrapolated weights: ', sm_xsec)
    print('SM xsec from MG weights: ', Extrapolate_EFTxsecs_fromWeights(list_SMweights))

    # print("-----------")
    # baseline_WCs = Translate_EFTpointID_to_WCvalues(operatorNames, "rwgt_ctZ_5_ctW_5_cpQM_5_cpQ3_5_cpt_5")
    # effWC_components_baseline = Get_EffectiveWC_eachComponent(n_components, components, baseline_WCs) #Get corresponding matrix
    # weights_baseline = Extrapolate_EFTweights(effWC_components_baseline, fit_coeffs) #Get corresponding event weights
    # weights_baseline = np.squeeze(weights_baseline) #2D -> 1D (single point)
    # print('\nBaseline weight from MG [:5] --> ', list_weights[:5])
    # print('Baseline weight from extrapol [:5] --> ', weights_baseline[:5])

    print("-----------")
    idx_refPoint = EFT_weightIDs.shape[1]-1 #Last rwgt point
    print('RefPoint ID : ', EFT_weightIDs[0,idx_refPoint])
    # lasPoint_WCs = Translate_EFTpointID_to_WCvalues(operatorNames, "rwgt_ctZ_0_ctW_0_cpQM_0_cpQ3_0_cpt_5")
    lasPoint_WCs = Translate_EFTpointID_to_WCvalues(operatorNames, "rwgt_ctZ_4.07_ctW_4.43_cpQM_3.7_cpQ3_-1.16_cpt_-1.82 ")
    effWC_components_refPoint = Get_EffectiveWC_eachComponent(n_components, components, lasPoint_WCs) #Get corresponding matrix
    weights_refPoint = Extrapolate_EFTweights(effWC_components_refPoint, fit_coeffs) #Get corresponding event weights
    weights_refPoint = np.squeeze(weights_refPoint) #2D -> 1D (single point)
    print('\nLast point weight from MG [:5] --> ', EFT_weights[:5,idx_refPoint])
    print('Ref point weight from extrapol [:5] --> ', weights_refPoint[:5])
    print('Ref point xsec from extrapol [:5] --> ', Extrapolate_EFTxsecs_fromWeights(weights_refPoint))

    print("-----------")
    print(colors.fg.red, 'Exiting CrossCheck_FewEvents(), which is for debugging only. Stopping here', colors.reset); exit(1)
    return

# //--------------------------------------------
# //--------------------------------------------

######## ##     ## ######## ######## ##    ## ########
##        ##   ##     ##    ##       ###   ## ##     ##
##         ## ##      ##    ##       ####  ## ##     ##
######      ###       ##    ######   ## ## ## ##     ##
##         ## ##      ##    ##       ##  #### ##     ##
##        ##   ##     ##    ##       ##   ### ##     ##
######## ##     ##    ##    ######## ##    ## ########

########     ###    ########    ###     ######  ######## ########
##     ##   ## ##      ##      ## ##   ##    ## ##          ##
##     ##  ##   ##     ##     ##   ##  ##       ##          ##
##     ## ##     ##    ##    ##     ##  ######  ######      ##
##     ## #########    ##    #########       ## ##          ##
##     ## ##     ##    ##    ##     ## ##    ## ##          ##
########  ##     ##    ##    ##     ##  ######  ########    ##

#Main EFT function, called in GetData.py. Extend the data (x, weights, ...) at all the EFT points theta at which the NN will get trained.
def Extend_Augment_Dataset(opts, list_labels, list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses, list_SMweights_allClasses):
    """
    Extend the original data so that the NN can be parameterized on WCs. Events need to be reweighted at many different EFT points for the NN to learn how to interpolate between them.
    First, define the points 'thetas' on which the NN will get trained. Then, for the same input features x, will duplicate events N times at these points, with corresponding reweights.
    Can also 'augment' the data, i.e. compute the JLR and score quantities on which the NN could then regress (or use as inputs).
    Samples are first 'unweighted', i.e. events are drawn from samples with probability corresponding to their relative weights, and then al attributed a weight of 1.

    Parameters:
    listOperatorsParam : list of operators considered in NN training
    list_labels, list_x_allClasses, list_weights_allClasses, list_EFTweights_allClasses, list_EFTweightIDs_allClasses, list_EFT_FitCoeffs_allClasses, list_SMweights_allClasses : obtained from previous functions
    maxEvents : number of events to sample for each EFT point. Events will be drawn randomly, with replacement, from the entire class sample with the proper probability.
    nPointsPerOperator, minWC, maxWC : define the interval [min, max, step] in which points will be drawn uniformly, separately for each operator ; these are the EFT points on which the NN will get trained

    Returns:
    Extended/augmented lists
    """

    parameterizedNN = opts["parameterizedNN"]

    need_jlr = (opts["strategy"] in ["ROLR", "RASCAL"])
    need_score = (opts["strategy"] in ["RASCAL"])

# //--------------------------------------------
# Sanity checks

    if parameterizedNN == False and opts["strategy"] != "CARL_singlePoint": #Return empty lists #Exception: for strategy 'CARL_singlePoint', still need to extend the lists
        list_thetas_allClasses = []; list_targetClass_allClasses = [];
        return list_x_allClasses, list_weights_allClasses, [], [], [], []

    if DEBUG_>1: #Debug printout
        for iclass in range(len(list_x_allClasses)):
            list_x_allClasses[iclass]=list_x_allClasses[iclass][:5]; list_weights_allClasses[iclass]=list_weights_allClasses[iclass][:5]; list_EFTweights_allClasses[iclass]=list_EFTweights_allClasses[iclass][:5]; list_EFTweightIDs_allClasses[iclass]=list_EFTweightIDs_allClasses[iclass][:5]; list_EFT_FitCoeffs_allClasses[iclass]=list_EFT_FitCoeffs_allClasses[iclass][:5]

# //--------------------------------------------
# Define hypotheses to train on, loop on process classes

    #Extended lists to return
    extendedList_x_allClasses = []
    extendedList_weights_allClasses = []
    extendedList_targetClass_allClasses = []
    extendedList_thetas_allClasses = []
    extendedList_jointLR_allClasses = []
    extendedList_score_allClasses_allOperators = [ [] for iclass in range(len(list_labels))] #List of list because score is a vector --> will store values for each class, for each score component

    rng = np.random.default_rng() #Init random generator

    #-- Loop on process classes (each expected to contain a single EFT process)
    for iclass in range(len(list_labels)):

        #-- Sanity checks
        operatorNames, _, _ = Parse_EFTpoint_IDs(list_EFTweightIDs_allClasses[iclass][0,0]) #Asumme that benchmark points are identical for all events in the sample
        operatorNames = operatorNames[0] #Single-element list -> array
        for op1 in opts["listOperatorsParam"]:
            if "PrivMC" in list_labels[iclass] and op1 not in operatorNames and op1.lower() not in operatorNames: print(colors.fg.red, 'Error : parameterized operator ', op1,' not found in class', list_labels[iclass], colors.reset, ' (Operators found : ', operatorNames, ')'); exit(1)

        nEventsPerPoint_class = opts["maxEvents"]
        if nEventsPerPoint_class > len(list_x_allClasses[iclass]):
            print(colors.fg.orange, 'Warning : requiring more events per EFT point(', nEventsPerPoint_class, ') than available in this class (', len(list_x_allClasses[iclass]), ') ! Setting param. \'maxEvents\' accordingly...', colors.reset)
            nEventsPerPoint_class = len(list_x_allClasses[iclass])
        elif nEventsPerPoint_class ==-1: nEventsPerPoint_class = len(list_x_allClasses[iclass]) #Use entire dataset at each point

        #-- Determine the components necessary to parameterize the event weights
        n_components, components = Find_Components(operatorNames) #Determine the components required to parameterize the event weight #NB: assumes that they are identical for all events in process class
        # print(components)

        #-- Define the hypotheses theta on which the NN will get trained. Draw values uniformly in given interval for each operator, translate into array. Also define the corresponding target to train on (<-> the operator which is activated at a given point) #NB: points theta are defined according to user-defined list 'listOperatorsParam', but must be shaped according to total nof operators present in samples
        #NB: could do that before class loop (also Find_Components), since expect all classes to share exact same parameterization. But just in case different classes have different parameterizations, and we are only considering a common subset of operators, define everything separately for each class.
        thetas, targetClass = Get_ThetaParameters(opts, operatorNames)
        # print(thetas); print(targetClass)

        #-- For given class, append arrays for all points thetas
        list_x_allThetas_class = []
        list_weights_allThetas_class = []
        list_thetas_allThetas_class = []
        list_targetClass_allThetas_class = []
        list_jointLR_allThetas_class = []
        list2D_score_allOperators_allThetas_class = [ [] for iop in range(len(opts['listOperatorsParam']))] #List of list because score is a vector (first index corresponds to score component <-> operator) #Initialize an empty list for each operator

# //--------------------------------------------
# Start by computing necessary quantities (reweights, augmented data, etc.) for entire sample

        #-- Extrapolate the event weights at the reference point, for all events
        if opts["refPoint"] == "SM":
            weights_refPoint = list_SMweights_allClasses[iclass] #Choose SM as reference point
        else:
            refPointWCs = Translate_EFTpointID_to_WCvalues(operatorNames, opts["refPoint"]) #Encode reference point name into WC values
            effWC_components_refPoint = Get_EffectiveWC_eachComponent(n_components, components, refPointWCs) #Get corresponding matrix
            weights_refPoint = Extrapolate_EFTweights(effWC_components_refPoint, list_EFT_FitCoeffs_allClasses[iclass]) #Get corresponding event weights
            weights_refPoint = np.squeeze(weights_refPoint) #2D -> 1D (single point)

        if parameterizedNN == True:

            #-- Get the 'effective WC' values scaling each fit component
            effWC_components_thetas_class = Get_EffectiveWC_eachComponent(n_components, components, thetas) #Determine the 'effective WC' values associated with each component, for each benchmark point
            # print(effWC_components_thetas_class)

            #-- Extrapolate the event weights at the new points thetas, for all events
            newWeights = Extrapolate_EFTweights(effWC_components_thetas_class, list_EFT_FitCoeffs_allClasses[iclass])

        #-- For non-parameterized NN, only 1 EFT point to separate from SM, build lists differently
        else:
            newWeights = list_SMweights_allClasses[iclass] #Choose SM as reference point
            list_weights_allThetas_class.append(np.concatenate( (newWeights, weights_refPoint)) )
            list_x_allThetas_class.append(np.concatenate( (list_x_allClasses[iclass], list_x_allClasses[iclass])) )
            list_targetClass_allThetas_class.append(np.concatenate( (np.ones(len(list_x_allClasses[iclass])), np.zeros(len(list_x_allClasses[iclass])) ) ) ) #EFT ref point=0, SM=1

        #-- Extrapolate the cross section at the reference point and the new points thetas
        if need_jlr or need_score:
            newXsecs = Extrapolate_EFTxsecs_fromWeights(newWeights)
            xsecs_refPoint = Extrapolate_EFTxsecs_fromWeights(weights_refPoint)

            #-- Compute joint likelihood ratio at the new points thetas, for all events
            if need_jlr: jointLR_class = Compute_JointLR(weights_refPoint, xsecs_refPoint, newWeights, newXsecs)

            list_gradNewWeights_operators = []; list_gradNewXsecs_operators = []; list_score_allOperators = []
            if need_score:

                #-- Get the gradients of the 'effective WC' values scaling each fit component (3D array with shape (n_operators, n_points, n_components) )
                gradEffWC_components_operators_thetas_class = Get_Gradients_EffectiveWC_eachComponent(n_components, components, thetas)

                #-- Differentiating w.r.t. 1 operator at once, get the corresponding gradients of weights and xsecs at all points for all events, and use them to get the corresponding score component
                for iop in range(len(opts['listOperatorsParam']) ):
                    list_gradNewWeights_operators.append(Extrapolate_EFTweights(gradEffWC_components_operators_thetas_class[iop,...], list_EFT_FitCoeffs_allClasses[iclass]))
                    list_gradNewXsecs_operators.append(Extrapolate_EFTxsecs_fromWeights(list_gradNewWeights_operators[iop]))
                    list_score_allOperators.append(Compute_Score_Component(newWeights, newXsecs, list_gradNewWeights_operators[iop], list_gradNewXsecs_operators[iop]) )

        #-- Debugging
        # CrossCheck_FewEvents(opts, n_components, components, operatorNames, list_EFT_FitCoeffs_allClasses[iclass], list_EFTweights_allClasses[iclass], list_EFTweightIDs_allClasses[iclass], list_weights_allClasses[iclass], list_SMweights_allClasses[iclass])

# //--------------------------------------------
# Then, for each new point theta, build the corresponding dataset (features, weights, etc.), and append it to total dataset
# /!\ Reminder: both for classifying H0 vs H1 and LR regression r(H0,H1), need to provide events *both at theta (H1) and at a reference point (H0)* ! The reference point is arbitrary, it can be fixed or not, and is a hyperparameter. For convenience, we set it to the SM point.

        for itheta in range(len(thetas)):

            #-- Sample unweighting: translate the event weights (at point theta and for reference point) into probabilities to build the sample's PDF
            #NB: ignoring negative weights for now
            probas_point = newWeights[:,itheta]
            probas_refPoint = weights_refPoint
            probas_refPoint[probas_refPoint < 0] = 0; probas_point[probas_point < 0] = 0 #Ignore neg weights
            probas_refPoint /= probas_refPoint.sum(); probas_point /= probas_point.sum() #Normalize to 1

            #-- Draw N events from sample according to PDF (either at theta or ref point), with replacement
            #NB: events with large weights may get drawn many times, degrading the NN performance
            indices_refPoint = rng.choice(len(list_x_allClasses[iclass]), size=nEventsPerPoint_class, p=probas_refPoint)
            indices_point = rng.choice(len(list_x_allClasses[iclass]), size=nEventsPerPoint_class, p=probas_point)
            # print(indices_refPoint); print(indices_point)

            #-- Get the features of selected events
            #NB: at this point, input features are still stored into 1D structured arrays. Reshaped in 2D later
            x_class_refPoint = list_x_allClasses[iclass][indices_refPoint]
            x_class_point = list_x_allClasses[iclass][indices_point]

            #-- Both events at theta and ref point will share the same 'theta' given as input to the NN, i.e. they are both fed to the NN for the same value of theta. Hence, duplicate theta (<-> the WC values defining the current point) for all selected events
            thetas_class_point = np.tile(thetas[itheta,:], (nEventsPerPoint_class*2, 1) ) #Once for events at theta, once for events at ref point
            #Each original event gets repeated at all the different points theta

            #-- For classification, only 1 operator is 'activated' (non-zero) at once. Keep track of which operator it is, for event labelling
            targetClass_class_point = np.tile(targetClass[itheta], (nEventsPerPoint_class, 1) )
            if opts["nofOutputNodes"] == 1 or opts["regress"] == True: targetClass_class_refPoint = np.zeros((targetClass_class_point.shape))
            else: targetClass_class_refPoint = np.zeros((targetClass_class_point.shape)); targetClass_class_refPoint[:,0] = 1 #Ref point <-> first row in multiclass
            targetClass_class_point = np.squeeze(targetClass_class_point); targetClass_class_refPoint = np.squeeze(targetClass_class_refPoint)

            #-- After unweighting, all events are attributed a weight of 1.
            weights_class_point = np.ones(nEventsPerPoint_class*2)
            # print(thetas_class); print(weights_class)

# //--------------------------------------------
# Concatenate data for selected events at ref point and at theta

            list_x_allThetas_class.append(np.concatenate((x_class_refPoint, x_class_point)) )
            list_weights_allThetas_class.append(weights_class_point)
            list_targetClass_allThetas_class.append(np.concatenate((targetClass_class_refPoint, targetClass_class_point)) )
            list_thetas_allThetas_class.append(thetas_class_point)
            if need_jlr: list_jointLR_allThetas_class.append(np.concatenate((jointLR_class[indices_refPoint,itheta], jointLR_class[indices_point,itheta])))
            if need_score:
                for iop in range(len(opts['listOperatorsParam'])): list2D_score_allOperators_allThetas_class[iop].append(np.concatenate((list_score_allOperators[iop][indices_refPoint,itheta], list_score_allOperators[iop][indices_point,itheta])) )

# //--------------------------------------------
# Append arrays for all classes
        extendedList_x_allClasses.append(np.concatenate(list_x_allThetas_class))
        extendedList_weights_allClasses.append(np.concatenate(list_weights_allThetas_class))
        extendedList_targetClass_allClasses.append(np.concatenate(list_targetClass_allThetas_class))
        if parameterizedNN == True:
            extendedList_thetas_allClasses.append(np.concatenate(list_thetas_allThetas_class))
            if need_jlr: extendedList_jointLR_allClasses.append(np.concatenate(list_jointLR_allThetas_class))
            if need_score:
                for iop in range(len(opts['listOperatorsParam'])): extendedList_score_allClasses_allOperators[iclass].append(np.concatenate(list2D_score_allOperators_allThetas_class[iop][:])) #Concatenate different points thetas together, but retain correct operator ordering (for each process class)

    return extendedList_x_allClasses, extendedList_weights_allClasses, extendedList_targetClass_allClasses, extendedList_thetas_allClasses, extendedList_jointLR_allClasses, extendedList_score_allClasses_allOperators

# //--------------------------------------------
# //--------------------------------------------










# //--------------------------------------------
