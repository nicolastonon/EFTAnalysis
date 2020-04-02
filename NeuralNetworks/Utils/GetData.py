# Read ROOT files, shape/transform the data (x : input features, y : labels), compute event reweights, ...
#NB : if want to add a validation test, could simply split the test sample again with train_test_split...

# //--------------------------------------------
#Filtering out manually some unimportant warnings
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
warnings.filterwarnings("ignore", message="Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated")
# --------------------------------------------
import ROOT
from ROOT import TMVA, TFile, TTree, TCut, gROOT, TH1, TH1F
import numpy as np
from root_numpy import root2array, tree2array, array2root
import pandas as pd
import tensorflow
import keras
from pathlib import Path
from sklearn.utils import class_weight
from sklearn.feature_selection import RFE, SelectKBest, chi2
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, Normalizer, StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold
from tensorflow.keras import utils

from Utils.ColoredPrintout import colors
from Utils.Helper import get_normalization_iqr

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
def Get_Data_For_DNN_Training(weight_dir, lumi_years, ntuples_dir, processClasses_list, labels_list, var_list, cuts, nof_output_nodes, maxEvents_perClass, splitTrainEventFrac, nEventsTot_train, nEventsTot_test, lumiName, transfType='quantile'):

    #Get data from TFiles
    list_x_allClasses, list_weights_allClasses = Read_Store_Data(lumi_years, ntuples_dir, processClasses_list, labels_list, var_list, cuts)

    #Shape the data arrays properly #Also read arguments which are modified by the function
    x, y, list_weights_allClasses = Shape_Data(list_x_allClasses, list_weights_allClasses, maxEvents_perClass, nof_output_nodes)

    #Before we randomize the events, store the input values of the very first events (which belong to first process) --> Can use these known events for later validation/comparison
    x_control_firstNEvents = x[0:10,:]

    #Compute event weights considering only abs(weights) => duplicate weight arrays to hold absolute values
    list_weights_allClasses_abs = []
    for weights_class in list_weights_allClasses:
        list_weights_allClasses_abs.append(np.absolute(weights_class))

    #Get event-per-event reweights (for training phase)
    LearningWeights_allClasses, PhysicalWeights_allClasses = Get_Events_Weights(processClasses_list, labels_list, list_weights_allClasses_abs, list_weights_allClasses)

    #-- Define training & testing subsamples (takes care of splitting + shuffling)
    #-- http://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html
    #-- Default args : shuffle=True <-> shuffle events ; random_state=None <-> random seed ; could also use stratify=y so that the final splitting respects the class proportions of the array y, if desired (else : random)
    if (nEventsTot_train is not -1) and (nEventsTot_test is not -1): #Specify nof train/test events
        x_train, x_test, y_train, y_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test = train_test_split(x, y, PhysicalWeights_allClasses, LearningWeights_allClasses, train_size=nEventsTot_train, test_size=nEventsTot_test, shuffle=True)
    else: #Specify train/test relative proportions
        x_train, x_test, y_train, y_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test = train_test_split(x, y, PhysicalWeights_allClasses, LearningWeights_allClasses, test_size=1-splitTrainEventFrac, shuffle=True) #X% train, Y% test


    #Get rescaling parameters for each input feature, given to first DNN layer to normalize features -- derived from train data alone
    xTrainRescaled, shifts, scales = Transform_Inputs(weight_dir, x_train, var_list, lumiName, transfType)

    print(colors.fg.lightblue, "===========", colors.reset)
    print(colors.fg.lightblue, "-- Will use " + str(x_train.shape[0]) + " training events !", colors.reset)
    print(colors.fg.lightblue, "-- Will use " + str(x_test.shape[0]) + " testing events !", colors.reset)
    print(colors.fg.lightblue, "===========\n", colors.reset)

    return x_train, y_train, x_test, y_test, PhysicalWeights_train, PhysicalWeights_test, LearningWeights_train, LearningWeights_test, x, y, PhysicalWeights_allClasses, LearningWeights_allClasses, shifts, scales, x_control_firstNEvents, xTrainRescaled
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
def Read_Store_Data(lumi_years, ntuples_dir, processClasses_list, labels_list, var_list, cuts):

    testdirpath = Path(ntuples_dir)
    if not testdirpath.is_dir():
        print('Ntuple dir. '+ntuples_dir+' not found ! Abort !')
        exit(1)

    list_x_allClasses = [] #List of x-arrays storing the values of input variables for all events, for all considered years
    list_weights_allClasses = [] #Idem to store central event weights
    list_EFTreweights_allClasses = [] #Idem to store EFT reweights

    for procClass, label in zip(processClasses_list, labels_list): #NB : can not use keyword 'class'

        print('\n', colors.fg.purple, colors.underline, '* Class :', label, colors.reset)

        list_x_proc = []
        list_weights_proc = []
        list_EFTreweights_proc = []
        for process in procClass:

            print('\n', colors.fg.pink, '* Process :', colors.reset, process)

            for year in lumi_years:

                filepath = ntuples_dir + year + '/' + process + '.root'
                testfilepath = Path(filepath)
                if not testfilepath.is_file():
                    print('File '+filepath+' not found ! Abort !')
                    exit(1)

                # print(colors.fg.lightgrey, '* Opening file:', colors.reset, ' ', filepath)
                file = TFile.Open(filepath)
                tree = file.Get('result')
                print(colors.fg.lightgrey, '* Opened file:', colors.reset, filepath, '(', tree2array(tree, branches="eventWeight", selection=cuts).shape[0], 'entries )') #Dummy variable, just to read the nof entries

                list_x_proc.append(tree2array(tree, branches=var_list, selection=cuts))
                list_weights_proc.append(tree2array(tree, branches="eventWeight*eventMCFactor", selection=cuts))
                # list_weights_proc.append(tree2array(tree, branches="eventWeight", selection=cuts))
                # if procname contains privprod inly
                # list_EFTreweights_proc.append(tree2array(tree, branches="mc_EFTweights", selection=cuts))

        list_x_allClasses.append(np.concatenate(list_x_proc))
        list_weights_allClasses.append(np.concatenate(list_weights_proc))
        # list_EFTreweights_allClasses.append(np.concatenate(list_EFTreweights_proc))

    print('\n\n')

    return list_x_allClasses, list_weights_allClasses









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

#Shape the data into arrays (per process, lumiYear, etc.)
def Shape_Data(list_x_allClasses, list_weights_allClasses, maxEvents_perClass, nof_output_nodes):

    #Reshape as normal arrays (root_numpy uses different format) : 1 column per variable, 1 line per event
    list_x_arrays_allClasses = []
    for x_class in list_x_allClasses:
        list_x_arrays_allClasses.append(x_class.view(np.float32).reshape(x_class.shape + (-1,)) )

    #--- Get nof entries for each class
    list_nrows_class = []
    total_rows = 0
    for i in range(len(list_x_arrays_allClasses)):
        list_nrows_class.append(list_x_arrays_allClasses[i].shape[0])
        total_rows+= list_x_arrays_allClasses[i].shape[0]
    # print(total_rows)

#--- Max nof events for train.test phases

    if maxEvents_perClass is not -1:

        #Skim each process class (keep only maxEvents_perClass events) ; skim all relevant arrays coherently
        for i in range(len(list_x_arrays_allClasses)):
            if list_nrows_class[i] > maxEvents_perClass:
                list_nrows_class[i] = maxEvents_perClass
                list_x_arrays_allClasses[i] = list_x_arrays_allClasses[i][0:maxEvents_perClass]
                list_x_allClasses[i] = list_x_allClasses[i][0:maxEvents_perClass] #Don't forget to also skim this list, used again later
                list_weights_allClasses[i] = list_weights_allClasses[i][0:maxEvents_perClass]

    x = np.concatenate(list_x_arrays_allClasses, 0)

    #Create array of labels (1 row per event, 1 column per class)

    if nof_output_nodes == 1: #binary, single column => sig 1, bkg 0
        y_integer_sig = np.ones(list_nrows_class[0]) #'1' = signal
        y_integer_bkg = np.zeros(list_nrows_class[1]) #'0' = bkg
        y = np.concatenate((y_integer_sig, y_integer_bkg), axis=0)
    else: #multiclass, n columns => sig 1, bkg 0
        list_y_integer_allClasses = []
        for i in range(len(list_nrows_class)): #Concatenate subsequent classes
            list_y_integer_allClasses.append(np.full((list_nrows_class[i]), i) )
        y_integer = np.concatenate(list_y_integer_allClasses, 0)
        y = utils.to_categorical(y_integer, num_classes=nof_output_nodes) #Converts a class vector (integers) to binary class matrix #One-hot encode the integers

    return x, y, list_weights_allClasses











# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
########  ######## ##      ## ######## ####  ######   ##     ## ########  ######
##     ## ##       ##  ##  ## ##        ##  ##    ##  ##     ##    ##    ##    ##
##     ## ##       ##  ##  ## ##        ##  ##        ##     ##    ##    ##
########  ######   ##  ##  ## ######    ##  ##   #### #########    ##     ######
##   ##   ##       ##  ##  ## ##        ##  ##    ##  ##     ##    ##          ##
##    ##  ##       ##  ##  ## ##        ##  ##    ##  ##     ##    ##    ##    ##
##     ## ########  ###  ###  ######## ####  ######   ##     ##    ##     ######
# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------

#Compute and apply weights to training dataset to balance the training
def Get_Events_Weights(processClasses_list, labels_list, list_weights_allClasses_abs, list_weights_allClasses):

    #Compute 'yields' (from *absolute* weights) to reweight classes
    list_yields_abs_allClasses = []
    yield_abs_total = 0
    for i in range(len(processClasses_list)):
        list_yields_abs_allClasses.append(list_weights_allClasses_abs[i].sum())
        yield_abs_total+= list_weights_allClasses_abs[i].sum()

    #Compute scale factors to rescale each class to 'yield_abs_total'
    list_SFs_allClasses = []
    for i in range(len(processClasses_list)):
        list_SFs_allClasses.append(yield_abs_total / list_yields_abs_allClasses[i])
        print('Class', labels_list[i], ' / Scale factor = ', round(list_SFs_allClasses[i], 2), '===> Rescaled yield :', round(list_yields_abs_allClasses[i]*list_SFs_allClasses[i], 1) )
    print('\n')

    #Get array of reweighted 'training' weights, i.e. used for training only and which are not physical
    list_LearningWeights_allClasses = []
    for i in range(len(processClasses_list)):
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
def Transform_Inputs(weight_dir, x_train, var_list, lumiName, transfType='quantile'):

    np.set_printoptions(precision=3)
    # print('Before transformation :', x_train[0:5,:])

    if transfType not in ['quantile', 'range', 'gauss']:
        print('\n', colors.fg.red, 'Warning : transfType', colors.reset, transfType, colors.fg.red, 'not known ! Using default [quantile]', colors.reset)
        transfType = 'quantile'

    #--- QUANTILE RESCALING
    # a = median ; b = scale
    if transfType == 'quantile':
        shift_, scale_ = get_normalization_iqr(x_train, 0.95)
        xTrainRescaled = (x_train - shift_) / scale_ #Apply transformation on train events -- for control

    #--- RANGE SCALING
    # a = min ; b = scale
    elif transfType == 'range':
        scaler = MinMaxScaler(feature_range=(-1, 1)).fit(x_train) #Conpute macro parameters
        shift_ = scaler.min_
        scale_ = scaler.scale_
        x_train = scaler.fit(x_train) #Get params
        xTrainRescaled = scaler.transform(x_train) #Apply transformation on train events -- for control

    #--- RESCALE TO UNIT GAUSSIAN -- https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
    # a = mean ; b = stddev
    if transfType == 'quantile':
        scaler = StandardScaler().fit(x_train) #Get params
        shift_ = scaler.mean_
        scale_ = scaler.scale_ # = np.sqrt(var_)
        xTrainRescaled = scaler.transform(x_train) #Apply transformation on train events -- for control

    text_file = open(weight_dir + "DNN_infos.txt", "w")

    #Dump shift_ and scale_ params into txtfile
    for ivar in range(len(var_list)):
        text_file.write(var_list[ivar]); text_file.write(' ')
        text_file.write(str(shift_[ivar])); text_file.write(' ')
        text_file.write(str(scale_[ivar])); text_file.write('\n')

    text_file.close()
    print(colors.fg.lightgrey, '\n===> Saved DNN infos (input/output nodes names, rescaling values, etc.) in : ', weight_dir + "DNN_infos.txt \n", colors.reset)

    # print('After transformation :', x_train[0:5,:])

    return xTrainRescaled, shift_, scale_
