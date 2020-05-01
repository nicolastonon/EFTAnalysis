# OUTDATED code to perform NN hyperparameters optimization, using the Hyperas library

from __future__ import print_function
from hyperopt import Trials, STATUS_OK, tpe
from hyperas import optim
from hyperas.distributions import choice, uniform, conditional

import ROOT
from ROOT import TMVA, TFile, TTree, TCut, gROOT
import numpy as np
from root_numpy import root2array, tree2array
from keras.utils import np_utils
from sklearn.model_selection import train_test_split

import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.models import model_from_json
from keras.optimizers import SGD, Adam
from keras.regularizers import l2
from sklearn.preprocessing import MinMaxScaler, Normalizer

# //--------------------------------------------
#Filtering out manually some unimportant warnings
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
warnings.filterwarnings("ignore", message="Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated")

# //--------------------------------------------
# //--------------------------------------------





########     ###    ########    ###
##     ##   ## ##      ##      ## ##
##     ##  ##   ##     ##     ##   ##
##     ## ##     ##    ##    ##     ##
##     ## #########    ##    #########
##     ## ##     ##    ##    ##     ##
########  ##     ##    ##    ##     ##



# def data(ttH, ttV, is_trigger==1, 2l, CR_ttW):
def data():

    #Write helper inner function here. Because hyperas functions don't take args, etc.
    def Get_Boolean_Categ(nLep, region):
        if nLep == "3l":
            if region == "SR":
                return "is_tHq_3l_SR"
            elif region == "CR_ttZ":
                return "is_tHq_ttZctrl"
        elif nLep == "2l":
            if region == "SR":
                return "is_tHq_2lSS_SR"
            elif region == "CR_ttW":
                return "is_tHq_ttWctrl"

        return ""

    # Analysis options
    # //--------------------------------------------
    signal = "ttH" #Choose "tHq" or "ttH"
    bkg_type = "ttV" # Choose "ttV" or "ttbar"
    nLep = "2l" # Choose "3l" or "2l"
    region = "CR_ttW" #SR, CR_ttW, CR_ttZ
    use_2016_ntuples = False
    include_MEM_variables = False
    cuts = "is_trigger==1" #"1" <-> no cut

    dir_samples = "../input_ntuples/" + nLep + "/" + region + "/"

    # Get files
    file_ttH = TFile.Open(dir_samples + signal + '.root')

    print()
    print("===========")
    print("--- Signal = " + signal)

    file_bkg1 = TFile.Open(dir_samples + "ttZ.root")
    file_bkg2 = TFile.Open(dir_samples + "ttW.root")


    sel = cuts + "==1 && "+ Get_Boolean_Categ(nLep, region) + "==1"
    # sel = "is_tHq_3l_Training==1" #To use "Training" categ. events

    #Get trees
    tree_sig = file_ttH.Get('Tree')
    tree_bkg1 = file_bkg1.Get('Tree')
    tree_bkg2 = file_bkg2.Get('Tree')



  ####  #    #   ##   #####  ######    #####    ##   #####   ##
 #      #    #  #  #  #    # #         #    #  #  #    #    #  #
  ####  ###### #    # #    # #####     #    # #    #   #   #    #
      # #    # ###### #####  #         #    # ######   #   ######
 #    # #    # #    # #      #         #    # #    #   #   #    #
  ####  #    # #    # #      ######    #####  #    #   #   #    #


    # Define list of input variables
    # //--------------------------------------------
    var_list = []
    if use_2016_ntuples == True:
        var_list.append('nJet25')
        var_list.append('maxEtaJet25')
        var_list.append('lepCharge')
        var_list.append('nJetEta1')
        var_list.append('dEtaFwdJetBJet')
        var_list.append('dEtaFwdJet2BJet')
        var_list.append('dEtaFwdJetClosestLep')
        var_list.append('dPhiHighestPtSSPair')
        var_list.append('minDRll')
    else:
        if bkg_type == "ttV":
            var_list.append("lep1_conePt")
            var_list.append("lep2_conePt")
            var_list.append("mindr_lep1_jet")
            var_list.append("mindr_lep2_jet")
            var_list.append("mT_lep1")
            var_list.append("mT_lep2")
            var_list.append("max_lep_eta")
        else:
            var_list.append("mindr_lep1_jet")
            var_list.append("mindr_lep2_jet")
            var_list.append("mT_lep1")
            var_list.append("mT_lep2")
            var_list.append("max_lep_eta")


    if include_MEM_variables == True:
        var_list.append('log_ttZ_weight')
        var_list.append('log_ttZ_weight_kinmaxint')
        var_list.append('log_ttW_weight')
        var_list.append('log_ttW_weight_kinmaxint')
        var_list.append('log_tHq_weight')
        var_list.append('log_tHq_weight_kinmaxint')
    # //--------------------------------------------


    #Use 'root_numpy' to get numpy arrays from root files
    x_sig = tree2array(tree_sig, branches=var_list, selection=sel)
    x_bkg1 = tree2array(tree_bkg1, branches=var_list, selection=sel)
    x_bkg2 = tree2array(tree_bkg2, branches=var_list, selection=sel)
    # print(x_sig.shape)

    # exit(1)

    #Reshape as normal arrays (root_numpy uses different format) : 1 column per variable, 1 line per event
    x_sig_array = x_sig.view(np.float32).reshape(x_sig.shape + (-1,))
    x_bkg1_array = x_bkg1.view(np.float32).reshape(x_bkg1.shape + (-1,))
    x_bkg2_array = x_bkg2.view(np.float32).reshape(x_bkg2.shape + (-1,))
    # print(x_sig_array.shape)

    #Concatenate all 3 processes into 1 single array containing both ttH+bkgs
    x = np.concatenate((x_sig_array, x_bkg1_array, x_bkg2_array), 0)
    # print(x.shape)

    #Create array of labels, 1 column, 1 line per event
    #'1' = ttH, '0' = bdf.
    #NB : for multiclass classification, should add more labels !
    y_integer = np.ones((x.shape[0]), dtype = int)
    y_integer[0:x_sig_array.shape[0]] = 0 #ttH events <-> 0
    y_integer[x_sig.shape[0]:x.shape[0]] = 1 #Bkg events <-> 1
    # print(y_integer.shape)
    # print(y_integer)

    #Convert integers to dummy variables (i.e. one hot encoded) => Use 2 classes (sig, bkg)
    y = np_utils.to_categorical(y_integer, 2) #to_categorical(y, nb_classes=None) #Convert class vector to binary class matrix, for use with categorical_crossentropy.

    # Rescale -- Rescale data inputs between e.g. -1 and 1
    # scaler = MinMaxScaler(feature_range=(0, 1))
    # x = scaler.fit_transform(x)
    # np.set_printoptions(precision=3)
    # print(x[0:5,:])


    X_train, X_val, Y_train, Y_val = train_test_split(x, y, test_size=0.2)

    return X_train, Y_train, X_val, Y_val




##     ##  #######  ########  ######## ##
###   ### ##     ## ##     ## ##       ##
#### #### ##     ## ##     ## ##       ##
## ### ## ##     ## ##     ## ######   ##
##     ## ##     ## ##     ## ##       ##
##     ## ##     ## ##     ## ##       ##
##     ##  #######  ########  ######## ########


def create_model(X_train, Y_train, X_val, Y_val):

    model = Sequential()
    model.add(Dense({{choice([64, 128])}}, input_dim=len(var_list) ) )
    model.add(Activation({{choice(['relu', 'sigmoid'])}}))
    # model.add(Dropout({{uniform(0, 1)}}))
    model.add(Dense({{choice([64, 128])}}))
    model.add(Activation({{choice(['relu', 'sigmoid'])}}))
    # model.add(Dropout({{uniform(0, 1)}}))

    if conditional({{choice(['two', 'three'])}}) == 'three':
        model.add(Dense({{choice([64, 128])}}))
        model.add(Activation({{choice(['relu', 'sigmoid'])}}))
        # model.add(Dropout({{uniform(0, 1)}}))

    model.add(Dense(2))
    model.add(Activation({{choice(['relu', 'sigmoid', 'softmax'])}}))

    #Different choices of optimizers and LR values
    adam = keras.optimizers.Adam(lr={{choice([10**-3, 10**-2, 10**-1])}})
    rmsprop = keras.optimizers.RMSprop(lr={{choice([10**-3, 10**-2, 10**-1])}})
    sgd = keras.optimizers.SGD(lr={{choice([10**-3, 10**-2, 10**-1])}})

    choiceval = {{choice(['adam', 'sgd', 'rmsprop'])}}
    if choiceval == 'adam':
        optim = adam
    elif choiceval == 'rmsprop':
        optim = rmsprop
    else:
        optim = sgd


    model.compile(loss='categorical_crossentropy', metrics=['accuracy'],optimizer=optim)

    model.fit(X_train, Y_train,
              batch_size={{choice([32,64,128])}},
              nb_epoch=20,
              verbose=2,
              validation_data=(X_val, Y_val))
    score, acc = model.evaluate(X_val, Y_val, verbose=0)
    print('Test accuracy:', acc)
    return {'loss': -acc, 'status': STATUS_OK, 'model': model}








# //--------------------------------------------
######## ##     ## ##    ##  ######      ######     ###    ##       ##        ######
##       ##     ## ###   ## ##    ##    ##    ##   ## ##   ##       ##       ##    ##
##       ##     ## ####  ## ##          ##        ##   ##  ##       ##       ##
######   ##     ## ## ## ## ##          ##       ##     ## ##       ##        ######
##       ##     ## ##  #### ##          ##       ######### ##       ##             ##
##       ##     ## ##   ### ##    ##    ##    ## ##     ## ##       ##       ##    ##
##        #######  ##    ##  ######      ######  ##     ## ######## ########  ######
# //--------------------------------------------

X_train, Y_train, X_val, Y_val = data()

best_run, best_model = optim.minimize(model=create_model, data=data, algo=tpe.suggest, max_evals=30, trials=Trials() )

print("Evalutation of best performing model:")
print(best_model.evaluate(X_val, Y_val))
print("Best performing model chosen hyper-parameters:")
print(best_run)
