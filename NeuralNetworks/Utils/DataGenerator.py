# Define a DataGenerator class which can cut an entire dataset into small batches, to be fed sequentially to fit the NN. Rendered necessary for huge training datasets, which do not fit the RAM.
# Basic implementation, following conventions necessary to be given as arg to keras.fit()

import tensorflow
import keras
import numpy as np
from tensorflow.keras.utils import Sequence

from Utils.ColoredPrintout import colors

# //--------------------------------------------
# //--------------------------------------------


########     ###    ########    ###
##     ##   ## ##      ##      ## ##
##     ##  ##   ##     ##     ##   ##
##     ## ##     ##    ##    ##     ##
##     ## #########    ##    #########
##     ## ##     ##    ##    ##     ##
########  ##     ##    ##    ##     ##

 ######   ######## ##    ##
##    ##  ##       ###   ##
##        ##       ####  ##
##   #### ######   ## ## ##
##    ##  ##       ##  ####
##    ##  ##       ##   ### ###
 ######   ######## ##    ## ###

# //--------------------------------------------
# //--------------------------------------------

class DataGenerator(Sequence):
    'Generates data for Keras'

    def __init__(self, x, y, weights, strategy, batch_size=512, returnWeights=True, shuffle=False):
        'Initialization'

        self.x = x
        self.y = y
        self.weights = weights
        self.batch_size = batch_size
        self.strategy = strategy
        self.shuffle = shuffle #Not implemented
        self.returnWeights = returnWeights #True <-> event weights will be used (not needed for parameterized NN, since samples are unweighted)
        # self.on_epoch_end()

    def __len__(self): #Must return an int
        'Denotes the number of batches per epoch'

        return int(np.ceil(len(self.x) / float(self.batch_size)))

    def __getitem__(self, idx):
        'Returns batches of data at each call to DataGenerator'

        #Cut the batches
        batch_x = self.x[idx * self.batch_size:(idx + 1) * self.batch_size]
        batch_y = self.y[idx * self.batch_size:(idx + 1) * self.batch_size]
        batch_weights = self.weights[idx * self.batch_size:(idx + 1) * self.batch_size]

        if self.strategy == "RASCAL": #Need 2 separate targets [r,t]
            if self.returnWeights == False: return batch_x, [batch_y[:,0], batch_y[:,1:]]
            else: batch_x, [batch_y[:,0], batch_y[:,1:]], batch_weights

        else: #Single target
            if self.returnWeights == False: return batch_x, batch_y
            else: batch_x, batch_y, batch_weights

        return batch_x, batch_y, batch_weights

    # def on_epoch_end(self):
    #     'Updates indexes after each epoch'
    #     self.indexes = np.arange(len(self.list_IDs))
    #     if self.shuffle == True:
    #         np.random.shuffle(self.indexes)


# //--------------------------------------------
# //--------------------------------------------
