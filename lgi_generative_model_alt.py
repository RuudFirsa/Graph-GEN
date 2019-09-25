"""
Copyright 2019 Ruud van Deursen, Firmenich SA.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

"""
File lgi_generative_model_alt defines an alternate
sequence for the architecture, applying a layer
normalization after every encoding layer 
"""

# Keras inputs
from __future__ import print_function
from keras import backend as K
from keras.models import Model
from keras.layers import Input,CuDNNLSTM,CuDNNGRU,Dense, Activation, Dropout,concatenate,average
from keras.layers.wrappers import Bidirectional
from layernormalization import LayerNormalization
from keras.optimizers import Adam,SGD
from keras.utils.data_utils import get_file
from keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from traininglossplot import TrainingLossPlot,Timer
from lgi_generator import Generator,OnlineGenerator
from keras_wavg import WeightedAverage # Weighted average (learnable)

# Miscellaneous inputs
import numpy as np
import random
import sys
from datetime import datetime
from numpy import random,zeros

# Section with a model trainer
# TODO: Integrate with the new architecture
class Trainer:
    """
    Class ErtlLSTMTraininer defines the training workflow to train
    the LSTM-SMILES generator.
    """
    
    def __init__(self,model,Utils,sanitycheck=None,minimum=.95,patience=10,verbose=False):
        """
        Constructor of ErtlLSTMTraininer defining
        the model to train.
        Input:
        model -- Model to be trained.
        Utils -- Utils
        ...
        """
        self.model = model
        if sanitycheck is not None:
            self.gen = Generator(model,Utils,sanitycheck=sanitycheck)
        else:
            self.gen = None
            
        self.model.compile(loss="categorical_crossentropy",optimizer=Adam(lr=0.003),metrics=[])
        self.verbose = verbose
        self.minimum = minimum
        self.patience = patience
        self.interactive = True

    def SetInteractive(self,interactive):
        """
        Method updates the flag for the interactive mode.
        Input:
        interactive -- New value for the flag.
        Return:
        Updated instance of Trainer.
        """
        self.interactive=interactive
        return self        
        
    def Fit(self,X,y,filepath,logfile,num_epochs=100,batch_size=256,ncollect=180,npop=None,mycallbacks=list(),verbose=0):
        """
        Method fits the dataset.
        Input:
        X              -- Vector with X
        y              -- Vector with y
        filepath       -- File path for the model files.
        logfile        -- File to save the log to.
        num_epochs     -- Number of epochs to train at each iteration (default = 100).
        batch_size     -- Batch sized used for training (default = 256).
        ncollect       -- Number to generate after each epoch (default = 300).
                          Used to determine the statistically robust window.
        npop           -- Population size (default = None, assumed very large).
                          Used to determine the statistically robust window.
        mycallbacks    -- Custom callback functions, will be added to the list
                          of default callback functions.
        verbose        -- Flag for verbose.
        """ 
        # Define a set of callback method to save models, stop early and monitor progress
        checkpoint = ModelCheckpoint(filepath, monitor='loss', verbose=1, save_best_only=True, mode='min')
        callbacks_list = list([checkpoint])
        callbacks_list.extend(mycallbacks)
        if self.gen is not None:
            # Follow percentage valid
            earlystopping = OnlineGenerator(self.gen,
                                            ncollect=ncollect,
                                            npop=npop,
                                            min_value=self.minimum,
                                            patience=self.patience,
                                            verbose=verbose,
                                            restore_best_weights=True,
                                            filename=logfile)
            callbacks_list.append(earlystopping)
            if self.interactive:
                monitoring = TrainingLossPlot(num_epochs=num_epochs,other=[earlystopping])
                callbacks_list.append(monitoring)
        else:
            # Follow regular loss
            earlystopping = EarlyStopping(monitor='loss', min_delta=0, patience=5, verbose=1, mode='auto', baseline=None, restore_best_weights=True)
            callbacks_list.append(earlystopping)
            if self.interactive:
                monitoring = TrainingLossPlot(num_epochs=num_epochs)
                callbacks_list.append(monitoring)
                
        # Fit
        history = self.model.fit(X, y, batch_size=batch_size, epochs=num_epochs, callbacks=callbacks_list,shuffle=True, verbose=verbose)
            
        # Done
        return self.model,history

    
# Class BaseModel with reusable methods
class BaseModel:
    """
    Class ErltLSTMModel defines a class to setup LSTM models
    to train SMILES-generators. This class defines the shared
    method and architecture between LSTMModel and GRUModel. 
    """
    
    def __init__(self,Utils,Unit=CuDNNLSTM,Layers=[256,256],Bidirectional=[True,False],Dropout=0.3,Optimizer=Adam(lr=0.002),Loss="categorical_crossentropy",minimodels=4,split=0,merge=0):
        """
        Constructor of ErtlLSTMModel.
        Input:
        Utils         -- ErtlLSTMUtils, defining MaxLen and NumChars.
        Unit          -- Network training layer (e.g. CuDNNLSTM or CuDNNGRU; default = CuDNNLSTM).
        Layers        -- Layer sizes for the LSTM layers (default = [256,256]).
        Bidirectional -- Flags to ask bidirectional LSTM layers (default = [True,True]).
        Dropout       -- Dropout (default = 0.2).
        Optimizer     -- Optimized (default = Adam with lr=0.002)
        Loss          -- Loss function (default = 'categorical crossentropy').
        num_encoding  -- Number of encoding layers (default = 4).
        split         -- Splitting scheme for minimodels (default = 0):
                         Choose from:
                         0: Encoding only
                         1: Embedding+encoding
        merge         -- Merge mode of the model (default = 0).
                         Choose from:
                         0: Concatenated
                         1: Average
                         2: Learnable weighted average
        """
        super(BaseModel,self).__init__()
        
        # Collect the number of characters and the maximum
        # observed length during data preparation.
        self.Utils = Utils
        self.Unit = Unit
        
        # Cache the values
        self.Layers = Layers
        self.bilstm = Bidirectional
        self.Dropout = Dropout
        self.Optimizer = Optimizer
        self.Loss = Loss
        self.model = None
        self.num_models = minimodels
        self.merge = merge
        self.split = split
        
    def Init(self,weightsfile=None,verbose=False):
        """
        Method initializes the model based on the specified
        parameters.
        Input:
        weightsfile  -- File with weights (default is None).
        sanitycheck  -- Validation check.
        verbose      -- Flag for verbose mode, printing architecture (default = False).
        Return:
        Initialized model based on the specified parameters.
        """
        # Load the variables
        num_chars,maxlen = self.Utils.NumChars(),self.Utils.MaxLen()
        l1,l2 = self.Layers[0],self.Layers[1]
        dropout = self.Dropout
        
        # Build the model using the specified unit
        # Example units are for instance CuDNNLSTM or CuDNNGRU
        # The unit has been introduced as variable to facilitate flexible modifications.
        comment_seq = Input(shape=[maxlen,num_chars],name="Input")

        # Define image
        minimodels = []
        if self.split == 1:
            # Apply scheme 1: Multiple embedding and multiple encoding            
            for idx in range(self.num_models):
                if self.bilstm[0]:
                    output_i = Bidirectional(self.Unit(l1,return_sequences=True),name="Embedding")(comment_seq)
                else:
                    output_i = self.Unit(l1,return_sequences=True,name="Embedding")(comment_seq)
                if self.bilstm[1]:
                    output_i = Bidirectional(self.Unit(l2),name="Latent_%s"%(idx))(output_i)
                else:
                    output_i = self.Unit(l2,name="Latent_%s"%(idx))(output_i)
                output_i = LayerNormalization(name="LayerNormm_%s"%(idx))(output_i)
                minimodels.append(output_i)    
                
        else:
            # Apply scheme 0: One embedding and multiple encoding
            if self.bilstm[0]:
                output = Bidirectional(self.Unit(l1,return_sequences=True),name="Embedding")(comment_seq)
            else:
                output = self.Unit(l1,return_sequences=True,name="Embedding")(comment_seq)

            # Create multiple encoding models
            minimodels = []
            for idx in range(self.num_models):
                if self.bilstm[1]:
                    output_i = Bidirectional(self.Unit(l2),name="Latent_%s"%(idx))(output)
                else:
                    output_i = self.Unit(l2,name="Latent_%s"%(idx))(output)
                output_i = LayerNormalization(name="LayerNormm_%s"%(idx))(output_i)
                minimodels.append(output_i)            
            
        # Combine to a single model using the selected mode
        if len(minimodels)==1:
            output = minimodels[0]
        elif self.merge == 1:
            output = average(minimodels)
        elif self.merge == 2:
            output = WeightedAverage()(minimodels)
        else:
            # Default is concatenate
            output = concatenate(minimodels)
            
        # Apply a dropout and compute the probabilities
        output = Dropout(self.Dropout)(output)
        output = Dense(num_chars,name="Output")(output)
        output = Activation("softmax")(output)
        
        # Compile the model using the specified methods
        self.model = Model([comment_seq],output)       
        if verbose:
            self.model.summary()        
        
        # Load the weight if specified
        if weightsfile is not None:
            self.model.load_weights(weightsfile)
        
        # Done
        return self.model
    
    def InitTrainer(self,
                    weightsfile=None,
                    sanitycheck=None,
                    minimum = 0.95,
                    patience=10,
                    verbose=False):
        """
        Method initializes a model and immediately
        generates a training instance for the model.
        Input:
        weightsfile  -- File with weights (default is None).
        sanitycheck  -- Validation check (default is None).
                        The sanity check validates the generated
                        objects using an on-training generator.
                        If set to None, no validation is carried
                        out during training.        
        minimum      -- Minimum value to activate earlystopping method.
                        The minimum value defines the minimum 
                        percentage that should be reached.
        verbose      -- Flag for verbose mode (default = False).
        Return:
        Initialized model based on the specified parameters.
        """        
        return Trainer(self.Init(weightsfile),self.Utils,sanitycheck,minimum,patience,verbose)
    
    def InitGenerator(self,weightsfile=None,
                      sanitycheck=lambda x: True):
        """
        Method generates an instance to generate SMILES.
        Input:
        weightsfile  -- File with weights for the network.
        regsmiles    --
        Return:
        Instance to generate SMILES using the trained model.
        """
        # Get the instance from the cache and load the weights if specified
        if self.model is not None:
            model = self.model
            if weightsfile is not None:
                self.model.load_weights(weightsfile)
                
        # Generate a new instance
        else:
            model = self.Init(weightsfile=weightsfile)
            
        # Construct the generator
        return Generator(model,self.Utils,sanitycheck)
        
        
# Import unittest
import unittest
class BaseModelTest(unittest.TestCase):
    """ """
    
    def setUp(self):
        pass

class BiLSTMModelTest(unittest.TestCase):
    """ """
    
    def setUp(self):
        pass

class BiGRUModelTest(unittest.TestCase):
    """ """
    
    def setUp(self):
        pass
    
