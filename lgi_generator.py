"""
Copyright 2019 Ruud van Deursen, Firmenich SA.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


"""
File lgi_generator defines a class to generate graphs in LGI-format
"""

# Keras inputs
from __future__ import print_function
from keras import backend as K
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.layers import CuDNNLSTM, CuDNNGRU # Import for both model types
from keras.optimizers import Adam,SGD
from keras.utils.data_utils import get_file
from keras.callbacks import Callback
from traininglossplot import TrainingLossPlot,Timer
from samplesize import SampleSize

# Miscellaneous inputs
import numpy as np
import random
import sys
from datetime import datetime
from numpy import random,zeros,array


# Define a special function for early-stopping
# This method will stop the training when a better percentage
# has been reached.
class OnlineGenerator(Callback):
    """
    The class OnlineGenerator is a special Callback function
    to keep track of the number of valid generated molecules.
    The Callback has been derived from EarlyStopping and stops
    as soon as a threshold is repeatedly exceeded.
    """

    def __init__(self,
                 gen,
                 min_value=.95,
                 ncollect=300,
                 npop=None,
                 patience=0,
                 verbose=0,
                 baseline=None,
                 targetdist=None,
                 distmethod=None,
                 warmup=3,
                 restore_best_weights=True,
                 filename=None):
        """
        Constructor of OnlineGenerator.
        Input:
        gen                  -- Generator.
        min_value            -- Minimum value to reach and repeatedly hold (default = 0.95).
        patience             -- Patience until stop.
        verbose              -- Value for verbose mode (default = 0).
        baseline             -- Value for baseline (default = None).
        targetdist           -- Target distribution (default = None).
                                If set to None, the distribution of online generated
                                entries is not compared to the reference distribution.
        distmethod           -- Distribution method (default = None).
                                If set to None, the distribution of online generated
                                entries is not compared to the reference distribution.
        restore_best_weights -- Value to restore best weights (default = True).
        filename             -- File name to store the output data to (default is None).
                                If set to None, nothing is stored.
        """
        super(OnlineGenerator, self).__init__()
        self.gen = gen
        self.ncollect = ncollect
        self.targetdist = targetdist
        self.method = distmethod
        self.baseline = baseline
        self.patience = patience
        self.pct,self.epoch = 0,0
        self.mode = 'max'
        self.verbose = verbose
        self.history = {}
        self.wait = 0
        self.stopped_epoch = 0
        self.restore_best_weights = restore_best_weights
        self.best_weights = None
        self.monitor_op = np.greater
        self.filename = filename
        
        # Define minvalue with statistical window
        self.minvalue = min_value
        self.lowervalue,self.uppervalue = SampleSize.ComputeInterval(min_value,ncollect,npopulation=npop)

    def on_train_begin(self, logs=None):
        # Allow instances to be re-used
        self.wait = 0
        self.stopped_epoch = 0
        self.loss_hist = list()
        self.pct_hist = list()
        if self.baseline is not None:
            self.best = self.baseline
        else:
            self.best = np.Inf if self.monitor_op == np.less else -np.Inf

    def generate(self):
        self.epoch += 1
        self.pct,self.H = self.gen.PercentValid(ncollect=self.ncollect,verbose=self.verbose) 
        self.pct_hist.append(self.pct)
        self.history[self.epoch] = (self.pct,self.H)
        if self.verbose:
            print("Updating score after epoch %s: %.1f%%"%(self.epoch,100.0*self.pct))

        # Check the output valid
        current = self.pct
        if current is None:
            return

        # Update the best weights
        if current > self.best:
            self.best = current
            self.best_weights = self.model.get_weights()

        # Continue until the lower confidence bound is repeatedly exceeded
        if self.minvalue is not None and current < self.lowervalue:
            self.wait = 0
            if self.restore_best_weights:
                self.best_weights = self.model.get_weights()
            print("Min value of %.4f not yet reached: %.4f"%(self.minvalue,current))
            return
        else:
            # Count and stop after 5 rounds.
            self.wait += 1
            print("Above min value of %.4f for %s consecutive epochs"%(self.minvalue,self.wait))
            if self.wait >= self.patience:
                self.stopped_epoch = self.epoch
                self.model.stop_training = True
                if self.restore_best_weights:
                    if self.verbose > 0:
                        print('Restoring model weights from the end of '
                              'the best epoch')
                    self.model.set_weights(self.best_weights)        
            
    def on_epoch_end(self, batch, logs={}):
        self.generate()
        self.loss_hist.append(logs["loss"])
        print("l",self.loss_hist)
        print("p",self.pct_hist)

    def on_train_end(self, logs={}):
        #self.generate()
        self.collected_logs = logs
        if self.stopped_epoch > 0 and self.verbose > 0:
            print('Epoch %05d: early stopping' % (self.stopped_epoch + 1))
            lower,upper = self.epoch-self.patience,self.epoch-self.patience+2
            print("Recommended model at epoch = [%s,%s]"%(lower,upper))
        else:
            print("No stable results observed for target %.3f [%.3f,%.3f]"%(self.minvalue,self.lowervalue,self.uppervalue))
            
        if self.filename is not None:
            with open(self.filename,"w") as f:
                f.write("epoch,loss,pct_valid\n")
                e = 0
                for loss,pct in zip(self.loss_hist,self.pct_hist):
                    e += 1
                    f.write("%s,%.4f,%.1f\n"%(e,loss,100.0*pct))
            print("Loss/Valid written to %s"%(self.filename))
            
    def get_monitor_value(self, logs={}):
        monitor_value = self.pct
        if monitor_value is None:
            warnings.warn(
                'Early stopping conditioned on metric `%s` '
                'which is not available. Available metrics are: %s' %
                (self.monitor, ','.join(list(logs.keys()))), RuntimeWarning
            )
        return monitor_value


# Section with a model generator
class Generator:
    """
    Class Generator defines a generator using an LSTM-model
    to generate SMILES. 
    """
    
    def __init__(self,model,Utils,sanitycheck=lambda x: True):
        """
        Constructor of ErtlLSTMGenerator.
        Input:
        model        -- Trained model.
        Utils        -- Utils used to generate the training set.
        """
        self.model = model
        self.Utils = Utils
        self.sanitycheck = sanitycheck
        
    def Sample(self,preds):
        """
        Method samples an index from the probability array.
        Input:
        preds -- Prediction array.
        Return:
        Randomly sampled index.
        """
        # sampling an index from a probability array
        preds = np.asarray(preds).astype('float64')
        preds = np.log(preds)
        exp_preds = np.exp(preds)
        preds = exp_preds / np.sum(exp_preds)
        probas = np.random.multinomial(1, preds, 1)
        c = np.argmax(probas)
        return c
        
    def DoSmiles(self,seedstring,maxlen=40,smi=""):
        """
        Method generates a SMILES using the model.
        Input:
        seedstring -- Seed.
        maxlen     -- Maximum length (default = 40).
        smi        -- Currently known SMILES String.
                      New characters will be added.
        Return:
        Updated SMILES String.
        """
        model,Utils = self.model,self.Utils
        while (True):
            x = zeros((1, Utils.maxlen, len(Utils.chars)))
            for t, char in enumerate(seedstring):
                x[0, t, Utils.char_indices[char]] = 1.

            preds = model.predict(x, verbose=0)[0]
            next_index = self.Sample(preds)
            next_char = Utils.indices_char[next_index]

            seedstring = seedstring[1:] + next_char

            if next_char == "\n":
                return smi,seedstring
            smi += next_char
            if len(smi) > 120: return None,None # new seed needed        
        
    def Predict(self,
                ncollect=1000,
                ncopies=20,
                verbose=False,
                sanitycheck=lambda x: True,
                standardize=lambda x: x):
        """
        Class Predict generates SMILES.
        Input:
        ncollect    -- Number of SMILES to collect (default = 1,000).
        ncopies     -- Number of copies for every random seed (default = 20).
        verbose     -- Flag for verbose mode (default = False).
        sanitycheck -- Method to check generated Strings for
                       validity (default = 'lambda x: True', accepting all).
        standardize -- Method to standardize the molecule upon completion 
                       (default = lambda x: x, keeping unchanged).
        """
        model,maxlen = self.model,self.Utils.MaxLen()
        text = self.Utils.Text()
        nsmi = 0
        mols = list()
        good,bad = 0,0
        starttime = datetime.now()
        
        # Run as long as we have too few SMILES on the list
        while len(mols)<ncollect:
            start = random.randint(maxlen, len(text) - 1)
            end = text.index("\n",start)
            seedstring= text[end-maxlen+1:end+1]

            # Run the number of copies on this seed
            for n in range(ncopies):
                smi,seedstring = self.DoSmiles(seedstring)
                if smi == None:
                    if seedstring == None:
                        break; # new seed needed

                # Decode to molecule and check if valid
                keep = True
                if sanitycheck is not None:
                    keep = sanitycheck(smi)
                else:
                    keep = self.sanitycheck(smi)
                if keep:
                    # Count the molecule as passed 
                    good += 1
                    # Append the standardized SMILES
                    mols.append(standardize(smi))
                    nsmi += 1
                    if verbose:
                        print(nsmi,"Rate G/B = %s/%s"%(good,bad),smi)
                    
                    # Stop on completion
                    if len(mols) == ncollect:
                        break
                else: 
                    bad += 1
                    
        # Compute the elapsed time and print
        if verbose:
            diff = datetime.now()-starttime
            diff_in_seconds = diff.days*24*60*60 + diff.seconds
            print("Generation time: %s seconds"%(diff_in_seconds))
        
        # Done
        return mols     

    def PercentValid(self,ncollect=180,ncopies=5,distmethod=None,verbose=False):
        """
        Method validates the generation rate of valid strings.
        Input:
        ncollect    -- Number to generate (default = 300).
        ncopies     -- Number of random copies to generate (default = 5).
        distmethod  -- 
        verbose     -- Flag for verbose mode (default = False).
        """
        # Call the parent class to generate and accept all
        # We accept all with True because we want to measure
        mols = self.Predict(ncollect=ncollect,ncopies=ncopies,sanitycheck=lambda x: True)
        num_gen = float(len(mols))
        # Filter the valid ones using our real measure
        valid = [smi for smi in mols if self.sanitycheck(smi)]
        num_valid = float(len(valid))
        ratio = num_valid/num_gen
        print("Validity: NumGen=%s NumValid=%s Ratio=%.4f"%(num_gen,num_valid,ratio))
        
        # Compute a histogram
        H = None
        if distmethod is not None:
            values = [distmethod(smi) for smi in mols]
            H = PosIntHist(values)
        
        # Compute the ratio
        return ratio,H

# Test time
import unittest
class GeneratorTest(unittest.TestCase):
    """ Method checks GeneratorTest for correct return values """
    
    def setUp(self):
        pass
    
class ValidationGeneratorTest(unittest.TestCase):
    """ Methods checks ValidationGenerator for correct return values """
    
    def setUp(self):
        pass
