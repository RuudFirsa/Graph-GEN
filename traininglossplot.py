"""
Copyright 2019 Ruud van Deursen, Firmenich SA.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


"""
File training loss plot defines a live plotter for 
values during training. Thie file takes the loss
by default and can plot other values.
The plotter is herein modified to work with
the LSTM-generator.
TODO: Small updates.
"""

import IPython
from keras.callbacks import Callback
import matplotlib.pyplot as plt
from IPython.display import clear_output
import numpy
from datetime import datetime

# Define a Callback to make a simple updating loss-plot
class TrainingLossPlot(Callback):
    """
    Class TrainingLossPlot defines a simple plotting
    function to interactively follow the loss function.
    """
    
    def __init__(self,num_epochs,plot=["loss"],other=[]):
        """
        Constructor of TrainingLossPlot defining
        the number of epochs.
        Input:
        num_epochs -- Number of epochs.
        """
        super().__init__()
        self.num_epochs = num_epochs
        self.plot = plot
        self.other = other
        
    def Plot(self):
        """
        Method plots the current status of the loss-function.
        """
        clear_output(wait=False)
        other=self.other
        plt.figure(figsize=(12,4))
        plt.subplot(1,len(other)+1,1)
        plt.title("Epoch %s/%s"%(self.epoch,self.num_epochs))
        if self.epoch > 1:
            X = numpy.arange(self.epoch)+1
            Y = numpy.array(self.losses["loss"])
            plt.plot(X,Y)
        plt.xlabel("Epoch")
        plt.ylabel("Training loss")
        
        # Add the percentage valid
        if len(other)>0:
            plt.subplot(1,len(other)+1,2)
            plt.title("Epoch %s/%s"%(self.epoch,self.num_epochs))
            hist = other[0].history
            X = sorted(hist.keys())
            Y = [hist[key] for key in sorted(hist.keys())]
            plt.plot(X,Y)
            Y = [other[0].minvalue for key in sorted(hist.keys())]
            Yup = [other[0].uppervalue for key in sorted(hist.keys())]
            Ylow = [other[0].lowervalue for key in sorted(hist.keys())]
            plt.plot(X,Ylow,'--',c='k')
            plt.plot(X,Y,'-',c='k')
            plt.plot(X,Yup,'--',c='k')
            plt.xlabel("Epoch")
            plt.ylabel("Pct valid")
        plt.show()
    
    def on_train_begin(self, logs={}):
        """
        Method defines the initial parameters.
        Input:
        logs -- Initial logs.
        """
        self.losses = {}
        self.epoch = 0
        for fct in self.plot:
            self.losses[fct] = []
        self.Plot()
        
    def on_epoch_end(self, batch, logs={}):
        """
        Method updates on epoch end.
        Input:
        batch -- Batch.
        logs  -- Log.
        """
        for fct in self.plot:
            self.losses[fct].append(logs.get(fct))
        self.epoch += 1
        self.Plot()


# Define a Callback to compute the elapsed time
class Timer(Callback):
    """
    Simple timer to indicate the time/epoch.
    """
    
    def on_epoch_begin(self,batch,logs={}):
        """
        Method stores the start time.
        """
        self.start = datetime.now()
        
    def on_epoch_end(self,batch,logs={}):
        """
        Method computes and prints the elapsed time.
        """
        self.end = datetime.now()
        diff = self.end-self.start
        diff_in_seconds = diff.days*24*60*60 + diff.seconds
        print("Epoch time: %.3f seconds"%(diff_in_seconds))
