"""
Copyright 2019 Ruud van Deursen, Firmenich SA.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

"""
File lgi_generative_model_utils.py defines a file
to prepare the dataset for use with the 
g6_generative_model.py

For the presented dataset the method prepares
the feature vector and labels to train from.
"""

import numpy
from numpy import random

#######################################
# Section with a class defining utils #
#######################################
class DataUtils:
    """
    Class DataUtils defines a class for data preparation
    to convert SMILES String to a training set.
    """
    
    def __init__(self,maxlen=42,step=3):
        """
        Constructor of ErtlLSTMUtils.
        Input:
        maxlen -- Maximum length of the g6-string (default = 42).
        step   -- Step size used (default = 3).
        """
        # Define the values
        self.maxlen = maxlen
        self.step = step
        
        # Store the allowed characters in g6
        self.okchars = "ABCDEF()%123456789\n"
        self.chars = list(self.okchars)
        
        # Define the defaults
        self.char_indices = dict([(c, i) for i, c in enumerate(self.chars)])
        self.indices_char = dict([(i, c) for i, c in enumerate(self.chars)])

    
    def Prepare(self,smilist,clear=True,updateChars=True,augment=None,naug=0,shuffle=True):
        """
        Method prepares a list with SMILES String and returns
        a list with line-separated SMILES Strings for training.
        The SMILES will be subjected to characters replacement
        Input:
        smilist      -- List with multiple SMILES.
                        This can be a list of SMILES or a string
                        defining a filename
        clear        -- Flag to clear input (default = True).
        augment      -- Method to augment the entry for an identical
                        entry represented with a different string.
                        Method syntax: augment(text,naug), where
                        text is the input string and naug the number
                        of augmented entries (default is None, not executed).
        naug         -- Number of augmentations (default = 0, not executed).
        shuffle      -- Flag to shuffle data (default = True).
        Return:
        List with fine-tuned SMILES.
        """
        # Check the type of smilist - if string consider this
        # this variable a file with SMILES and read all lines
        if type(smilist) == str:
            with open(smilist) as f:
                data = [line.strip() for line in f.readlines()]
        else:
            data = smilist
        
        # Augment/randomize if instructed to do so
        if naug>0 and augment is not None:
            data = [aug for smi in data for aug in augment(smi,naug)]  
        self.evaluated = len(data)

        # Extract the max size
        self.maxg6 = max([len(d) for d in data])
            
        # Shuffle the data
        if shuffle:
            random.shuffle(data)
            
        # Create the replacement and return as a single
        # string with a line separator between the words.
        keep = list(filter(lambda x: len(x)>0 and len(x)<=self.maxlen,data))
        self.kept = len(keep)
        self.text = "\n".join(keep)+"\n" # Add a breaker for the last line.
        
        # Clear input to save memory
        if clear:
            data.clear()

        # Cache the data and return the joined list
        self.data = data
        return self.text
    
    def Encode(self,text):
        """
        Method encodes a list of Strings to vectors for the LSTM.
        Input:
        text   -- Text to be translated.
        maxlen -- Maximum sentence length (default = 40).
                  This length will be set to the longest observed
                  length, if the longest observed length is lower
                  than maxlen.
        step   -- Step (default = 3).
        """
        # Update maxlen based on the longest observed word
        # This is important, otherwise no sentences will be presented
        # to the network.
        maxlen,step = self.maxlen,self.step
        maxword = numpy.max([len(x) for x in text.split("\n")])
        maxlen = numpy.min([maxword,maxlen])
        self.maxlen = maxlen
        
        # Translate the sentences to little pieces.
        sentences = []
        next_chars = []
        for i in range(0, len(text) - maxlen, step):
            sentences.append(text[i: i + maxlen])
            next_chars.append(text[i + maxlen])

        # Vectorize the sentences in boolean format
        X = numpy.zeros((len(sentences), maxlen, len(self.chars)), dtype=numpy.bool)
        y = numpy.zeros((len(sentences), len(self.chars)), dtype=numpy.bool)
        for i,sentence in enumerate(sentences):
            for t,char in enumerate(sentence):
                X[i,t,self.char_indices[char]] = True
            y[i,self.char_indices[next_chars[i]]] = True

        # Done: Return X,Y and maxlen
        return X,y,maxlen
        
    def NumChars(self):
        """
        Return:
        Number of cached characters.
        """
        if self.chars is None:
            return 0
        return len(self.chars)
    
    def MaxLen(self):
        """
        Return:
        Maximum observed length.
        """
        if self.maxlen is None:
            return 0
        return self.maxlen
    
    def Text(self):
        """
        Return:
        The text corpus.
        """
        return self.text
    
    def LenText(self):
        """
        Return:
        Length of the cached text.
        """
        return len(self.text)    

            
#############################################
# Section with main method to run tests cmd #
#############################################
import unittest
if __name__ == "__main__":
    unittest.main(argv=['first-arg-is-ignored'],  verbosity=2, exit=False)
