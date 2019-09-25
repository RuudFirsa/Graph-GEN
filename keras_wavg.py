"""
Copyright 2019 Ruud van Deursen, Firmenich SA.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

"""
File keras_wavg.py defines a layer to compute
the weighted average for a list of equally sized
layers. The method defines a trainable weight
with a single weight per layer. The output
is then simply the dot-product between the
vector with the trainable weights and the 
list of layers.
"""

from keras.layers import Layer,Concatenate,Reshape,Flatten
from keras import backend as K

class WeightedAverage(Layer):
    """ 
    Class WeightedAverage defines a layer of trainable weights for a list of layers.
    The class computes a vector with length of number of layers and then takes
    the dot-product to compute a single average output layer.
    """
    
    def __init__(self, **kwargs):
        """ 
        Constructor of WeightedAverage. The methods compute
        the shape based on the input shape.
        Input:
        kwargs -- Arguments.
        """
        super(WeightedAverage, self).__init__(name="WeightedAvg",**kwargs)

    def build(self, input_shape):
        """
        Method defines a trainable vector to compute
        the weighted average.
        Input:
        input_shape -- Input shape of the object.
                       This shape is converted to number
                       of weights and output dimension.
        """
        # Extract dimensions
        self.output_dim = input_shape[0][1]
        self.num = len(input_shape)
        # Define a vector with trainable weights
        self.kernel = self.add_weight(name='kernel', 
                                      shape=(self.num,1),
                                      initializer='uniform',
                                      trainable=True)
        # Build the thing
        super(WeightedAverage, self).build(input_shape)  # Be sure to call this at the end

    def call(self, x):
        """ 
        Method computes the output tensor.
        Input:
        x -- List with input tensors.
        Return:
        Layer with output values.
        """
        # Concatenate, Reshape, take the inner product and flatten to a single layer
        # (this is way easier than writing some sort of lambda thing)
        C = Concatenate()(x)
        R = Reshape(target_shape=(self.output_dim,self.num))(C)
        inner = K.dot(R,self.kernel)
        return Flatten()(inner) # We need to flatten to inner product to get correct dimensions

    def compute_output_shape(self, input_shape):
        """
        Method computes the output shape, which equals
        the shape of the input layers.
        Return:
        Output shape of the layer.
        """
        return (None,self.output_dim)
