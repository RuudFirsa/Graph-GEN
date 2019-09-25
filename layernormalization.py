"""
File defines the class LayerNormalization.
This class is an implementation of the method
proposed in the following publication:

https://arxiv.org/abs/1607.06450
"""

from keras.layers import Layer
import keras.backend as K
from keras.initializers import Ones, Zeros

class LayerNormalization(Layer):
    """
    Class LayerNormalization defines a normalization
    to normalize computing the mean and variance used for 
    normalization from all of the summed inputs to the neurons 
    in a layer on a single training case
    """
    
    def __init__(self, eps=1e-6, **kwargs):
        self.eps = eps
        super(LayerNormalization, self).__init__(**kwargs)
    
    def build(self, input_shape):
        self.gamma = self.add_weight(name='gamma', shape=input_shape[-1:],
                                     initializer=Ones(), trainable=True)
        self.beta = self.add_weight(name='beta', shape=input_shape[-1:],
                                    initializer=Zeros(), trainable=True)
        super(LayerNormalization, self).build(input_shape)
    
    def call(self, x):
        mean = K.mean(x, axis=-1, keepdims=True)
        std = K.std(x, axis=-1, keepdims=True)
        return self.gamma * (x - mean) / (std + self.eps) + self.beta
    
    def compute_output_shape(self, input_shape):
        return input_shape
