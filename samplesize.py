"""
Copyright 2019 Ruud van Deursen, Firmenich SA.
File is inspired by https://github.com/shawnohare/samplesize/blob/master/samplesize.py.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import argparse
from scipy.stats import norm
from numpy import ceil,sqrt
import numpy as np

class SampleSizeCalculator:
    """
    Class SampleSizeCalculator computes the minimal number of 
    datapoints needed to get a statistically significant result.
    """    
    
    def __init__(self):
        """
        Constructor of SampleSizeCalculator, defining
        the z-dictionary.
        """
        self.zdict = { .90: 1.645,.91: 1.695,.99: 2.576,.97: 2.17,.94: 1.881,.93: 1.812,.95: 1.96,.98: 2.326, .96: 2.054,.92: 1.751 }
        
    def GetZ(self,confidence_level):
        """
        Method looksup or computes the Z-value used
        to estimate the population size.
        Input:
        confidence_level -- Level of confidence wanted.
        Return:
        Z-value for the confidence level.
        """
        if confidence_level in self.zdict:
            z = self.zdict[confidence_level]
        else:
            alpha = 1 - (confidence_level)
            z = norm.ppf(1 - (alpha/2))
        return z
        
    def ComputeSampleSize(self,population_size,margin_error=.025,confidence_level=.99,sigma=.50):
        """
        Method computes the minimum population size needed.
        Input:
        population_size  -- Size of the population.
        margin_error     -- Error margin (default = 0.025)
        confidence_level -- Level of confidence wanted (default = .99)
        sigma            -- Assumed standard deviation in the population (default = .50)
        Return:
        Minimum population size needed, rounded up to the next integer.
        """
        # Make it fractional
        if margin_error > 1:
            margin_error /= 100.0
        if confidence_level > 1:
            confidence_level /= 100.0
        
        # Compute z and get the number output
        z = self.GetZ(confidence_level)
        N = population_size
        M = margin_error
        num = z**2 * sigma**2 * (N / (N-1))
        den = M**2 + ((z**2 * sigma**2)/(N-1))
        return int(ceil(num/den))
    
    def ComputeInterval(self,percentile,nsample,npopulation=None,confidence_level=.95):
        """
        Method computes the percentage interval for the
        given sample and population sizes at the requested
        confidence level. 
        Input:
        percentile       -- Desired percentile.
        nsample          -- Sample size.
        npopulation      -- Population size.
        confidence_level -- Confidence level (default = .99)
        """
        n,p = float(nsample),float(percentile)
        z = self.GetZ(confidence_level)
        percratio = p*(1-p)/float(n)
        if npopulation is not None:
            # Apply exact formula
            N = float(npopulation)
            popratio = (N-n)/(N-1)
            factor = z*sqrt(percratio*popratio)
        else:
            # Assume the population is very large and (N-n)/(N-1) ~ 1
            factor = z*sqrt(percratio)
            
        # Return the interval as p-factor,p+factor
        return np.max([p-factor,0.0]),np.min([p+factor,1.0])

    
# Define a statically constructed method
SampleSize = SampleSizeCalculator()

import unittest
class SampleSizeTest(unittest.TestCase):
    """
    Test class evaluating the correct outcome of the program.
    """
    
    def setUp(self):
        """ Setup """
        pass
    
    def test_Value(self):
        """ Method tests a value. """
        pop,exp = 215922,2623
        self.assertEqual(exp,SampleSize.ComputeSampleSize(pop,margin_error=0.025))
        self.assertEqual(exp,SampleSize.ComputeSampleSize(pop,margin_error=2.5))
        pop,exp = 1000000,999
        self.assertAlmostEqual(exp,SampleSize.ComputeSampleSize(pop,margin_error=0.031,confidence_level=.95),0)
        self.assertAlmostEqual(exp,SampleSize.ComputeSampleSize(pop,margin_error=3.1,confidence_level=95),0)
        
    def test_ComputeInterval(self):
        """ Method test the correct computation of lower and upper limit for 95% CI """
        nsample = 1406
        lower,upper = 0.4739,0.5261
        act_lower,act_upper = SampleSize.ComputeInterval(0.5,1406)
        self.assertAlmostEqual(lower,act_lower,4)
        self.assertAlmostEqual(upper,act_upper,4)
        
# Runner for tests
if __name__ == "__main__":
    unittest.main()
