"""
Source : http://sabermetricinsights.blogspot.com/2014/05/bayesian-linear-
         regression-with-pymc.html

Date : 11 August 2016
"""

### Imports 
import numpy as np
import matplotlib.pyplot as plt
import pymc as bae
import pandas as pd

### We are trying to model the mpg of a car as a function of the 
### car's weight. So data frame is called 'float_df'

### We are trying to solve simple linear regression given by:
### y = b0 + b1(x) + error
### where b0 is the intercept term, b1 is the slope, and error is error

### Model the intercept/slope terms of our model as normal random variables
### with comically large variances
b0 = bae.Normal('b0',0,0.0003)
b1 = bae.Normal('b1',0,0.0003)

### Model our error term as a uniform random variable
err = bae.Uniform('err',0,500)

### "model" observed x values as a normal random variable
### in reality because x is observed, it doesn't really matter
x_weight = bae.Normal('weight',0,1,value=np.array(float_df['weight']),
                      observed=True)