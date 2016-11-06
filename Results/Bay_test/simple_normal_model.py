import pymc as bay
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal

### Generate Random Data
data = [normal(100,20) for i in xrange(1000)]

#########EXAMPLE 1
#### Look to see if it is normally distributed
#plt.hist(data)
#
#### Don't know mean but believe it is between max and min of data
#mean = bay.Uniform('mean',lower=min(data),upper=max(data))
#
#### Don't know precision so prior belief is it is in between 0.0001 and 1.0
#### with belief values in this range are uniform
#precision = bay.Uniform('precision',lower=0.0001,upper=1.0)
#
#### Run a model around mean and precision
#process = bay.Normal('process',mu=mean,tau=precision,value=data,observed=True)



#########EXAMPLE 2
mean = bay.Uniform('mean',lower=min(data),upper=max(data))

std_dev = bay.Uniform('std_dev',lower=0,upper=50)

@bay.deterministic(plot=False)
def precision(std_dev = std_dev):
    return 1.0/(std_dev * std_dev)
    
process = bay.Normal('process',mu=mean,tau=precision,value=data,observed=True)