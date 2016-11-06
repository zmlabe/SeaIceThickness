"""
Practice bayesian modeling for March sit data using pymc3 round 2

example : http://twiecki.github.io/blog/2014/03/17/bayesian-glms-3/

 
Author : Zachary Labe
Date : 24 August 2016
"""

### Read in modules -----> using new PyMC3
import numpy as np
import pymc3 as pm
from mpl_toolkits.basemap import Basemap
import scipy.stats as sts
import matplotlib.pyplot as plt
from scipy import optimize

### Read in data reader function 
import sitReader as R

### Read in sea ice data (lat x lon are equal grids)
lat,lon,sitpi = R.piomasReader('icej') # Piomas on icesat domain
lat,lon,sitpc = R.piomasReader('cryo') # Piomas on cryosat domain
lat,lon,siti,sitc = R.icesatReader()   # siti = ICESat, sitc = CryoSat-2

siti = siti[0,:,:]
sitpi = sitpi[0,:,:]

i = np.ravel(siti)
p = np.ravel(sitpi)

#i[np.where(np.isnan(i))] = 0.
#p[np.where(np.isnan(p))] = 0.
#
#i1 = np.ma.masked_array(i,mask=i==np.nan)
#p1 = np.ma.masked_array(p,mask=i==np.nan)

#mask = np.isfinite(i) & np.isfinite(p)
#i = i[mask]
#p = p[mask]

numbers = np.arange(0,len(i),1)

### Create test model
basic_model = pm.Model()
with basic_model:

    ### priors for unkown model parameters
    alpha = pm.Normal('alpha',mu=2,sd=2,testval=0.1)
    beta = pm.Normal('beta',mu=2,sd=2,shape=2,testval=0.1)
    sigma = pm.Uniform('sigma')
    
    ### expected value of outcome
    mu = alpha + beta[0]*p + beta[1]*i
    
    ### likelihood (sampling distribution) of observations
    yobs = pm.Normal('yobs',mu=mu,sd=sigma,observed=i)
 
    start = pm.Metropolis()   
#    step = pm.NUTS(state=start)
    trace = pm.sample(2000,start,progressbar=True)
    pm.traceplot(trace)
    pm.summary(trace)
#
#a = trace['alpha']
#b = trace['beta']
#b = b[:,0]
#x = np.linspace(0,7.01,2000)
#
#plt.figure()
#plt.plot(x,a*x+b)
#
#plt.figure()
#m = 2
#b1 = 1
#line = m*x + b
#plt.plot(x,line)