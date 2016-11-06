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

###########################################################################
###########################################################################
###########################################################################
### hierarchical model example

with pm.Model() as hierarchical_model:
    
    ### a = PIOMAS, b = ICESat
    
    mu_a = pm.Normal('mu_alpha',mu=np.nanmean(sitpi),sd=1**2)
    sigma_a = pm.Uniform('sigma_alpha',lower=0,upper=2)
    mu_b = pm.Normal('mu_beta',mu=np.nanmean(siti),sd=1*2)
    sigma_b = pm.Uniform('sigma_beta',lower=0,upper=5)
    
    a = pm.Normal('alpha',mu=mu_a,sd=sigma_a,shape=i.shape[0])
    b = pm.Normal('beta',mu=mu_b,sd=sigma_b,shape=i.shape[0])
    
    eps = pm.Uniform('eps',lower=0,upper=5)
    
    estimate = a[0] + b[0]*p
    
    est_like = pm.Normal('est_like',mu=estimate,sd=eps,observed=p)
    
with hierarchical_model:
#    start = pm.find_MAP()
#    step = pm.NUTS(scaling=start)
    mu,sds,elbo = pm.variational.advi(n=100)
    step = pm.NUTS(scaling=sds.get('beta')**2,
                   is_cov=True)
    sit_trace = pm.sample(5000,step,start=start,progressbar=True)