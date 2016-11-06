"""
Practice bayesian modeling for March sit data using pymc3

 
Author : Zachary Labe
Date : 24 August 2016
"""

### Read in modules -----> using new PyMC3
import numpy as np
import matplotlib.pyplot as plt
from pymc3 import *
from mpl_toolkits.basemap import Basemap
import scipy.stats as sts

### Read in data reader function 
import sitReader as R

### Read in sea ice data (lat x lon are equal grids)
lat,lon,sitpi = R.piomasReader('icej') # Piomas on icesat domain
lat,lon,sitpc = R.piomasReader('cryo') # Piomas on cryosat domain
lat,lon,siti,sitc = R.icesatReader()   # siti = ICESat, sitc = CryoSat-2

###########################################################################
###########################################################################
###########################################################################
### Estimates posterior predictive in regression model
i = np.ravel(siti)
p = np.ravel(sitpi)
timex = np.arange(0,8,1)

mask = np.isfinite(i) & np.isfinite(p)
slope, intercept, r_value, p_value, std_err = sts.linregress(p[mask],i[mask])
line = slope*timex + intercept

data = dict(x=p,y=i)
with Model() as model:
    glm.glm('y ~ x',data)
    start = find_MAP()
    step = NUTS(scaling = start)
    trace = sample(300,step,progressbar=True)

plt.figure()   
plt.plot(timex,timex,color='k',linewidth=2,zorder=1) 
plt.plot(p,i,'x')
glm.plot_posterior_predictive(trace,samples=300,eval=p,label='posterior predictive')
plt.plot(timex,line,color='r',linewidth=2,label='True lin regress')
plt.legend(loc='upper left',fontsize=6)
plt.ylim([0,7])
plt.xlim([0,7])
plt.xlabel('PIOMAS')
plt.ylabel('ICESat')
plt.savefig('posterior_predict.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
### Is it normally distributed?

with Model() as model_robusts:
    family = glm.families.StudentT()
    glm.glm('y~x',data,family=family)
    start = find_MAP()
    step = NUTS(scaling=start)
    trace_robust = sample(2000,step,progressbar=True)
    
plt.figure()
plt.plot(timex,timex,color='k',linewidth=2,zorder=1)
plt.plot(p,i,'x')
glm.plot_posterior_predictive(trace_robust,eval=p,label='posterior predictive')
plt.plot(timex,line,color='r',linewidth=2,label='True lin regress')
plt.legend(loc='upper left',fontsize=6)
plt.ylim([0,7])
plt.xlim([0,7])
plt.xlabel('PIOMAS')
plt.ylabel('ICESat')
plt.savefig('posterior_predict_tdist.png',dpi=300)
