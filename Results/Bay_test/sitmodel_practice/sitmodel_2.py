"""
Practice bayesian modeling for March sit data.

check this example : https://healthyalgorithms.com/2010/08/27/mcmc-in-
                     python-global-temperature-reconstruction-with-pymc/
 
Author : Zachary Labe
Date : 22 August 2016
"""

### Read in modules
import numpy as np
import matplotlib.pyplot as plt
import pymc as bae 
from mpl_toolkits.basemap import Basemap

### Read in data reader function 
import sitReader as R

### Read in sea ice data (lat x lon are equal grids)
lat,lon,sitpi = R.piomasReader('icej') # Piomas on icesat domain
lat,lon,sitpc = R.piomasReader('cryo') # Piomas on cryosat domain
lat,lon,siti,sitc = R.icesatReader()   # siti = ICESat, sitc = CryoSat-2

### Create model
# priors
beta = bae.Normal('beta',mu=np.zeros(3),tau=.1,value=np.zeros(3))
sigma = bae.Uniform('sigma',lower=0,upper=10.,value=0.01)

# predictions
@bae.deterministic
def mu(beta=beta,model=sitpi,satellite1=siti):
    return beta[0] + beta[1]*satellite1 + beta[2]*model**2  
    
@bae.deterministic
def predicted(mu=mu,sigma=sigma):
    return bae.rnormal(mu,sigma**-2)
    
# likelihood
@bae.observed
def y(value=sitpi,mu=mu,sigma=sigma):
    return bae.normal_like(value,mu,sigma**-2)
    
### MCMC
var = [beta,sigma,mu,predicted,y]
m = bae.MCMC(var)
m.use_step_method(bae.Metropolis,beta)
m.sample(iter=20000,thin=10,burn=10000,verbose=1)

quantiles = predicted.stats()['quantiles']
q90l = quantiles[2.5]
smooth = quantiles[50]
q90h = quantiles[97.5]

### Define figure
fig = plt.figure()
ax = plt.subplot(111)

style='ortho'

if style == 'ortho':
    mm = Basemap(projection='ortho',lon_0=-90,
                lat_0=70,resolution='l',round=True)
elif style == 'polar':
    mm = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    mm.drawparallels(parallels,labels=[False,False,False,False],
                    linewidth=0.5,color='k',fontsize=9)
    mm.drawmeridians(meridians,labels=[True,True,False,False],
                    linewidth=0.5,color='k',fontsize=9)
    
mm.drawmapboundary(fill_color='white')
mm.drawcoastlines(color='k',linewidth=0.5)
mm.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

info = 'Sea Ice Thickness'

### Plot filled contours    
cs = mm.contourf(lon[:,:],lat[:,:],smooth[0],
                latlon=True,extend='max')
                  
### Set colormap                              
cs.set_cmap('viridis')

### Set colorbar
cbar = mm.colorbar(cs,drawedges=True,location='right',pad = 0.55)

