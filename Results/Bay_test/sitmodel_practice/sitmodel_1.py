"""
Practice bayesian modeling for March sit data.

check this example : https://healthyalgorithms.com/2010/08/27/mcmc-in-
                     python-global-temperature-reconstruction-with-pymc/
 
Author : Zachary Labe
Date : 22 August 2016
"""

import numpy as np
import matplotlib.pyplot as plt
import pymc as bae
import sitReader as R
from mpl_toolkits.basemap import Basemap

lat,lon,sitpi = R.piomasReader('icej')
lat,lon,sitpc = R.piomasReader('cryo')
lat,lon,siti,sitc = R.icesatReader()

pi = np.ravel(sitpi)
si = np.ravel(siti)
#mask = np.isfinite(pi) & np.isfinite(si)
#pi = pi[mask]
#si = si[mask]

### Attempt 1

# define priors
beta = bae.Normal('beta',mu=np.zeros(2),tau=0.001,value=np.zeros(2))
sigma = bae.Uniform('sigma',lower=0.,upper=1.5,value=0.1)

# define predictions 
@bae.deterministic
def mu(beta=beta,sat=si):
    return beta[0]+beta[1]*pi

@bae.deterministic
def predicted(mu=mu,sigma=sigma):
    return bae.rnormal(mu,sigma**-2.)
    
# define likelihood
def y(value=pi,mu=mu,sigma=sigma):
    bae.normal_like(value,mu,sigma**-2.)
    
var = [beta,sigma,mu,predicted,y]
m = bae.MCMC(var)
m.use_step_method(bae.Metropolis,beta)
m.sample(iter=20000,thin=10,burn=10000,verbose=1)

quantiles = predicted.stats()['quantiles']
q90l = quantiles[2.5]
smooth = quantiles[50]
q90h = quantiles[97.5]

q90l = np.reshape(q90l,(sitpi.shape))
q90h = np.reshape(q90h,(sitpi.shape))
smooth = np.reshape(smooth,(sitpi.shape))

smooth=smooth

smooth[np.where(smooth<0.)]=np.nan

plt.figure()
plt.plot(q90l[0],color='grey',linestyle='--')
plt.plot(q90h[0],color='grey',linestyle='--')
plt.plot(smooth[0],color='r',linestyle='-')

### Define figure
fig = plt.figure()
ax = plt.subplot(111)

style='polar'

if style == 'ortho':
    mm = Basemap(projection='ortho',lon_0=-90,
                lat_0=70,resolution='l',round=True)
elif style == 'polar':
    mm = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
    
mm.drawmapboundary(fill_color='white')
mm.drawcoastlines(color='k',linewidth=0.5)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
mm.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.5,color='k',fontsize=9)
mm.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.5,color='k',fontsize=9)
mm.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

info = 'Sea Ice Thickness'

### Plot filled contours    
cs = mm.contourf(lon[:,:],lat[:,:],q90h[0],
                latlon=True,extend='max')
                  
### Set colormap                              
cs.set_cmap('viridis')

### Set colorbar
cbar = mm.colorbar(cs,drawedges=True,location='right',pad = 0.55)