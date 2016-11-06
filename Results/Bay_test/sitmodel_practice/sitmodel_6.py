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

#pi = np.ravel(sitpi)
#si = np.ravel(siti)
pi = sitpi
si = siti

### Attempt 1

# define priors
alpha = bae.Normal('alpha',mu=np.nanmean(si),tau=0.1,
                   observed=False)
beta = bae.Normal('beta',mu=np.nanmean(pi),tau=0.1,value=np.zeros(1),
                  observed=False)
sigma = bae.Uniform('sigma',lower=0.,upper=1.5,value=0.1)

tau_c = bae.Gamma('tau_c',alpha=1.0,beta=1.0,value=1)

# define predictions 
@bae.deterministic
def mu(alpha=alpha,beta=beta,sat=si):
    return alpha + beta[0]*pi

@bae.stochastic
def mu_phi(tau=tau_c,value=np.zeros((si.shape))):
    mu = np.empty((si.shape))
    for yr in xrange(si.shape[0]):
        for i in xrange(si.shape[1]):
            for j in xrange(si.shape[2]):
                if i == 0:
                    mu[yr,i,j] = pi[yr,i,j]
                elif i == si.shape[1]-1:
                    mu[yr,i,j] = pi[yr,i,j]
                elif j == 0:
                    mu[yr,i,j] = pi[yr,i,j]
                elif j == si.shape[2]-1:
                    mu[yr,i,j] = pi[yr,i,j]
                else:
                    mu[yr,i,j] = (pi[yr,i,j] + pi[yr,i-1,j] \
                    + pi[yr,i,j-1] + pi[yr,i+1,j] + pi[yr,i,j+1])/5.
    taux = tau*5.
    return bae.normal_like(pi,mu,taux)
    
phi = bae.Lambda('phi',lambda mu=mu_phi: mu-np.nanmean(mu))

@bae.deterministic
def predicted(mu=mu,sigma=sigma):
    return bae.rnormal(mu,sigma**-2.)
    
# define likelihood
def y(value=pi,mu=mu,sigma=sigma):
    bae.normal_like(value,mu,sigma**-2.)
    
var = [alpha,beta,sigma,mu,mu_phi,phi,predicted,y]
m = bae.MCMC(var)
m.use_step_method(bae.Metropolis,beta)
m.sample(iter=2000,thin=10,burn=300,verbose=1)

bae.Matplot.plot(m)
plt.savefig('sitmodel6_tests.png',dpi=300)

quantiles = predicted.stats()['quantiles']
q90l = quantiles[2.5]
smooth = quantiles[50]
q90h = quantiles[97.5]

q90l = np.reshape(q90l,(sitpi.shape))
q90h = np.reshape(q90h,(sitpi.shape))
smooth = np.reshape(smooth,(sitpi.shape))
smooth[np.where(smooth<0.)]=np.nan

plt.figure()
plt.plot(q90l[0],color='grey',linestyle='--')
plt.plot(q90h[0],color='grey',linestyle='--')
plt.plot(smooth[0],color='r',linestyle='-')
plt.savefig('sitmodel_6_testss.png',dpi=300)

### Define figure
fig = plt.figure()
ax = plt.subplot(111)

#sit = q90h[0]
diff = smooth[-1]-smooth[0]
sit = smooth[0]
print '\n change in SIT', round(np.nanmean(diff),2),'m'
print '\n mean SIT', round(np.nanmean(sit),2),'m'
print '\n max SIT', round(np.nanmax(sit),2),'m'
print '\n min SIT', round(np.nanmin(sit),2),'m'

style='polar'

if style == 'ortho':
    mm = Basemap(projection='ortho',lon_0=-90,
                lat_0=70,resolution='l',round=True)
elif style == 'polar':
    mm = Basemap(projection='npstere',boundinglat=66,lon_0=270,
                 resolution='l',round =True)
    
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
cs = mm.contourf(lon[:,:],lat[:,:],sit,
                latlon=True,extend='max')
                  
### Set colormap                              
cs.set_cmap('viridis')

### Set colorbar
cbar = mm.colorbar(cs,drawedges=True,location='right',pad = 0.55)
plt.savefig('sitmodel_6_test.png',dpi=300)