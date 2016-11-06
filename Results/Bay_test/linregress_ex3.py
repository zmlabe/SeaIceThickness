"""
Source : https://scrogster.wordpress.com/2011/04/05/pymc-for-bayesian-models/

Date : 11 August 2016
"""

### Imports 
import numpy as np
import matplotlib.pyplot as plt
import pymc as bae

yy = np.array([-19.23776, 1.559197, 27.90364, -14.94222, -41.34614, 5.857922,
               -26.24492, -1.670176, -8.349098, -24.91511, 63.86167, 20.87778, 
               4.176622, -35.65956, 4.482383, 36.20763, 33.60314, 23.25372,
               -15.52639, -25.59295, 42.48803, -29.46465, 30.25402, -5.66534, 
               -20.92914, 44.87109, 19.07603, 22.19699, 18.89613, 2.835296, 
               12.68109, -17.19655, 26.60962, -28.74333, -24.69688,  -19.02279,
               -31.39471, -17.83819, 15.389, 40.41935, 0.972758, -36.49488,
               -2.041068, 23.22597, 1.226252, 11.87125, 36.32597, 29.20536,
               16.24043, -0.8978296])
xx = np.array([-14, -6, 19, -12, -16, 1, -15, -13, 0, -6, 15, 8, 1, -16, -5,
               19, 8, 7, -11, -13, 13, -18, 10, -1, -13, 13, 13, 17, 13, 11,
               4, -6, 14, -14, 3, -3, -18, -11, 6, 13, -10, -12, -2, 9, -7,
               -1, 14, 15, 6, -2])

plt.scatter(xx,yy)

### priors
sigma = bae.Uniform('sigma',0.0,200.0,value=20)
alpha = bae.Normal('alpha',0.0,0.001,value=0)
beta = bae.Normal('beta',0.0,0.001,value=0)

### model
@bae.deterministic(plot=False)
def modelled_yy(XX=xx,beta=beta,alpha=alpha):
    return beta*XX + alpha

### likelihood
y=bae.Normal('y',mu=modelled_yy,tau=1.0/sigma**2,value=yy,observed=True)

### mcmc sample
m = bae.MCMC()
m.sample(10000,burn=5000)
bae.Matplot.plot(m)