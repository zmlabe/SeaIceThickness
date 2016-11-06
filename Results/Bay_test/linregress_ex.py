"""
Source : https://users.obs.carnegiescience.edu/cburns/ipynbs/PyMC.html

Date : 11 August 2016
"""

### Imports 
import numpy as np
import matplotlib.pyplot as plt
import pymc as bae

### Create fake data
### Linear regression model with two predictor valariables (x,y) and 
### one outcome variable (z).
### Relationship between them is [equation]

Nobs = 20 
x_true = np.random.uniform(0,10,size=Nobs)
y_true = np.random.uniform(-1,1,size=Nobs)

### Values for relationship [equation]
alpha_true = 0.5
beta_x_true = 1.0
beta_y_true = 10.0
eps_true = 0.5
z_true = alpha_true + (beta_x_true*beta_x_true) + (beta_y_true*y_true)
z_obs = z_true + np.random.normal(0,eps_true,size=Nobs)

####### build model in hierachal model
### add parameters
alpha = bae.Uniform('alpha',-100,100,value=np.median(z_obs))
betax = bae.Uniform('betax',-100,100,value=np.std(z_obs)/np.std(x_true))
betay = bae.Uniform('betay',-100,100,value=np.std(z_obs)/np.std(y_true))
eps = bae.Uniform('eps',0,100,value=0.01)

### define model
### Deterministic decorator tells pymc that model is a function
### that depends on the stochastic (random variables) objects, but is itself
### a determistic function of them. That is --> a function of random
### variables not a random variable itself. 
@bae.deterministic
def model(alpha=alpha,betax=betax,betay=betay,x=x_true,y=y_true):
    return alpha + betax*x + betay*y
    
### pymc parameterizes the width of the normal distribution by
### tau=1/sigma**2
### tau converts standard deviation parameter to precisison parameter
@bae.deterministic
def tau(eps=eps):
    return np.power(eps,-2)
    
### lastly relate the model/parameters to the data
### z_obs drawn from normal distribution with mean equation to the model,
### special observed=True parameter tells pymc that this stochastic's
### value should remain constant
data = bae.Normal('data',mu=model,tau=tau,value=z_obs,observed=True)

### now lets feed data into pymc.MCMC
sampler = bae.MCMC([alpha,betax,betay,eps,model,tau,z_obs,x_true,y_true])
sampler.use_step_method(bae.AdaptiveMetropolis,[alpha,betax,betay,eps],
                       scales={alpha:0.1,betax:0.1,betay:1.0,eps:0.1})
sampler.sample(iter=10000)

### plot
#bae.Matplot.plot(sampler)

alpha.summary()

### take median of each variable from sampler
m_alpha = np.median(alpha.trace())
m_betax = np.median(betax.trace())
m_betay = np.median(betay.trace())
m_eps = np.median(eps.trace())