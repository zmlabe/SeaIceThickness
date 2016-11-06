# run_mcmc.py

import pymc as bay
import simple_normal_model
import numpy as np
import matplotlib.pyplot as plt

### build MCMC model from our previous normal distribution model
model = bay.MCMC(simple_normal_model)

### run the model 500x
model.sample(iter=500)

### print stats about model
print '\n \n',model.stats()

#print model.trace('precision')[:]

mean = model.trace('mean')[:]
std = model.trace('std_dev')[:]

plt.plot(mean,label='mean')
plt.plot(std,label='std')
plt.legend()
plt.xlim([-5,5])