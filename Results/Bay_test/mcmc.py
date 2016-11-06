from pymc import Metropolis
import disaster_model
from pymc import MCMC
M = MCMC(disaster_model)
M.use_step_method(Metropolis, disaster_model.late_mean, proposal_sd=2.)