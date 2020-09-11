######################################################################
# This scipt generates all the related figures for the nitrate
# leaching project Team 12:
# Anton Aish, Jenny Choi, William Laughton, Annie Li, Sidharth Varma
######################################################################
import numpy as np
import matplotlib
from model_functions import *
from testing_functions import *
from benchmark import *
from plotting_functions import *
# Getting posterior samples and best fit parameters
pos, pars = posterior_pars()

# Plotting the data
plot_data()

# unit tests
test_ie() 
test_dPdt()
test_dCdt()


# Convergence Analysis
convergence_analysis()


# Benchmarks
# Pressure benchmark
pressure_benchmark()

# Concentration benchmark
concentration_benchmark1()
concentration_benchmark2()

# Initial model
initial_model()

# Improved model
improved_model()


#forecasting different mar pressures and without active carbon sink
what_ifs()
without_acs()


#Forecasts with uncertainty
what_ifs_uncertainty()
without_acs_uncertainty()
