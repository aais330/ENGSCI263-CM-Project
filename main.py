######################################################################
# This scipt generates all the related figures for the nitrate
# leaching project Team 12:
# Anton Aish, Jenny Choi, William Laughton, Annie Li, Sidharth Varma
######################################################################
import numpy as np
import matplotlib
from model_functions import *
from testing_functions import *

# Getting posterior samples and best fit parameters
pos, pars = posterior_pars()

# Plotting the data
data_cow = np.genfromtxt("nl_cows.txt", delimiter = ',', skip_header = 1)
x_cow = data_cow[:,0]
y_cow = data_cow[:,1]

data_nconc = np.genfromtxt("nl_n.csv", delimiter = ',', skip_header = 1)
x_nconc = data_nconc[:,0]
y_nconc = data_nconc[:,1]

fig1 = plt.figure(figsize = (20,10))
ax = fig1.add_subplot(111)
cows = ax.plot(x_cow, y_cow, color = 'red', label = 'Cows')
ax.set_xlabel('Years', fontsize = 20)
ax.set_ylabel('Number of Cows', fontsize = 20)
ax2 = ax.twinx()
n_conc = ax2.plot(x_nconc, y_nconc, color = 'blue', label = 'Nitrate Concentration')
ax2.set_ylabel('Nitrate Concentration(mg/L)', fontsize = 20)
plt.title('Annual Cattle Numbers in Southland and Nitrate Concentrations of Edendale Farm', fontsize = 20)
ax.plot([],[], 'b', label='Nitrate Concentration')
ax.legend()
fig1.savefig('nc_data.png', dpi = 200)

# unit tests
test_ie()
test_dPdt()
test_dCdt()


# Convergence Analysis
step_sizes = np.arange(0.05, 0.1, 0.05)
step_solutions = convergence_sols(step_sizes, pars)
fig2 = plt.figure(figsize= (20,10))
ax2 = fig2.add_subplot(111)
ax2.plot(step_sizes, step_solutions, 'mo--', label = 'Nitrate concentration in 2012')
ax2.legend() 
ax2.set_title('Convergence of nitrate concentration step sizes')
ax2.set_xlabel('Step size (1/h)')
ax2.set_ylabel('Nitrate Concentration in 2012')
fig2.savefig('convergence.png')



# Benchmarks



# Intial model



# Improved model




# Use forecasting and without active carbon sink




# Forecasts with uncertainty

