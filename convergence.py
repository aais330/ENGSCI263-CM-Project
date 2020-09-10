from model_functions_posterior import *

# Set parameters of the models
pos, p = posterior_pars()


# Choose the concentrations in year 2012
atyear = 2012
# Create array of stepsizes ranging from 0.1 to 3
stepsizes=np.arange(0.05,1,0.05)
# Pre-allocating solutions array
sols=[]

# Loop through each step size and solve for concentration
for step in 1/stepsizes:
    # create t range with step sixe
    tv = np.arange(1998,2020,step)

    C = LPM_Model(tv,*p)
    C = np.interp(atyear,tv,C)
    # add concentration to 2010 solutions array
    sols.append(C)

# Plot concentration at 2012 vs 1/stepsize
f,ax = plt.subplots(1,1)
ax.plot(stepsizes, sols,'mo--', label = 'Nitrate concentration in 2012 using Improved Eulers') 
ax.legend() 
ax.set_title('Convergence of nitrate concentration step sizes')
ax.set_xlabel('Step size (1/h)')
ax.set_ylabel('Nitrate Concentration in 2012')
plt.show()
#plt.savefig('Convergence.png')

