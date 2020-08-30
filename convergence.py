from model_functions import *

# Set parameters of the models
b=0.5
b1=0.5
alpha=0
bc=1
tau = 5
pi = 0
# Choose the concentrations in year 2010
atyear = 2010
# Create array of stepsizes ranging from 0.1 to 3
stepsizes=np.arange(0.1,3,0.15)
# Pre-allocating solutions array
sols=[]

# Loop through each step size and solve for concentration
for step in 1/stepsizes:
    # create t range with step sixe
    tv = np.arange(1998,2020,step)
    # calculate pressure
    P = solve_dPdt(dPdt,tv,pi,[b])
    # calculate concentration
    C = solve_dCdt(dCdt,tv,P,b1,alpha,bc,tau)
    # interpolate concentration at year 2010
    C = np.interp(atyear,tv,C)
    # add concentration to 2010 solutions array
    sols.append(C)

# Plot concentration at 2010 vs 1/stepsize
f,ax = plt.subplots(1,1)
ax.plot(stepsizes, sols,'mo--', label = 'Nitrate concentration in 2010 using Improved Eulers') 
ax.legend() 
ax.set_title('Convergence of nitrate concentration step sizes')
ax.set_xlabel('Step size (1/h)')
ax.set_ylabel('Rate of change in concentration in 2010')
plt.show()
plt.savefig('convergence.png')

