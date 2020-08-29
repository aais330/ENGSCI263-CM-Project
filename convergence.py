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
# Create array of stepsizes ranging from 0.1 to 4
stepsizes=np.arange(0.1,4.1,0.01)
# Pre-allocating arrays 
sols=[]
h=[]

# Loop through each step size and solve for concentration
for j in range(len(stepsizes)):
    tv = np.arange(1998,2020,step=stepsizes[j])
    P = solve_dPdt(dPdt,tv,pi,[b])
    C = solve_dCdt(dCdt,tv,P,b1,alpha,bc,tau)
    for i in range(len(tv)):
        # Find solution at 2010 if there is one for each step size
        if tv[i]==atyear:
            # Store in 2010 solutions array
            sols.append(C[i])
            # Store 1/h in corresponding stepsize array
            h.append(1/stepsizes[j])

# Plot concentration at 2010 vs 1/stepsize
f,ax = plt.subplots(1,1)
ax.plot(h, sols,'mo--', label = 'Nitrate concentration in 2010 using Improved Eulers') 
ax.legend() 
ax.set_title('Convergence of nitrate concentration step sizes')
ax.set_xlabel('Step size (1/h)')
ax.set_ylabel('Rate of change in concentration in 2010')
plt.show()
# plt.savefig('convergence.png')

