import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from model_functions_posterior import *


# Pressure ODE
def dPdt_forecast(P, t, b, dP_Mar):
    '''
    Parameters:
    -----------
    P : float
        pressure value (dependent variable)
    t : float
        time value (independent variable)
    b : float
        Recharge strength parameter 
        
    Returns:
    -------
     dPdt : float
        rate of change of pressure in aquifer
    '''
    
    dP_a = 0.1 # Pressure difference across aquifer(given in project doc)
    dP_a1 = 0.1
    

    
    t_mar = 2020 # Time when MAR begins

    if (t>=t_mar): # FIX THIS IF ELSE STATEMENT
        dP_a1  += dP_Mar


        
    return -b*(P + dP_a/2) -b*(P-dP_a1/2)

# Concentration ODE
def dCdt_forecast(ci, t, P, b1, alpha, bc, tau, dP_Mar):
    '''
    Parameters
    ----------
    ci : float
        Concentration of nitrate (dependent variable)
    t : float
        time value (independent variable)
    P : float
        Current pressure of aquifer
    b1 : float
        infliltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        dilution parameter
    tau : float
        time lag paramter
    Returns
    -------
    dCdt : float
        rate of change of concentration in aquifer
    '''
    dP_a = 0.1 # Pressure difference across aquifer
    dP_surf = 0.05 # Oversurface pressure
    t_mar = 2020 # Time when MAR begins
    t_acs = 2010 # Time active carbon sink was installed


    # number of cows
    tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T
    
    if ((t-tau) >= 1990.5): # THINK ABOUT THIS!
        ni = np.interp((t-tau),tn,n) #interpolating number of cows
    else:
        ni = 10000
    
    
    # Active carbon sink
    if ((t-tau)>t_acs):
        b1 = alpha*b1

    # MAR
    if (t>t_mar):
        dP_a += dP_Mar # Pressure difference increase due to MAR
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))


def solve_dPdt_forecast(f, t, pi, b, dP_Mar):
    '''
    Parameters:
    -----------
    f  : callable
        Function that returns dxdt given variable and parameter inputs.
    t  : array-like
        array of floats containing time value.
    pi : float
        Initial value of solution.
    b  : float
        recharge strength parameter
        
    Returns: 
    --------
    P : array like
        array of floats containing pressure in aquifer
    '''
    P = np.zeros(t.shape) # intialising pressure vector
    P[0] = pi 
    dt = t[1] - t[0]

    # Solve using Improved Euler Method
    for i in range(len(t)-1):
        P[i+1] = improved_euler_step(f,t[i],P[i],dt,[b, dP_Mar])
        
    return P

# Solve concentration ODE
def solve_dCdt_forecast(f,t,P,ci, b1, alpha, bc, tau, dP_Mar):
    '''
    Parameters:
    -----------
    f : callable
        function that returns dxdt given variable and parameter inputs.
    t : array-like
        array of floats containing time value.
    P : array-like
        array of floats containing the pressure within the aquifer
    b1 : float
        infliltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        dilution parameter
    tau : float
        time lag paramter
    Returns:
    -------
    C : array-like
        An array of concentrations
    '''
    C = np.zeros(t.shape) # intialising concentration array
    C[0] = ci
    dt = t[1]-t[0]

    # Solve using improved euler method
    for i in range(len(t)-1):
        pars = [P[i],b1,alpha,bc,tau, dP_Mar]
        output = improved_euler_step(f,C[i],t[i],dt,pars)
        C[i+1] = output

    return C

'''
LPM_Model is a single function that solves the LPM for nitrate concentration in the aquifer
'''

def LPM_Model_forecast(t, ci, b, pi ,b1,alpha, bc,tau, dP_Mar):
    '''
    Parameters
    ----------
    t : array-like
        Array of time values to solve LPM
    b : float
        Recharge strength parameter
    b1 : float
        infliltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
            dilution parameter
    tau : float
        time lag
    dP_Mar: float
        Mar induced pressure difference
    Returns
    -------
    P : array-like  
        array of pressure values in aquifer
    C : array-like
        array of concentration values in the aquifer  
    t : array-like
        array of times that correspond to pressure and concentration values
    '''
    
    
    #Define tv
    tv = np.arange(1980,2030,step = 0.25)
    
    #mean_dP_Mar = dP_Mar
    #sd_dP_Mar = 0.2*dP_Mar
    #dP_Mar = sd_dP_Mar * np.random.randn() + mean_dP_Mar

    P = solve_dPdt_forecast(dPdt_forecast,tv,pi,b, dP_Mar)

    # Solve concentration ODE 
    C = solve_dCdt_forecast(dCdt_forecast,tv,P,ci,b1,alpha,bc,tau, dP_Mar)

    # INTERPOLATE to T (from Tcon)
    C_interp = np.interp(t,tv,C)
    
    return C_interp




# Testing curve_fit
# load in cow data and concentration data
tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
t = np.arange(1980,2030,step = 0.25)
'''
pars = curve_fit(LPM_Model_forecast,tcon,c,[0.2,0.5,0.5,1,1,5])
print(pars)
b=pars[0][0]
b1=pars[0][1]
alpha=pars[0][2]
bc=pars[0][3]
'''
tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
#ci = c[0]

dP_Mar = 0.01

pos, p = posterior_pars()
print(2)
ci = p[0]
b = p[1]
pi = p[2]
b1 = p[3]
alpha = p[4]
bc = p[5]
tau = p[6]


v=0.3




C_Out = LPM_Model_forecast(t,ci,b,pi,b1,alpha, bc,tau, dP_Mar)
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.plot(t, C_Out, 'b-', label='best-fit')

#ax.errorbar(tcon,c,yerr=v,fmt='ro', label='data', markersize=2.2)
plt.show()