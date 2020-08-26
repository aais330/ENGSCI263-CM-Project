import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


def improved_euler_step(f, tk, xk, h, pars = []): 

    ''' Computes a single Improved Euler step.

	Paramters
	-----------
	f: callable
		Derivate function.
	xk: float
		Independent variable value at begining of step.
	yk: float
		Solution at begining of step.
	h: float
		step size.
	pars: iterable
		Optional parameters to pass to derivate function.

	Returns
	-----------
	xk1: float
		Solution at end of Improved Euler step

	Notes
	---------
	Assumes the order of inputs to f is f(y,t,*pars).
    '''
    xk1_prediction = xk + h*f(xk,tk,*pars)
    return xk + (h/2)*(f(xk,tk,*pars) + f(xk1_prediction, tk+h, *pars))

# Pressure ODE
def dPdt(P, t, b):
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
    dP_a = 0.5 # Pressure difference across aquifer
    t_mar = 2020 # Time when MAR begins

    if (t>t_mar):
        dP_a += 0.5 # Pressure difference increase due to MAR
        
    return -b*(P + dP_a/2) -b*(P-(dP_a)/2) 

# Concentration ODE
def dCdt(ci, t, P, b1, alpha, bc, tau):
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
    
    if ((t-tau) > 1990.5): # THINK ABOUT THIS!
        ni = np.interp((t-tau),tn,n) #interpolating number of cows
    else:
        ni = n[0]
    
    
    # Active carbon sink
    if ((t-tau)>t_acs):
        b1 = alpha*b1

    # MAR
    if (t>t_mar):
        dP_a += 0.5 # Pressure difference increase due to MAR
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))


def solve_dPdt(f, t, pi, b):
    '''
    Parameters:
    -----------
    f : callable
        Function that returns dxdt given variable and parameter inputs.
    t : array-like
        array of floats containing time value.
    P0 : float
        Initial value of solution.
     b : float
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
        P[i+1] = improved_euler_step(f,t[i],P[i],dt,b)
        
    return P

# Solve concentration ODE
def solve_dCdt(f,t,P, b1, alpha, bc, tau):
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
    C[0] = 0.2
    dt = t[1]-t[0]

    # Solve using improved euler method
    for i in range(len(t)-1):
        pars = [P[i],b1,alpha,bc,tau]
        output = improved_euler_step(f,t[i],C[i],dt,pars)
        C[i+1] = output

    return C


'''
LPM_Model is a single function that solves the LPM for nitrate concentration in the aquifer
'''

def LMP_Model(t, b ,b1, alpha, bc,tau):
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

    Returns
    -------
    P : array-like  
        array of pressure values in aquifer
    C : array-like
        array of concentration values in the aquifer  
    t : array-like
        array of times that correspond to pressure and concentration values
    '''
    # intial pressure
    pi = 0
      
    # Solve pressure ODE
    P = solve_dPdt(dPdt,t,pi,[b])

    # Solve concentration ODE 
    C = solve_dCdt(dCdt,t,P,b1,alpha,bc,tau)

    return C


# Solve pressure ODE
pi = 0
t = np.arange(1980,2020,step = 0.1)
#parameters
b=0.5
b1=0.5
alpha=0
bc=1
tau = 5
# LMP_Model(t,b,b1,alpha,bc,tau)


# Testing curve_fit
# load in cow data and concentration data
tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T
tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

pars = curve_fit(LMP_Model,tcon,c,[1,1,1,1,15])
# print(pars)
b=pars[0][0]
b1=pars[0][1]
alpha=pars[0][2]
bc=pars[0][3]
tau = pars[0][4]
C = LMP_Model(t,b,b1,alpha,bc,tau)


f,ax = plt.subplots(1,1)

ax.plot(t,C,'k', label = 'Numeric Solution')
ax.plot(tcon,c,'r+', label = 'Data')
ax.set_title('Numerical Solution and data')
plt.ylabel('Concentration')
plt.xlabel('Time')

'''
# Plotting cows for reference
ax2 = ax.twinx()
ax2.plot(tn,n,'b', label = 'Data cows')
ax2.set_ylabel('Number of cows')
ax.legend() 
'''
    
# plt.show()
