#import this file to any other file to access the functions

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
    
    dP_a = 0.1 # Pressure difference across aquifer(given in project doc)
    dP_a1 = 0.1
    dP_mar = 0.05

    t_mar = 2030 # Time when MAR begins

    if (t>t_mar): # FIX THIS IF ELSE STATEMENT
        dP_a1  += dP_mar


        
    return -b*(P + dP_a/2) -b*(P-(dP_a1)/2) 

# Concentration ODE
def dCdt(ci, t, P,b1,alpha, bc, tau):
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
    t_mar = 2030 # Time when MAR begins
    t_acs = 2010 # Time active carbon sink was installed
    

    # number of cows
    tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T #move outside
    '''
    if ((t-tau) >= 1990.5): # THINK ABOUT THIS!
        ni = np.interp((t-tau),tn,n) #interpolating number of cows
    else:
        ni = 10000
   
    if ((t-tau) <= tn[0]):
        ni = n[0]
    elif ((t-tau)>=tn[-1]):
        ni = n[-1]
    else:
    '''
    
    ni = np.interp((t-tau),tn,n)
    
    # Active carbon sink
    if ((t-tau)>t_acs):
        b1 = alpha*b1

    # MAR
    if (t>t_mar):
        dP_a += 0.05 # Pressure difference increase due to MAR
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))


def solve_dPdt(f, t, b):
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
    #P[0] = pi 
    dt = t[1] - t[0]

    # Solve using Improved Euler Method
    for i in range(len(t)-1):
        P[i+1] = improved_euler_step(f,t[i],P[i],dt,b)
        
    return P

# Solve concentration ODE
def solve_dCdt(f,t,P, b1,alpha, bc, tau):
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
    #C[0] = 0.2  OLD INTIAL CONCENTRAION (DO NOT DELETE)
    #C[0] = 4.2
    dt = t[1]-t[0]

    tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T

    # Solve using improved euler method
    for i in range(len(t)-1):
        pars = [P[i],b1,alpha,bc,tau]
        output = improved_euler_step(f,t[i],C[i],dt,pars)
        C[i+1] = output

    return C


'''
LPM_Model is a single function that solves the LPM for nitrate concentration in the aquifer
'''

def LPM_Model(t, b,b1,alpha, bc,tau):
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
    
    #Define tv
    # tv = np.arange(1980,2020,step = 0.1) ORIGNAL TIME PERIOD (DO NOT DELETE)
    tv = np.arange(1980,2030,step = 0.25) # New Time period
      
    # Solve pressure ODE
    P = solve_dPdt(dPdt,tv,[b])

    # Solve concentration ODE 
    C = solve_dCdt(dCdt,tv,P,b1,alpha,bc,tau)

    # INTERPOLATE to T (from Tcon)
    C_interp = np.interp(t,tv,C)
    return C_interp

def posterior_pars(sigma):
    '''
    Parameter
    ---------
    sigma : array
        Variance limit of pars
    Returns
    -------
    pos : ndarray
        Array of shape (N,) containing spread of parameter values
    p : array
        Optimal values of the parameters with minimal variance from data
    
    Notes
    -----
    Utilises other functions in the same file and requires no inputs
    '''
    # reading data for curve_fit calibration
    t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

    #sigma = [0.1]*len(c0) # variance limit of pars


    sigma = [0.1]*len(c0) # variance limit of pars

    # calibrating model to data and creating covariance matrix
    p, cov = curve_fit(LPM_Model,t0,c0,bounds=((0,0,0,0,0),(np.inf,np.inf,1,np.inf,5))) 

    pos = np.random.multivariate_normal(p, cov, 100) # random variates of the calibrated pars
    
    return pos, p