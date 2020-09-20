# Implementation of LPM model of Edendale aquifer

import os
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from scipy import stats

# Cow and nitrate concentration data
tn, n = np.genfromtxt("Data"+ os.sep +'nl_cows.txt', delimiter=',', skip_header=1).T
t0, c0 = np.genfromtxt("Data"+ os.sep +'nl_n.csv', delimiter=',', skip_header=1).T


def improved_euler_step(f, tk, xk, h, pars = []): 
    ''' 
    Computes a single Improved Euler step.

	Parameters
	-----------
	f: callable
		Derivate function.
	tk: float
		Independent variable value at beginning of step.
	xk: float
		Solution at beginning of step.
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


def dPdt(P, t):
    '''
    Pressure ODE, rate of change of pressure within the aquifer

    Parameters:
    -----------
    P : float
        pressure value (dependent variable)
    t : float
        time value (independent variable)
        
    Returns:
    -------
     dPdt : float
        rate of change of pressure in aquifer
    '''
    
    dP_a = 0.1 # Pressure difference across aquifer(given in project doc)
    dP_a1 = 0.1
    dP_mar = 0.0

    t_mar = 2030 # Time when MAR begins

    if (t>t_mar): 
        dP_a1  += dP_mar # Pressure increase at high pressure boundary due to MAR


        
    return -1*(P + dP_a/2) -1*(P-(dP_a1)/2) 


def dCdt(ci, t, P, b1, alpha, bc, tau):
    '''
    Concentration ODE, rate of change of nitrate concentration in aquifer

    Parameters
    ----------
    ci : float
        Concentration of nitrate (dependent variable)
    t : float
        time value (independent variable)
    P : float
        Current pressure of aquifer
    b1 : float
        infiltration parameter
    alpha : float
        Active carbon sink infiltration modification parameter
    bc : float
        dilution parameter
    tau : float
        time lag parameter

    Returns
    -------
    dCdt : float
        rate of change of concentration in aquifer
    '''
    dP_a = 0.1 # Pressure difference across aquifer
    dP_surf = 0.05 # Oversurface pressure
    t_mar = 2030 # Time when MAR begins
    t_acs = 2010 # Time active carbon sink was installed
    
    ni = np.interp((t-tau),tn,n) # Interpolating number of cows
    
    # Active carbon sink
    if ((t-tau)>t_acs):
        b1 = alpha*b1 

    # MAR
    if (t>t_mar):
        dP_a += 0.0 # Pressure increase at high pressure boundary due to MAR (for increase see dCdt_forecast)
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))


def solve_dPdt(f, t):
    '''
    Solves the ODE for pressure within the aquifer

    Parameters:
    -----------
    f  : callable
        Function that returns dPdt given variable and parameter inputs.
    t  : array-like
        array of floats containing time value.
        
    Returns: 
    --------
    P : array-like
        array of floats containing pressure in aquifer
    '''
    P = np.zeros(t.shape) # intialising pressure vector
    dt = t[1] - t[0]

    # Solve using Improved Euler Method
    for i in range(len(t)-1):
        P[i+1] = improved_euler_step(f,t[i],P[i],dt)
        
    return P


def solve_dCdt(f,t,P, b1,alpha, bc, tau):
    '''
    Solves the ODE for nitrate concentration within the aquifer

    Parameters:
    -----------
    f : callable
        function that returns dCdt given variable and parameter inputs.
    t : array-like
        array of floats containing time value.
    P : array-like
        array of floats containing the pressure within the aquifer
    b1 : float
        infiltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        dilution parameter
    tau : float
        time lag parameter

    Returns:
    -------
    C : array-like
        An array of concentrations
    '''
    C = np.zeros(t.shape) # intialising concentration array
    dt = t[1]-t[0]

    # Solve using improved euler method
    for i in range(len(t)-1):
        pars = [P[i],b1,alpha,bc,tau]
        output = improved_euler_step(f,t[i],C[i],dt,pars)
        C[i+1] = output

    return C


def LPM_Model(t, b1, alpha, bc, tau):
    '''
    Solves the pressure then concentration ODEs in serial

    Parameters
    ----------
    t : array-like
        Array of time values to solve LPM
    b1 : float
        Infiltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        Dilution parameter
    tau : float
        Time lag

    Returns
    -------
    C_interp : array-like
        Array of concentration values in the aquifer corresponding to times in input array t
    '''
    

    tv = np.arange(1980,2030,step = 0.1) # time array to solve for concentration
      
    # Solve pressure ODE
    P = solve_dPdt(dPdt,tv)

    # Solve concentration ODE 
    C = solve_dCdt(dCdt,tv,P,b1,alpha,bc,tau)

    # interpolate concentration at times corresponding to the times in input array t
    C_interp = np.interp(t,tv,C)
    return C_interp


def posterior_pars():
    '''
    Finds the parameter values of the new best fit model

    Returns
    -------
    pos : ndarray
        Array of shape (N,) containing alternative of parameter values
    p : array
        Optimal values of the parameters with minimal variance from data
    
    Notes
    -----
    parameter values are found using built-in curve_fit function
    '''

    sigma = [0.15]*len(c0) # uncertainty in data

    # calibrating model to data and creating covariance matrix
    p, cov = curve_fit(LPM_Model,t0,c0, sigma = sigma,bounds=((0,0,0,0.5),(1e-04,np.inf,np.inf,np.inf)), absolute_sigma=True)
    
    pos = np.random.multivariate_normal(p, cov, 100) # random variates of the calibrated parameters
    
    return pos, p


# Functions for forecasting 
def dPdt_forecast(P, t, dP_Mar):
    '''
    Pressure ODE, rate of change of pressure within the aquifer with managed aquifer recharge

    Parameters:
    -----------
    P : float
        pressure value (dependent variable)
    t : float
        time value (independent variable)
    dP_Mar : float
        Pressure increase due to managed aquifer recharge
        
    Returns:
    -------
     dPdt : float
        rate of change of pressure in aquifer
    '''
    
    dP_a = 0.1 # Pressure difference across aquifer(given in project doc)
    dP_a1 = 0.1
    
    t_mar = 2020 # Time when MAR begins

    if (t>=t_mar):
        dP_a1  += dP_Mar # Pressure increase at high pressure boundary due to MAR
        
    return -1*(P + dP_a/2) -1*(P-dP_a1/2)


def dCdt_forecast(ci, t, P, b1, alpha, bc, tau, dP_Mar):
    '''
    Concentration ODE, rate of change of nitrate concentration in aquifer

    Parameters
    ----------
    ci : float
        Concentration of nitrate (dependent variable)
    t : float
        time value (independent variable)
    P : float
        Current pressure of aquifer
    b1 : float
        infiltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        dilution parameter
    tau : float
        time lag parameter
    dP_Mar : float
        Pressure increase due to managed aquifer recharge

    Returns
    -------
    dCdt : float
        rate of change of concentration in aquifer
    '''
    dP_a = 0.1 # Pressure difference across aquifer
    dP_surf = 0.05 # Oversurface pressure
    t_mar = 2020 # Time when MAR begins
    t_acs = 2010 # Time active carbon sink was installed
    
    ni = np.interp((t-tau),tn,n) # interpolating number of cows
  
    # Active carbon sink
    if ((t-tau)>t_acs):
        b1 = alpha*b1

    # MAR
    if (t>t_mar):
        dP_a += dP_Mar # Pressure difference increase high pressure due to MAR
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))


def solve_dPdt_forecast(f, t, dP_Mar):
    '''
    Solves the ODE for pressure within the aquifer with managed aquifer recharge

    Parameters:
    -----------
    f  : callable
        Function that returns dxdt given variable and parameter inputs.
    t  : array-like
        array of floats containing time value.
    pi : float
        Initial value of solution.
    dP_Mar  : float
        Pressure increase due to managed aquifer recharge
        
    Returns: 
    --------
    P : array-like
        array of floats containing pressure in aquifer
    '''
    P = np.zeros(t.shape) # intialising pressure vector
    dt = t[1] - t[0]

    # Solve using Improved Euler Method
    for i in range(len(t)-1):
        P[i+1] = improved_euler_step(f,t[i],P[i],dt,[dP_Mar])
        
    return P


def solve_dCdt_forecast(f,t,P, b1, alpha, bc, tau, dP_Mar):
    '''
    Solves the ODE for nitrate concentration within the aquifer with managed aquifer recharge 

    Parameters:
    -----------
    f : callable
        function that returns dxdt given variable and parameter inputs.
    t : array-like
        array of floats containing time value.
    P : array-like
        array of floats containing the pressure within the aquifer
    b1 : float
        infiltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        dilution parameter
    tau : float
        time lag parameter
    dP_Mar  : float
        Pressure increase due to managed aquifer recharges
    Returns:
    -------
    C : array-like
        An array of concentrations
    '''
    C = np.zeros(t.shape) # intialising concentration array
    dt = t[1]-t[0]

    # Solve using improved euler method
    for i in range(len(t)-1):
        pars = [P[i],b1,alpha,bc,tau, dP_Mar]
        output = improved_euler_step(f,t[i],C[i],dt,pars)
        C[i+1] = output

    return C


def LPM_Model_forecast(t,b1,alpha, bc,tau, dP_Mar):
    '''
    Parameters
    ----------
    t : array-like
        Array of time values to solve LPM
    b1 : float
        infiltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        dilution parameter
    tau : float
        time lag
    dP_Mar: float
        MAR induced pressure difference

    Returns
    -------
    C : array-like
        array of concentration values in the aquifer corresponding to times in input array t 
    '''
    

    tv = np.arange(1980,2030,step = 0.1) # time array of to solve for concentration

    P = solve_dPdt_forecast(dPdt_forecast,tv,dP_Mar) 
    C = solve_dCdt_forecast(dCdt_forecast,tv,P,b1,alpha,bc,tau, dP_Mar)

    # interpolate concentration at times corresponding to the times in input array t
    C_interp = np.interp(t,tv,C)
    
    return C_interp

def confidence_int(data, message): 
    '''
    Computes and prints a 90% confidence interval for the data

    Parameters:
    ----------
    data: array-like
        Data set for which to find  90% confidence interval for
    message: string
        string that describes confidence interval
        

    Returns:
    ----------
    ci: 2 element tuple
        lower and upper bounds of confidence interval
        

    Notes:
    ------
        Prints confidence interval to screen with message
    '''

    std = np.std(data) #mean and standard deviation
    mean = np.mean(data)
    ci = stats.norm.interval(0.9,loc = mean, scale = std)
    print("The 90 percent confidence interval for " + message, "is: ", ci) #printing confidence interval
    
    return ci
