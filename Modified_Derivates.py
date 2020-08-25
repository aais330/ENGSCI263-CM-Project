import numpy as np
from Improved_euler_model import *

'''
LPM_Model is a single function that solves the LPM for nitrate concentration in the aquifer
'''

def LMP_Model(t,pi, b,b1,alpha,bc):
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

    Returns
    -------
    P : array-like  
        array of pressure values in aquifer
    C : array-like
        array of concentration values in the aquifer  
    t : array-like
        array of times that correspond to pressure and concentration values
    '''

    # Reading in data
    # load in cow data and concentration data
    tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T
    tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

    #Pressure ODE
    def dpPdt(pi, t, b):
        '''
        Parameters:
        -----------
        pi : float
            intial pressure value (dependent variable)
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
            dP_mar = 0.5 # Pressure difference increase due to MAR
        else:
            dP_mar = 0
        
        dPdt = -b*(pi + dP_a/2) -b*(pi-(dP_a + dP_mar)/2) 
        return dPdt

   
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

    # Solve pressure ODE
    P = solve_dPdt(dpPdt,t,pi,[b])

    print(P)

    # Concentration ODE
    def dCdt(ci, t, b1, alpha, bc):
        '''
        Parameters
        ----------
        ci : float
            Concentration of nitrate (dependent variable)
        t : float
            time value (independent variable)
        b1 : float
            infliltration parameter
        alpha : float
            Active carbon sink infiltration modication parameter
        bc : float
            dilution parameter

        Returns
        -------
        dCdt : float
            rate of change of concentration in aquifer
        '''
        ni = np.interp(t,tn,n) # interpolate number of cows

    
        return  0 #-ni*b1*(P-Psurf)+bc*c*(P-(Pa/2))



    return 0


# Solve pressure ODE
pi = 0
t = np.arange(1980,2020,step = 0.5)
#parameters
b=0.5
b1=0
alpha=0
bc=0
LMP_Model(t,pi,b,b1,alpha,bc)