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



'''
LPM_Model is a single function that solves the LPM for nitrate concentration in the aquifer
'''

def LMP_Model(t,pi, b ,b1, alpha, bc, tau):
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

    # Reading in data
    # load in cow data and concentration data
    tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T
    tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

    #Pressure ODE
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
        
        dPdt = -b*(P + dP_a/2) -b*(P-(dP_a)/2) 
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
    P = solve_dPdt(dPdt,t,pi,[b])
    
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
        t_acs = 2020 # Time active carbon sink was installed


        # number of cows
        tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T
        '''
        if ((t-tau) > 1990.5):
            ni = np.interp(t-tau,tn,n) #interpolating number of cows
        else:
            ni = n[0]
        '''
        ni = 0

        # Active carbon sink
        if (t>t_acs):
            b1 = alpha*b1

        # MAR
        if (t>t_mar):
            dP_a += 0.5 # Pressure difference increase due to MAR

        change = -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))
        return  change

    # Solve concentration ODE
    def solve_dCdt(f,t,P, b1, alpha, bc, tau):
        '''
        Parameters:
        -----------
        f : callable
            Function that returns dxdt given variable and parameter inputs.
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

    # Solve concentration ODE
    C = solve_dCdt(dCdt,t,P,b1,alpha,bc,tau)
    plt.plot(t,C)
    plt.show()


    return 0


# Solve pressure ODE
pi = 0
t = np.arange(1980,2020,step = 0.5)
#parameters
b=0.5
b1=0.5
alpha=0
bc=1
tau = 5
LMP_Model(t,pi,b,b1,alpha,bc,tau)