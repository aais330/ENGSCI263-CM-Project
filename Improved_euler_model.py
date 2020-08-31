import numpy as np
from matplotlib import pyplot as plt
import math
from concmodel import *

def ode_model(t,P,b,dP_aq,dP_Mar=0):
    ''' Return the derivative dx/dt at time, t, for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable.
        P : float
            Dependent variable.
       
        b : float
            Recharge strength parameter.
        dP_aq : float
            Pressure difference of acquifer
        dP_Mar: optional float (if not put in, defaults to 0)
            MAR pressure added

        Returns:
        --------
        dPdt : float
            Derivative of dependent variable with respect to independent variable.

       

    '''
    
    DPdt = -b * (P + dP_aq/2) -b * (P-(dP_aq + dP_Mar)/2) 
    
    return DPdt
    

    


def euler_step(f, tk, yk, dt, pars):
	''' Compute a single Euler step.
	
		Parameters
		----------
		f : callable
			Derivative function.
		tk : float
			Independent variable at beginning of step.
		yk : float
			Solution at beginning of step.
		dt : float
			Step size.
		pars : iterable
			Optional parameters to pass to derivative function.
			
		Returns
		-------
		yk1 : float
			Solution at end of the Euler step.
	'''
	# use *pars to pass an array of parameters into a function
	# e.g. f(xk, yk, *pars)
	dPdt = f(tk, yk, *pars)         #compute dy/dx at (xk, yk)
	y1 = yk + dt * dPdt              #Compute single euler step using stepsize h
	return y1


def improved_euler_step(f, tk, yk, dt, pars):
	'''
	Performs 1 step of improved euler method

	Parameters:
	---------
	f: callable
		The ODE for which a solution must be computed
	tk: float
		Independent variable value at beginning of step
	yk: float
		Solution value at beginning of step
	dt: float
		Step size to take
	pars: iterable
		Additional parameters to pass to function apart from t and y

	Returns:
	---------
	y2 : float
		Solution at end of the improved-euler step.

	Notes:
	-------
	f must be such that order of inputs is f(x,y,*pars).

	'''

	y1 = euler_step(f, tk, yk, dt, pars) #Calling the euler_step, to compute the predictor step
	x1 = tk + dt #Increase time by h 
	y2 = yk + (dt/2) * (f(tk, yk, *pars) + f(x1, y1, *pars)) #Combining predictor and corrector steps to get y2 (the output)

	return y2


def solve_ode(f, t0, t1, dt, P0, t_Mar, dP_Mar, pars):
    ''' Solve an ODE numerically.

        Parameters:
        -----------
        f : callable
            Function that returns dxdt given variable and parameter inputs.
        t0 : float
            Initial time of solution.
        t1 : float
            Final time of solution.
        dt : float
            Time step length.
        P0 : float
            Initial value of solution.
        t_Mar: int
            Year at which to start MAR
        dP_Mar: float
            Added pressure due to Mar
        pars : array-like
            List of parameters passed to ODE function f.
            These include: 
            b, dP_aq, dP_Mar (in order) 


        Returns:
        --------
        xs : array-like
            Independent variable solution vector.
        ys : array-like
            Dependent variable solution vector.
        '''
    nx = int(np.ceil((t1-t0)/dt))  # compute number of improved-Euler steps to take
    xs = t0 + np.arange(nx+1)*dt
    ys = 0.*xs
    ys[0] = P0  #First element of y array = y0
    
    
    
    for i in range(nx):
        if (xs[i] >= t_Mar):
            if (len(pars) != 3):
                pars.append(dP_Mar)
        ys[i+1] = improved_euler_step(f,xs[i], ys[i], dt, pars) #perform nx improved euler steps, and store them in ys
    

    return xs, ys

def plot_benchmark():

    t, P = solve_ode(f=ode_model, t0=1980, t1=2030, dt=1, P0=0, t_Mar=2020, dP_Mar=100000, pars=[0.5, 0.1*10**6])

    #pars are b and dP_aq for above

    fig, ax = plt.subplots()


    

    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('Pressure (Pa)')
    ax.legend('Analytic solution')
    
    

    ax.plot(t,P, 'b')
    ax.legend('Improved Euler solution')
    
    # plt.show()

if __name__ == "__main__":
    plot_benchmark()

    #order of calling functions:

    #plot_benchmark() --> solve_ode() --> improved_euler_step --> euler_step --> ode_model()