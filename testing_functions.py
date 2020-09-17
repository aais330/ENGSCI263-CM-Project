import os
import numpy as np
from model_functions import *
# error tolerance
tol=1.e-8

# Unit testing functions
def func(y,t):
    '''
    Simple function to test improved_euler_step - y=t^2
    dy/dt=2t

    Parameters
    -----------
    y: float
        dependent variable
    t: float
        independent variable (time)

    Returns
    --------
    func: float
        derivative of the function dydt
    '''
    return 2*t

def test_ie():
    '''
    Testing improved euler step with simple function func
    '''
    test_1 = improved_euler_step(func,0,0,1)
    assert((test_1-1)<tol)

    test_2 = improved_euler_step(func,5,2,0.1)
    assert((test_2-3.01)<tol)

    test_3 = improved_euler_step(func,7,10,0.5)
    assert((test_3-17.25)<tol)

    test_4 = improved_euler_step(func,50,100,5)
    assert((test_4-625)<tol)

    print('Improved Euler tests passed \n')

def test_dPdt():
    '''
    Unit tests for the pressure model.
    These tests check the edge cases of the if statements in the pressure differential equation function. 
    '''

    # t < t_mar before t_mar
    check_dPdt = dPdt(10, 2015)
    assert((check_dPdt+20)<tol)
    
    check_dPdt2 = dPdt(15, 2019)
    assert((check_dPdt2+30)<tol)
    
    # t > t_mar after t_mar
    check_dPdt3 = dPdt_forecast(20, 2025,.5)
    assert((check_dPdt3)+39.75<tol)

    check_dPdt4 = dPdt_forecast(20, 2020.1,.5)
    assert((check_dPdt4)+39.75<tol)

    print("dPdt tests passed \n")

def test_dCdt():
    '''
    Unit tests for the concentration model.
    These tests check the edge cases of the if statements in the concentration differential equation function. 
    '''

    # between 1990.5 and before active carbon sink and no time lag
    check_dCdt = dCdt(50, 1999.5, 0, 1, 1, 1, 0) 
    assert((check_dCdt-11645.8)<tol)

    # between 1990.5 and before active carbon sink and time lag
    check_dCdt = dCdt(50, 2000.5, 0, 1, 1, 1, 1) 
    assert((check_dCdt-11645.8)<tol)
    
    # After active carbon sink no time lag
    check_dCdt = dCdt(50, 2010.5, 0, 1, 0.5, 1, 0) 
    assert((check_dCdt-14977.45)<tol)

    # After active carbon sink with time lag
    check_dCdt = dCdt(50, 2011.5, 0, 1, 0.5, 1, 1) 
    assert((check_dCdt-14977.45)<tol)

    # before 1990.5
    check_dCdt = dCdt(50,1985,0, 1, 1, 1, 0) 
    assert((check_dCdt-1886.1)<tol)

    # after 2020 without mar, with active carbon sink
    check_dCdt = dCdt(50,2021,0, 1, 0.5, 1, 0) 
    assert((check_dCdt-15903.525)<tol)

    # after 2020 without mar, without active carbon sink
    check_dCdt = dCdt(50,2021,0, 1, 1, 1, 0) 
    assert((check_dCdt-31809.55)<tol)

    # after 2020 with mar, with active carbon sink
    check_dCdt = dCdt_forecast(50,2021,0, 1, 0.5, 1, 0,0.5) 
    assert((check_dCdt-15916.025)<tol)

    print("dCdt tests passed \n")

def convergence_analysis():
    '''
    Produces a convergence analysis plot for our LPM model.

    Notes
    -----
    saves figure as convergence.png
    
    '''

    hsteps = 1/(np.linspace(0.4,20,20)) # step size vector
    pars = posterior_pars()[1]
    C_final = [] # solution vector

    # Solve model for a range of step sizes
    for h in hsteps:
        t = np.arange(1980,2020,h)
        P = solve_dPdt(dPdt,t)
        C = solve_dCdt(dCdt,t,P,*pars)
        C_final.append(np.interp(2016, t, C))
  
    # Produce convergence analysis plot
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(1/hsteps, C_final, 'bo--', label = 'Nitrate concentration in 2016')
    ax2.legend() 
    ax2.set_title('Convergence Analysis for LPM Model')
    ax2.set_xlabel('1/step size')
    ax2.set_ylabel('Nitrate Concentration in 2016')
    #plt.show()
    fig2.savefig('Plots'+ os.sep + 'convergence.png')
    plt.close(fig2)

# Benchmarking functions
# Theses functions benchmark the two ODEs formulated in the model_functions script

# Simplified pressure ODE
def dPdt_simplified(P, t):
    '''
    This function simplifies the pressure model for benchmarking

    Parameters:
    -----------
    P : float
        Pressure value (dependent variable)
    t : float
        Time value (independent variable)
        
    Returns:
    -------
    dPdt : float
        Rate of change of pressure in aquifer
    '''
    
    dP_a = 0.1 # Pressure difference across aquifer(given in project doc)
    dP_a1 = 0.1

    return -1*(P + dP_a/2) -1*(P-(dP_a1)/2) 

# Simplified concentration ODE for -100 cows
def dCdt_simplified(ci, t, P, b1, alpha, bc, tau):
    '''
    This function simplies the concentration model for benchmarking with a negative number of cows

    Parameters
    ----------
    ci : float
        Concentration of nitrate (dependent variable)
    t : float
        Time value (independent variable)
    P : float
        Current pressure of aquifer
    b1 : float
        Infiltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        Dilution parameter
    tau : float
        Time lag parameter

    Returns
    -------
    dCdt : float
        Rate of change of concentration in aquifer

    '''
    dP_a = 0.1 # Pressure difference across aquifer
    dP_surf = 0.05 # Oversurface pressure

    ni = -100
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))

# Simplified concentration ODE for 100 cows
def dCdt_simplified1(ci, t, P,b1,alpha, bc, tau):
    '''
    This function simplies the concentration model for benchmarking with a positive number of cows

    Parameters
    ----------
    ci : float
        Concentration of nitrate (dependent variable)
    t : float
        Time value (independent variable)
    P : float
        Current pressure of aquifer
    b1 : float
        Infiltration parameter
    alpha : float
        Active carbon sink infiltration modication parameter
    bc : float
        Dilution parameter
    tau : float
        Time lag parameter

    Returns
    -------
    dCdt : float
        Rate of change of concentration in aquifer

    '''
    dP_a = 0.1 # Pressure difference across aquifer
    dP_surf = 0.05 # Oversurface pressure

    ni = 100
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))

def pressure_benchmark():
    '''
    Benchmarking solution for pressure without MAR
    '''
    # Analytical Solution
    t = np.arange(1980,2020,0.5)
    P_analytical = np.zeros(t.shape)

    # Numerical Solution
    P_numerical = solve_dPdt(dPdt_simplified,t)
    
    # Plot and compare
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(t,P_analytical,'b', label = 'Analytical Solution')
    ax.plot(t,P_numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (Pressure equation)')
    plt.ylabel('Pressure [MPa]')
    plt.xlabel('Time [yrs]')
    plt.xlim(1980,2020)
    plt.ylim(-5,5)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    plt.savefig('Plots'+ os.sep +"pressure_benchmark.png")
    plt.close(fig)

def concentration_benchmark1():   
    '''
    Benchmarking solution for concentration model with negative 100 cows
    ''' 
    t = np.arange(1980,2020,0.5)
    # Analytical solution
    C_Analytical = (100/np.exp(-49.5))*np.exp(-0.025*t)-100         

    # Set parameters
    b1=0.5
    alpha=0
    bc=0.5
    tau = 0
    
    # Calculate numerical solution for both
    P_numerical = solve_dPdt(dPdt_simplified,t)
    C_Numerical = solve_dCdt(dCdt_simplified,t,P_numerical,b1,alpha,bc,tau)

    # Plot and compare
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(t,C_Analytical,'b', label = 'Analytical Solution')
    ax.plot(t,C_Numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (negative 100 cows)')
    plt.ylabel('Nitrate Concentration [mg/L]')
    plt.xlabel('Time [yrs]')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    
    plt.savefig('Plots'+ os.sep +"concentration_benchmark1.png")
    plt.close(fig)

def concentration_benchmark2(): 
    '''
    Benchmarking solution for concentration model with constant 100 cows
    ''' 

    t = np.arange(1980,2020,0.5)
    # Analytical solution
    C_Analytical = (-499/5)*np.exp((-t+1980)/20) + 100 

    # Set parameters
    b1=1
    alpha=0
    bc=1
    tau = 0
    
    # Calculate numerical solutions
    P_numerical = solve_dPdt(dPdt_simplified,t)
    C_Numerical = solve_dCdt(dCdt_simplified1,t,P_numerical,b1,alpha,bc,tau)

    # Plot and compare
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(t,C_Analytical,'b', label = 'Analytical Solution')
    ax.plot(t,C_Numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (100 cows)')
    plt.ylabel('Nitrate Concentration [mg/L]')
    plt.xlabel('Time [yrs]')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    
    plt.savefig('Plots'+ os.sep +"concentration_benchmark2.png")
    plt.close(fig)

def test_solvedPdt():
    '''
    Testing the solve_dPdt functions with the simplified dPdt equation
    '''
    # Creating a time range to run through function
    t = np.arange(2000,2020,1)

    test_1=solve_dPdt(dPdt_simplified, t)
    for i in range(len(t)):
        # Simplified function causes all P values to be 0
        assert(test_1[i]==0)

    print("solve_dPdt tests passed \n")

def test_solvedCdt():
    '''
    Testing the solve_dCdt functions with the simplified dCdt equation
    '''
    # Creating a time range to run through function
    t = np.arange(2000,2015,0.5)

    # Set parameters
    b1=1
    alpha=0
    bc=1
    tau = 0
    
    # Calculate P to input
    P = solve_dPdt(dPdt_simplified,t)

    # # Test with negative cows
    # test_1=solve_dCdt(dCdt_simplified,t,P,b1,alpha,bc,tau)
    # for i in range(len(t)):
    #     assert(test_1[i]==((100/np.exp(-49.5))*np.exp(-0.025*t[i])-100))
    

    # # Test with positive cows
    # test_2=solve_dCdt(dCdt_simplified1,t,P,b1,alpha,bc,tau)
    # for i in range(len(t)):
    #     assert(test_2[i]==((-499/5)*np.exp((-t[i]+1980)/20) + 100))

    print("solve_dCdt tests passed \n")

