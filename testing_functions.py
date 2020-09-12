import os
import numpy as np
from model_functions import *
tol=1.e-8

# Unit testing functions
def func(x,y,pars=[]):
    pass

def test_ie():
    #simple function testing
    pass

def test_dPdt():
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

    # test when switches to tmar on
    # THINK THIS IS TESTING IMPROVED EULER NOT dPdt
    # # t < t_mar
    # check_dPdt4 = improved_euler_step(dPdt, 2019, -6, 1, pars=[0.2])
    # assert((check_dPdt4+4.08)<tol) #calculation gave me 11.975, maybe i missed something?

    # check_dPdt5 = improved_euler_step(dPdt, 1998, 6, 0.2, pars=[0.2])
    # assert((check_dPdt5-5.54)<tol) # Keep getting 5.54 even though the year changes

    # # t > t_mar
    # check_dPdt6 = improved_euler_step(dPdt, 2021, 6, 0.2, pars=[0.2])
    # assert((check_dPdt6-5.55)<tol)


    print("dPdt passed \n")


# dCdt(ci, t, P, b1, alpha, bc, tau):

def test_dCdt():

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

    # Testing Improved euler not dCdt
    # # Tests when t-tau > 1990.5 
    # # Test when t-tau>t_acs 2015>2010
    # check_dCdt5=improved_euler_step(dCdt,2020,0.2,1,pars=[0,b1,alpha,bc,tau])
    # assert((check_dCdt5-0.1665)<tol)

    # # Test when t-tau<t_acs 2009<2010
    # check_dCdt6=improved_euler_step(dCdt,2014,0.2,1,pars=[0,b1,alpha,bc,tau])
    # assert((check_dCdt6-13870.7)<tol)

    # # Test when t>t_mar
    # check_dCdt7=improved_euler_step(dCdt,2021,0.2,1,pars=[0,b1,alpha,bc,tau]) 
    # assert((check_dCdt7-0.149)<tol) 

    # # Test when t<t_mar
    # check_dCdt8=improved_euler_step(dCdt,2018,0.2,1,pars=[0,b1,alpha,bc,tau])
    # assert((check_dCdt8-0.1903)<tol)

    print("dCdt passed \n")


# convergence testing solutions for different step sizes
def convergence_analysis():
    '''
    Convergence analysis for our model
    '''

    # pre-allocating solutions array
    sols = []
    step_sizes = np.arange(0.1, 1, 0.1)

    pars = posterior_pars()[1]

    for step in 1/step_sizes:
        # time range with step size
        t = np.arange(1998, 2020, step)

        C = LPM_Model(t, *pars)
        C = np.interp(2012, t, C)
        # add concentration solution at each step size to array
        sols.append(C)
    

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(step_sizes, sols, 'mo--', label = 'Nitrate concentration in 2012')
    ax2.legend() 
    ax2.set_title('Convergence of nitrate concentration step sizes')
    ax2.set_xlabel('Step size (1/h)')
    ax2.set_ylabel('Nitrate Concentration in 2012')
    #plt.show()
    fig2.savefig('Plots'+ os.sep + 'convergence.png')
    plt.close(fig2)


# Benchmarking functions
# Theses functions benchmark the two ODEs formulated in the Modified_Derivates script

# Simplfied pressure ODE
def dPdt_simplified(P, t):
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

    return -1*(P + dP_a/2) -1*(P-(dP_a1)/2) 

# Simplified concentration ODE for -100 cows
def dCdt_simplified(ci, t, P,b1,alpha, bc, tau):
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

    ni = -100
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))

# Simplified concentration ODE for 100 cows
def dCdt_simplified1(ci, t, P,b1,alpha, bc, tau):
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

    ni = 100
             
    return  -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))

# Pressure benchmark without MAR
def pressure_benchmark():
    # Analytical Solution
    t = np.arange(1980,2020,0.5)
    P_analytical = np.zeros(t.shape)

    # Numerical Solution
    P_numerical = solve_dPdt(dPdt_simplified,t)

    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)

    ax.plot(t,P_analytical,'b', label = 'Analytical Solution')
    ax.plot(t,P_numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution')
    plt.ylabel('Pressure')
    plt.xlabel('Time')
    plt.xlim(1980,2020)
    plt.ylim(-5,5)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    #plt.show()
    plt.savefig('Plots'+ os.sep +"pressure_benchmark.png")
    plt.close(fig)

# Concentration benchmark with negative 100 cows
def concentration_benchmark1():    
    t = np.arange(1980,2020,0.5)
    C_Analytical = (100/np.exp(-49.5))*np.exp(-0.025*t)-100            # FIX THIS

    
    b1=0.5
    alpha=0
    bc=0.5
    tau = 0
    
    P_numerical = solve_dPdt(dPdt_simplified,t)
    C_Numerical = solve_dCdt(dCdt_simplified,t,P_numerical,b1,alpha,bc,tau)

    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)

    ax.plot(t,C_Analytical,'b', label = 'Analytical Solution')
    ax.plot(t,C_Numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (negative 100 cows)')
    plt.ylabel('Concentration')
    plt.xlabel('Time')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    
    #plt.show()
    plt.savefig('Plots'+ os.sep +"concentration_benchmark1.png")
    plt.close(fig)

# Constant number of cows (100)
def concentration_benchmark2(): 
    t = np.arange(1980,2020,0.5)
    C_Analytical = (-499/5)*np.exp((-t+1980)/20) + 100 

    
    b1=1
    alpha=0
    bc=1
    tau = 0
    
    P_numerical = solve_dPdt(dPdt_simplified,t)
    C_Numerical = solve_dCdt(dCdt_simplified1,t,P_numerical,b1,alpha,bc,tau)

    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)

    ax.plot(t,C_Analytical,'b', label = 'Analytical Solution')
    ax.plot(t,C_Numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (100 cows)')
    plt.ylabel('Concentration')
    plt.xlabel('Time')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    
    #plt.show()
    plt.savefig('Plots'+ os.sep +"concentration_benchmark2.png")
    plt.close(fig)
