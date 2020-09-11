# This script benchmarks the two ODEs formulated in the Modified_Derivates script

from model_functions import *

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

    f,ax = plt.subplots(1,1)

    ax.plot(t,P_analytical,'b', label = 'Analytical Solution')
    ax.plot(t,P_numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution')
    plt.ylabel('Pressure')
    plt.xlabel('Time')
    plt.xlim(1980,2020)
    plt.ylim(-5,5)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    plt.show()

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

    f,ax = plt.subplots(1,1)

    ax.plot(t,C_Analytical,'b', label = 'Analytical Solution')
    ax.plot(t,C_Numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (negative 100 cows)')
    plt.ylabel('Concentration')
    plt.xlabel('Time')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    
    plt.show()

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

    f,ax = plt.subplots(1,1)

    ax.plot(t,C_Analytical,'b', label = 'Analytical Solution')
    ax.plot(t,C_Numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (100 cows)')
    plt.ylabel('Concentration')
    plt.xlabel('Time')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    
    plt.show()






