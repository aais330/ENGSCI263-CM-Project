# This script benchmarks the two ODEs formulated in the Modified_Derivates script

from Modified_Derivates import *

# Pressure benchmark without MAR
if False:
    # Analytical Solution
    t = np.arange(1980,2020,0.5)
    P_analytical = np.zeros(t.shape)

    # Numerical Solution
    P_numerical = solve_dPdt(dPdt,t,0,[1])

    f,ax = plt.subplots(1,1)

    ax.plot(t,P_analytical,'b', label = 'Analytical Solution')
    ax.plot(t,P_numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution')
    plt.ylabel('Pressure')
    plt.xlabel('Time')
    plt.xlim(1980,2020)
    plt.ylim(-5,5)

    plt.show()

# Concentration benchmark with no cows
# NOTE: to run this test you must go into the concentration derivative function (dCdt) and 
# set the number of cows to 0.
if False:    
    t = np.arange(1980,2020,0.5)
    C_Analytical = np.exp((-0.05*t + 99))*0.2

    b=0.5
    b1=0.5
    alpha=0
    bc=1
    tau = 5
    
    #P_numerical = solve_dPdt(dPdt,t,0,[1])
    #C_Numerical = solve_dCdt(dCdt,t,P_numerical,b1,alpha,bc,tau)
    C_Numerical = LMP_Model(t,b,b1,alpha,bc,tau)

    f,ax = plt.subplots(1,1)

    ax.plot(t,C_Analytical,'b', label = 'Analytical Solution')
    ax.plot(t,C_Numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (no cows)')
    plt.ylabel('Concentration')
    plt.xlabel('Time')
    
    plt.show()

# Constant number of cows (100)
# NOTE: to run this test you must go into the concentration derivative function (dCdt) and 
# set the number of cows to 100
if False:
    t = np.arange(1980,2020,0.5)
    C_Analytical = (-499/5)*np.exp((-t+1980)/20) + 100 

    b=1
    b1=1
    alpha=0
    bc=1
    tau = 0
    
    P_numerical = solve_dPdt(dPdt,t,0,[1])
    C_Numerical = solve_dCdt(dCdt,t,P_numerical,b1,alpha,bc,tau)

    f,ax = plt.subplots(1,1)

    ax.plot(t,C_Analytical,'b', label = 'Analytical Solution')
    ax.plot(t,C_Numerical,'r+', label = 'Numeric Solution')

    ax.set_title('Analytical vs Numeric solution (100 cows)')
    plt.ylabel('Concentration')
    plt.xlabel('Time')
    
    plt.show()






