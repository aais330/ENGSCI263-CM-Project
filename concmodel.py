# imports
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

def conc_model(t,c,tc,n,b,P,Psurf,bc,Pa,PMAR,tMAR,alpha):
    ''' Return the derivative dc/dt at time, t, for given parameters.

        Parameters:
        -----------
        t: float

        c : float
            Concentration of nitrate across aquifer in mg/L**
        tc: float

        n : float
            Number of cows
        b : float
            Carbon sink parameter
        P : float
            Pressure across aquifer in MPa**
        Psurf : float
            Change in surface pressure **
        bc : float
            Recharge parameter **
        Pa : float
            Pressure drop across aquifer in MPa**
        PMAR: float

        tMAR: float
        
        alpha: float


        Returns:
        --------
        dcdt : float
                Rate of change of nitrate concentration

    '''
    if t>=tc:
        b=alpha*b
    if t>=tMAR:
        Pa+=PMAR

    return -n*b*(P-Psurf)+bc*c*(P-(Pa/2))

tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T

ni = np.interp(2000,tn,n) #interpolating number of cows
print(ni)