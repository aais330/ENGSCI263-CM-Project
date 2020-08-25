
import numpy as np
from matplotlib import pyplot as plt
import math
from concmodel import *
import pandas as pd

def IEStep(f, t, c, step, tc,n,b,P,Psurf,bc,Pa,PMAR,tMAR,alpha,n1,b1,P1):

    # Predictor step using the standard euler method
    # Finding x*[n+1]
    cStark1 = c + step*f(t,c,tc,n,b,P,Psurf,bc,Pa,PMAR,tMAR,alpha)
    # Original function with t[n] and x[n]
    fn = f(t, c, tc,n,b,P,Psurf,bc,Pa,PMAR,tMAR,alpha)
    # Corrector step
    # Finding f*[n+1] with t[n+1] and x*[n+1]
    fStark1 = f(t+step,cStark1, tc,n1,b1,P1,Psurf,bc,Pa,PMAR,tMAR,alpha)
    ck1 = c + 0.5*step*(fn+fStark1)

    return ck1

# load in cow data and concentration data
tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T
tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

# # initialise constants
step=0.5
tc=2020
b=0.5
P=0
Psurf=0.05
bc=1
Pa=0.1
PMAR=0
tMAR=2020
alpha=1

nx = int(np.ceil((tcon[-1]-tcon[0])/step))
# compute number of Improved Euler steps to take
ts = tcon[0]+np.arange(nx+1)*step
# t array
cs = 0.*ts
# array to store solution
cs[0] = 0.2
# set initial value

# interpolate cow values and conc values
n_int=np.interp(ts,tn,n)
c_int=np.interp(ts,tcon,c)

for i in range(len(ts)-1):
    P1=0
    b1=1
    n1=n_int[i+1]
    cs[i+1]=IEStep(conc_model,ts[i],c_int[i],step,tc,n_int[i],b,P,Psurf,bc,Pa,PMAR,tMAR,alpha,n1,b1,P1)

#take away first 10 (as no data for those years)
# tn=tcon[10:]
# n_int=n_int[10:]

# plot
# plt.plot(tn,n,'or',label='cow data')
# plt.plot(tcon,c,'ok',label='conc data')
plt.plot(ts,cs,'b',label='IE solution')
# plt.legend()
plt.show()
        