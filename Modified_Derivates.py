import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from model_functions import *

# Solve pressure ODE
pi = 0
t = np.arange(1998,2020,step = 1)
#parameters
b=0.5
b1=0.5
alpha=0
bc=1
tau = 5
# LMP_Model(t,b,b1,alpha,bc,tau)


# Testing curve_fit
# load in cow data and concentration data
tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T
tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

pars = curve_fit(LMP_Model,tcon,c,[1,1,1,1,15])
# print(pars)
b=pars[0][0]
b1=pars[0][1]
alpha=pars[0][2]
bc=pars[0][3]
tau = pars[0][4]
C = LMP_Model(t,b,b1,alpha,bc,tau)


f,ax = plt.subplots(1,1)

ax.plot(t,C,'k', label = 'Numeric Solution')
ax.plot(tcon,c,'r+', label = 'Data')
ax.set_title('Numerical Solution and data')
plt.ylabel('Concentration')
plt.xlabel('Time')

'''
# Plotting cows for reference
ax2 = ax.twinx()
ax2.plot(tn,n,'b', label = 'Data cows')
ax2.set_ylabel('Number of cows')
ax.legend() 
'''
    
# plt.show()
