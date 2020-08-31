import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from model_functions import *

# Solve pressure ODE
pi = 0
# t = np.arange(1999,2019.25,step =0.25)
t = np.arange(1980,2020,step =0.25)

# Testing curve_fit
# load in cow data and concentration data
tn, n = np.genfromtxt('nl_cows.txt', delimiter=',', skip_header=1).T
tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

# pars, cov = curve_fit(LMP_Model,tcon,c,[1,1,1,1,15], sigma= sigma)
pars= curve_fit(LMP_Model,tcon,c,[1,1,1,1,15])
#  print(pars)
b=pars[0]
b1=pars[1]
alpha=pars[2]
bc=pars[3]
tau = pars[4]
C = LMP_Model(t,b,b1,alpha,bc,tau)

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)

ax.plot(t,C,'k', label = 'Numeric Solution')
ax.plot(tcon,c,'ro', markersize=2.0, label = 'Data')
ax.set_title('Numerical Solution vs Data')
plt.ylabel('Concentration(mg/L)')
plt.xlabel('Time(Years)')
plt.legend()
plt.show()

# Producing misfit plot
# Pre-allocating array
misfit=[]
print(len(t))
print(len(tcon))
# i=0
# Only include data for both cows and concentration
# while tcon[i]<1999:
#     i+=1
# tcon=tcon[i:]
# c=c[i:]
# Interpolate concentration at data points for model
# C=np.interp(tcon,t,C)
# Append differences to misfit array
# for i in range(len(C)):
#     misfit.append(c[i]-C[i])

# Plot
# f, ax2 = plt.subplots(1,1)
# ax2.plot(tcon,misfit,'bx', label = 'Misfit')
# ax2.axhline(0., c='k', ls=':')
# ax2.set_title('Best fit LPM model')
# plt.ylabel('Concentration misfit (mg/L)')
# plt.xlabel('Time (Years)')
# plt.show()

'''
# Plotting cows for reference
ax2 = ax.twinx()
ax2.plot(tn,n,'b', label = 'Data cows')
ax2.set_ylabel('Number of cows')
ax.legend() 
'''

    
# plt.show()
#plt.savefig("model.png")
