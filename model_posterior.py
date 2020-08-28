import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from Modified_Derivates import LMP_Model, solve_dCdt, solve_dPdt, dCdt, dPdt, improved_euler_step

#reading data
t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
t = np.arange(1980,2020,step = 0.5)
'''
#using curve_fit() to find best parameters to input
pars = curve_fit(LMP_Model,t0,c0,[1,1,1,1,15])
b = pars[0][0]
b1 = pars[0][1]
alpha = pars[0][2]
bc = pars[0][3]
tau = pars[0][4]
C = LMP_Model(t,b,b1,alpha,bc,tau)
'''

#setting variance
v = 0.5
#calculate the sum-of-squares objective function 
#S = np.sum(np.square(np.subtract(c0, C)/v))

# plotting commands
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.plot(t0,c0,'b-', label='model')
ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data')
ax.set_xlabel('time')
ax.set_ylabel('pressure')
#ax.set_title('objective function: S={:3.2f}'.format(S))
ax.legend()
plt.show()
#plt.savefig('task1.png')