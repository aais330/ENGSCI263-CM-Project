import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from model_functions_posterior import *
from sklearn.linear_model import BayesianRidge

t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
t = np.arange(1980,2030,step = 0.05)
'''
n_order = 1
cv, tv = np.array(c0), np.array(t0)
Tv = np.vander(tv, n_order+1, increasing=True)

br = BayesianRidge(fit_intercept=False, tol=1e-5)
br.fit(Tv, cv)

T0 = np.vander(t, n_order+1, increasing=True)
c, c_var = br.predict(T0, return_std=True)
'''
sigma = [1.e-14]*len(c0) # variance limit of pars
ps, p = posterior_pars(sigma)
v=0.3
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data', markersize=2.2)
ax.plot(t, LPM_Model(t, *p), 'b-', label='best-fit')

for pi in ps:
    ax.plot(t, LPM_Model(t, *pi), 'k-', alpha=0.3, lw=0.5)

#ax.fill_between(t, c-2*c_var, c+2*c_var, color='r', alpha=0.3, label='2$\sigma$ envelope')
ax.plot([], [], lw=0.5, label='posterior samples')
ax.set_ylim(0,16)
ax.legend()
plt.show()
#plt.savefig("posterior.png")


'''
C = LPM_Model(t, *p)
tcon, c = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

misfit=[]
# Only include data for both cows and concentration
i=0 # start commenting here for original model
while tcon[i]<1999:
    i+=1
tcon=tcon[i:]
c=c[i:] # for original model stop commenting here

# Interpolate concentration at data points for model
C=np.interp(tcon,t,C)

# Append differences to misfit array
for i in range(len(C)):
    misfit.append(c[i]-C[i])

# Plot
f, ax2 = plt.subplots(1,1)
ax2.plot(tcon,misfit,'bx', label = 'Misfit')
ax2.axhline(0., c='k', ls=':')
ax2.set_title('Best fit LPM model')
plt.ylabel('Concentration misfit (mg/L)')
plt.xlabel('Time (Years)')
plt.show()
'''