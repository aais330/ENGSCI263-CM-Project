import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from model_functions_posterior import *
from sklearn.linear_model import BayesianRidge

t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
t = np.arange(1999,2019.5,step = 0.05)

n_order = 1
cv, tv = np.array(c0), np.array(t0)
Tv = np.vander(tv, n_order+1, increasing=True)

br = BayesianRidge(fit_intercept=False, tol=1e-5)
br.fit(Tv, cv)

T0 = np.vander(t, n_order+1, increasing=True)
c, c_var = br.predict(T0, return_std=True)

sigma = [1.e-14]*len(c0)

p, cov = curve_fit(LPM_Model,t0,c0, sigma=sigma, p0=[0.2,0.5,0.5,1,1,5])
ps = np.random.multivariate_normal(p, cov, 100)
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