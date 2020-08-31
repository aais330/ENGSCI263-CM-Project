import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from model_functions_posterior import *

t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
print(t0)
t = np.arange(1980,2020,step = 0.5)

sigma = [0.5]*len(c0)

p, cov = curve_fit(LPM_Model,t0,c0, sigma=sigma, p0=[0.2,0.5,0.5,1,1,5])
ps = np.random.multivariate_normal(p, cov, 75)
v=0.5
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data', markersize=1.8)
ax.plot(t, LPM_Model(t, *p), 'b-', label='best-fit')
for pi in ps:
    ax.plot(t, LPM_Model(t, *pi), 'k-', alpha=0.3, lw=0.5)
ax.plot([], [], lw=0.5, label='posterior samples')
ax.set_ylim(0,16)
ax.legend()
plt.show()