import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from model_functions_posterior import *
from model_ensemble_project import *
from sklearn.linear_model import BayesianRidge
import cProfile, pstats

t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
t = np.arange(1980,2030,step = 0.1)

t_pos = np.arange(1980,2020,step = 0.1)

ps, p = posterior_pars()
print(p)

b1 = p[0]
alpha = p[1]
bc = p[2]
tau = p[3]

t_forecast = np.arange(2020,2030,step=0.05)

v=0.3
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
#ax.plot(t0,c0, 'ro', label="data", markersize=2.5)
ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data', markersize=2.5)

#ax.plot(t, LPM_Model_forecast(t, b1, alpha, bc, tau, 0), 'k-', label='Forecast best-fit', alpha=0.5)

#posterior plotting
# for i in np.linspace(0,0.005,10):
#     ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, i), 'm-', lw=0.5)

# for i in np.linspace(0.007,0.013,10):
#     ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, i), 'g-', lw=0.5)

# for i in np.linspace(0.04,0.06,10):
#     ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, i), 'b-', lw=0.5)

# for i in np.linspace(0.09,0.11,10):
#     ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, i), 'r-', lw=0.5)

# for pi in range(0,ps.shape[0]):
#     ax.plot(t, LPM_Model_forecast(t, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.05), 'g-', lw=0.5)

# for pi in range(0,ps.shape[0]):
#     ax.plot(t, LPM_Model_forecast(t, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.1), 'b-', lw=0.5)

# for pi in range(0,ps.shape[0]):
#     ax.plot(t, LPM_Model_forecast(t, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.15), 'r-', lw=0.5)

# for pi in range(0,ps.shape[0]):
#     ax.plot(t, LPM_Model_forecast(t, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0), 'k-', lw=0.5)

#effects of acs
ax.plot(t, LPM_Model(t, b1,1,bc,tau), 'b-', label='Best fit', alpha=1)
ax.plot(t, LPM_Model(t, b1, alpha, bc, tau), 'k-', label='Forecast Best-fit')
ax.plot([], [], 'k-', lw=0.5, label='posterior samples')
ax.plot([], [], 'm-', label='$dP_{mar}$ = 0.0 MPa')
ax.plot([], [], 'g-', label='$dP_{mar}$ = 0.01 MPa')
ax.plot([], [], 'b-', label='$dP_{mar}$ = 0.05 MPa')
ax.plot([], [], 'r-', label='$dP_{mar}$ = 0.1 MPa')
ax.legend(loc=2)
plt.show()
#plt.savefig("posterior_ensemble.png")