import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from model_functions import *

t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
t = np.arange(1980,2020,step = 0.1)

t_pos = np.arange(1980,2019.75,step = 0.1)

ps, p = posterior_pars()

b1 = p[0]
alpha = p[1]
bc = p[2]
tau = p[3]

t_forecast = np.arange(2020,2030,step=0.05)

v=0.3
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.plot(t0,c0, 'ro', label="data", markersize=2.5)
# ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data', markersize=2.2)

#ax.plot(t, LPM_Model_forecast(t, b1, alpha, bc, tau, 0), 'k-', label='Forecast best-fit', alpha=0.5)

ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.0), 'm-')

ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.01), 'g-')

ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.05), 'b-')


ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.1), 'r-')

# for pi in ps:
#     ax.plot(t_pos, LPM_Model(t_pos, *pi), 'k-', alpha=0.5, lw=0.5)

#ax.plot(t, LPM_Model(t, *p), 'k-', label='Best fit', alpha=1)
ax.plot(t, LPM_Model_forecast(t, b1, alpha, bc, tau, 0), 'k-', label='Best-fit')
#ax.plot([], [], lw=0.5, label='posterior samples')
ax.plot([], [], 'm-', label='$dP_{mar}$ = 0.0 MPa')
ax.plot([], [], 'g-', label='$dP_{mar}$ = 0.01 MPa')
ax.plot([], [], 'b-', label='$dP_{mar}$ = 0.05 MPa')
ax.plot([], [], 'r-', label='$dP_{mar}$ = 0.1 MPa')
ax.legend(loc=2)
#plt.show()
plt.savefig("what_if_scenarios.png")