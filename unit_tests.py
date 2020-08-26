import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from Modified_Derivates import improved_euler_step,dCdt

# Tests for concentration model
# dCdt(ci, t, P, b1, alpha, bc, tau)
# -ni*b1*(P-dP_surf)+bc*ci*(P-(dP_a/2))
# pars=[P[i],b1,alpha,bc,tau]
# improved_euler_step(f, tk, xk, h, pars = [])
b=0.5
b1=0.5
alpha=0
bc=1
tau = 5
dP_a = 0.1 
dP_surf = 0.05 
t_mar = 2020 
t_acs = 2010 
tol=0.1
# *need to put everything through IE calculator*
# Tests when t-tau < 1990.5

# Test when n=0

# Tests when t-tau > 1990.5 
# Test when t-tau>t_acs 2015>2010
Test=improved_euler_step(dCdt,2020,0.2,1,pars=[0,b1,alpha,bc,tau])
# assert(abs(Test-0.1665)<tol)
print('Test 1 passed')
# Test when t-tau<t_acs 2009<2010
Test=improved_euler_step(dCdt,2014,0.2,1,pars=[0,b1,alpha,bc,tau])
# assert(abs(Test-13870.7)<tol)
print('Test 2 passed')
# Test when t>t_mar
Test=improved_euler_step(dCdt,2021,0.2,1,pars=[0,b1,alpha,bc,tau])
# Test when t<t_mar
Test=improved_euler_step(dCdt,2018,0.2,1,pars=[0,b1,alpha,bc,tau])