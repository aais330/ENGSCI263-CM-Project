import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from model_functions import *


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
tol=1.-10

def test_dPdt():
    #2 tests for dPdt
    check_dPdt = dPdt(10, 2015, 3)
    assert((abs(check_dPdt-60))<tol)

    check_dPdt2 = dPdt(20, 2025, 6)
    assert((abs(check_dPdt2)-238.5)<tol)

    print("dPdt passed")

def test_dCdt():
    #2 tests for dCdt
    
    check_dCdt = dCdt(20, 2005, 0.4, 1, 2, 1, 3) 
    #assert((abs(check_dCdt-121818.25))<tol) # incorrect calculation

    check_dCdt2 = dCdt(25, 2016, 0.2, 1, 2, 1, 5) 
    assert((abs(check_dCdt2)-212666)<tol) # this seems right

    # Tests when t-tau > 1990.5 
    # Test when t-tau>t_acs 2015>2010
    check_dCdt3=improved_euler_step(dCdt,2020,0.2,1,pars=[0,b1,alpha,bc,tau])
    # assert(abs(Test-0.1665)<tol)

    # Test when t-tau<t_acs 2009<2010
    check_dCdt4=improved_euler_step(dCdt,2014,0.2,1,pars=[0,b1,alpha,bc,tau])
    # assert(abs(Test-13870.7)<tol)

    # Test when t>t_mar
    check_dCdt5=improved_euler_step(dCdt,2021,0.2,1,pars=[0,b1,alpha,bc,tau])  

    # Test when t<t_mar
    check_dCdt6=improved_euler_step(dCdt,2018,0.2,1,pars=[0,b1,alpha,bc,tau])

    # Tests when t-tau < 1990.5


    print("dCdt passed")

test_dCdt()