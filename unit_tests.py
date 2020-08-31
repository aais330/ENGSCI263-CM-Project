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
tol=1.e-10

def func(x,y,pars=[]):
    pass

def test_ie():
    #simple function testing
    pass

def test_dPdt():
    # t < t_mar
    check_dPdt = dPdt(10, 2015, 3)
    assert((check_dPdt+60)<tol)
    
    check_dPdt2 = dPdt(15, 2019, 0.2)
    assert((check_dPdt2+6)<tol)
    
    # t > t_mar
    check_dPdt3 = dPdt(20, 2025, 6)
    assert((check_dPdt3)+238.5<tol)

    # test when switches to tmar on
    # t < t_mar
    check_dPdt4 = improved_euler_step(dPdt, 2019, -6, 1, pars=[0.2])
    assert((check_dPdt4+4.08)<tol) #calculation gave me 11.975, maybe i missed something?

    check_dPdt5 = improved_euler_step(dPdt, 1998, 6, 0.2, pars=[0.2])
    assert((check_dPdt5-5.54)<tol) # Keep getting 5.54 even though the year changes

    # t > t_mar
    check_dPdt6 = improved_euler_step(dPdt, 2021, 6, 0.2, pars=[0.2])
    assert((check_dPdt6-5.55)<tol)

    print("dPdt passed \n")

def test_dCdt():
    # t-tau >= 1990.5
    check_dCdt = dCdt(20, 2005, 0.4, 1, 2, 1, 3) 
    assert((check_dCdt+117480.18)<tol)
    
    # t-tau >= 1990.5 and t-tau>t_acs
    check_dCdt2 = dCdt(25, 2016, 0.2, 1, 2, 1, 5) 
    assert((check_dCdt2-212666)<tol)

    # t-tau < 1990.5
    check_dCdt3 = dCdt(25, 1992, 0.2, 1, 2, 1, 3) 
    assert((check_dCdt3+2996.25)<tol)

    # t>t_mar and t-tau >= 1990.5 and t-tau>t_acs
    check_dCdt4 = dCdt(25, 2021, 0.2, 1, 2, 1, 3)
    assert((check_dCdt4+200948.05)<tol)

    # Tests when t-tau > 1990.5 
    # Test when t-tau>t_acs 2015>2010
    check_dCdt5=improved_euler_step(dCdt,2020,0.2,1,pars=[0,b1,alpha,bc,tau])
    assert((check_dCdt5-0.1665)<tol)

    # Test when t-tau<t_acs 2009<2010
    check_dCdt6=improved_euler_step(dCdt,2014,0.2,1,pars=[0,b1,alpha,bc,tau])
    assert((check_dCdt6-13870.7)<tol)

    # Test when t>t_mar
    check_dCdt7=improved_euler_step(dCdt,2021,0.2,1,pars=[0,b1,alpha,bc,tau]) 
    assert((check_dCdt7-0.149)<tol) 

    # Test when t<t_mar
    check_dCdt8=improved_euler_step(dCdt,2018,0.2,1,pars=[0,b1,alpha,bc,tau])
    assert((check_dCdt8-0.1903)<tol)

    print("dCdt passed \n")

# Do we need to add tests for LMP_Model and solve_dCdt?

test_dPdt()
test_dCdt()