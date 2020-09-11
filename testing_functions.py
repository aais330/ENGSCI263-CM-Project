import numpy as np
from model_functions import *
tol=1.e-8

def func(x,y,pars=[]):
    pass

def test_ie():
    #simple function testing
    pass

def test_dPdt():
    # t < t_mar before t_mar
    check_dPdt = dPdt(10, 2015)
    assert((check_dPdt+20)<tol)
    
    check_dPdt2 = dPdt(15, 2019)
    assert((check_dPdt2+30)<tol)
    
    # t > t_mar after t_mar
    check_dPdt3 = dPdt_forecast(20, 2025,.5)
    assert((check_dPdt3)+39.75<tol)

    check_dPdt4 = dPdt_forecast(20, 2020.1,.5)
    assert((check_dPdt4)+39.75<tol)

    # test when switches to tmar on
    # THINK THIS IS TESTING IMPROVED EULER NOT dPdt
    # # t < t_mar
    # check_dPdt4 = improved_euler_step(dPdt, 2019, -6, 1, pars=[0.2])
    # assert((check_dPdt4+4.08)<tol) #calculation gave me 11.975, maybe i missed something?

    # check_dPdt5 = improved_euler_step(dPdt, 1998, 6, 0.2, pars=[0.2])
    # assert((check_dPdt5-5.54)<tol) # Keep getting 5.54 even though the year changes

    # # t > t_mar
    # check_dPdt6 = improved_euler_step(dPdt, 2021, 6, 0.2, pars=[0.2])
    # assert((check_dPdt6-5.55)<tol)


    print("dPdt passed \n")


# dCdt(ci, t, P, b1, alpha, bc, tau):

def test_dCdt():

    # between 1990.5 and before active carbon sink and no time lag
    check_dCdt = dCdt(50, 1999.5, 0, 1, 1, 1, 0) 
    assert((check_dCdt-11645.8)<tol)

    # between 1990.5 and before active carbon sink and time lag
    check_dCdt = dCdt(50, 2000.5, 0, 1, 1, 1, 1) 
    assert((check_dCdt-11645.8)<tol)
    
    # After active carbon sink no time lag
    check_dCdt = dCdt(50, 2010.5, 0, 1, 0.5, 1, 0) 
    assert((check_dCdt-14977.45)<tol)

    # After active carbon sink with time lag
    check_dCdt = dCdt(50, 2011.5, 0, 1, 0.5, 1, 1) 
    assert((check_dCdt-14977.45)<tol)

    # before 1990.5
    check_dCdt = dCdt(50,1985,0, 1, 1, 1, 0) 
    assert((check_dCdt-1886.1)<tol)

    # after 2020 without mar, with active carbon sink
    check_dCdt = dCdt(50,2021,0, 1, 0.5, 1, 0) 
    assert((check_dCdt-15903.525)<tol)

    # after 2020 without mar, without active carbon sink
    check_dCdt = dCdt(50,2021,0, 1, 1, 1, 0) 
    assert((check_dCdt-31809.55)<tol)


    # after 2020 with mar, with active carbon sink
    check_dCdt = dCdt_forecast(50,2021,0, 1, 0.5, 1, 0,0.5) 
    assert((check_dCdt-15916.025)<tol)

    # Testing Improved euler not dCdt
    # # Tests when t-tau > 1990.5 
    # # Test when t-tau>t_acs 2015>2010
    # check_dCdt5=improved_euler_step(dCdt,2020,0.2,1,pars=[0,b1,alpha,bc,tau])
    # assert((check_dCdt5-0.1665)<tol)

    # # Test when t-tau<t_acs 2009<2010
    # check_dCdt6=improved_euler_step(dCdt,2014,0.2,1,pars=[0,b1,alpha,bc,tau])
    # assert((check_dCdt6-13870.7)<tol)

    # # Test when t>t_mar
    # check_dCdt7=improved_euler_step(dCdt,2021,0.2,1,pars=[0,b1,alpha,bc,tau]) 
    # assert((check_dCdt7-0.149)<tol) 

    # # Test when t<t_mar
    # check_dCdt8=improved_euler_step(dCdt,2018,0.2,1,pars=[0,b1,alpha,bc,tau])
    # assert((check_dCdt8-0.1903)<tol)

    print("dCdt passed \n")

# convergence testing solutions for different step sizes
def convergence_sols(step_sizes, pars=[]):
    # pre-allocating solutions array
    sols = []

    for step in 1/step_sizes:
        # time range with step size
        t = np.arange(1998, 2020, step)

        C = LPM_Model(t, pars[0], pars[1], pars[2], pars[3])
        C = np.interp(2012, t, C)
        # add concentration solution at each step size to array
        sols.append(C)
    
    return sols 
