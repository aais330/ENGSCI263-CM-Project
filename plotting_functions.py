# This scipt contains the function that generate figures
import os
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from testing_functions import *

# Variables to be used by all functions
t0, c0 = np.genfromtxt("Data"+ os.sep +'nl_n.csv', delimiter=',', skip_header=1).T

def plot_data():
    '''
    Plots the given nitrate concentration and cow data

    Notes
    ------
    Saves the plot in file 'nc_data.png'
    '''

    # reading in cow data
    data_cow = np.genfromtxt("Data"+ os.sep +"nl_cows.txt", delimiter = ',', skip_header = 1)
    x_cow = data_cow[:,0]
    y_cow = data_cow[:,1]


    # plotting cow and concentration data
    fig1 = plt.figure(figsize = (15,7.5))
    ax = fig1.add_subplot(111)
   
    ax.plot(x_cow, y_cow, color = 'red', label = 'Cows')  # plotting cow data
    ax.set_xlabel('Time [yrs]', fontsize = 20)
    ax.set_ylabel('Number of Cows', fontsize = 20)

    ax2 = ax.twinx()
    ax2.plot(t0, c0, color = 'blue', label = 'Nitrate Concentration') # plotting nitrate data
    ax2.set_ylabel('Nitrate Concentration(mg/L)', fontsize = 20)

    plt.title('Annual Cattle Numbers in Southland and Nitrate Concentrations of Edendale Farm', fontsize = 20)
    ax.plot([],[], 'b', label='Nitrate Concentration')
    ax.legend()

    fig1.savefig('Plots'+ os.sep +'nc_data.png', dpi = 200)
    plt.close(fig1)



def initial_model():
    '''
    Plots the intial concentration model

    Notes
    ------
    Saves the plot in file 'initial_model.png'
    '''
    
    p = curve_fit(LPM_Model,t0,c0)[0] # solve model
    t = np.arange(1980,2020,step=0.05)
    C = LPM_Model(t,*p) 

    # plotting commands
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5) # plotting concentration data
    ax.plot(t, C, 'b-', label="Model") # plotting model
    plt.title("best fit LPM model")
    ax.set_xlabel('Time [yrs]')
    ax.set_ylabel('Nitrate Concentration [mg/L]')
    ax.legend()
    ax.legend(loc=2)

    fig.savefig('Plots'+ os.sep +'initial_model.png', dpi = 200)
    plt.close(fig)

def improved_model():
    '''
    Plots the improved concentration model

    Notes
    ------
    Saves the plot in file 'improved_model.png'
    '''
    p = posterior_pars()[1] # find improved parameters
    t = np.arange(1980,2020,step=0.05)
    C = LPM_Model(t,*p) 

    # plotting commands
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5) # plotting data
    ax.plot(t, C, 'b-', label="Model") # plotting model
    ax.set_xlabel('Time [yrs]')
    ax.set_ylabel('Nitrate Concentration [mg/L]')
    plt.title("best fit LPM model")
    ax.legend()
    ax.legend(loc=2)

  
    fig.savefig('Plots'+ os.sep +'improved_model.png', dpi = 200)
    plt.close(fig)


def what_ifs():
    '''
    Plots the several different 'what-if' scenarios for different managed aquifer recharge pressures.

    Notes
    ------
    Saves the plot in file 'what_if_scenarios.png'
    '''
    t_mas = np.arange(1980,2030,step = 0.1)
    MAS = [11.3]*len(t_mas)
    t = np.arange(1980,2020,step = 0.1) #time until 2020

    p = posterior_pars()[1] #fixed parameters with no uncertainty

    b1 = p[0]
    alpha = p[1]
    bc = p[2]
    tau = p[3]

    t_forecast = np.arange(2020,2030,step=0.05) #10 year projection time
    #LPM_Model_forecast(t, b1, alpha, bc, tau, 0.1)

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5)
    # ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data', markersize=2.2)

    #ax.plot(t, LPM_Model_forecast(t, b1, alpha, bc, tau, 0), 'k-', label='Forecast best-fit', alpha=0.5)

    ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.0), 'm-', label='$dP_{MAR}$ = 0.0 MPa')
    
    ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.05), 'g-', label='$dP_{MAR}$ = 0.05 MPa')
    
    ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.1), 'b-',label='$dP_{MAR}$ = 0.10 MPa' )
    
    ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.15), 'r-', label='$dP_{MAR}$ = 0.15 MPa')
    
    ax.plot(t, LPM_Model_forecast(t, b1, alpha, bc, tau, 0), 'k-', label='Best-fit model')
    ax.plot(t_mas, MAS, 'k--', alpha = 0.7, label='Maximum acceptable standard')
    plt.title("Edendale Aquifer LPM: different MAR pressures")
    ax.set_xlabel('Time [yrs]')
    ax.set_ylabel('Nitrate Concentration [mg/L]')
    ax.legend(loc=2)
    plt.savefig('Plots'+ os.sep +"what_if_scenarios.png")
    plt.close(fig)

def without_acs():
    '''
    Plots the comparision between nitrate concentrations with and without an active carbon sink.

    Notes
    ------
    Saves the plot in file 'without_acs.png'
    '''
    
    t = np.arange(1980,2020,step = 0.1)
    MAS = [11.3]*len(t)

    p = posterior_pars()[1]

    b1 = p[0]
    bc = p[2]
    tau = p[3]

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5)

    ax.plot(t, LPM_Model(t, b1, 1, bc, tau), 'g-', label='Forecast without ACS') #set alpha = 1 to nullify effect of ACS
    ax.plot(t, LPM_Model(t, *p), 'k-', label='Forecast with ACS') #model with ACS
    ax.plot(t, MAS, 'k--', alpha = 0.7, label='Maximum acceptable standard')
    ax.legend(loc=2)
    
    #printing 2020 concentration values to find difference in concentration with and without ACS
    #print(LPM_Model(t, b1, 1, bc, tau)[-1])
    #print(LPM_Model(t, *p)[-1])
    

    plt.title("Edendale Aquifer LPM: with and without ACS")
    ax.set_xlabel('Time [yrs]')
    ax.set_ylabel('Nitrate Concentration [mg/L]')
    
    fig.savefig('Plots'+ os.sep +'without_acs.png', dpi = 200)
    plt.close(fig)



def what_ifs_uncertainty():
    '''
    Plots the several different 'what-if' scenarios for different managed aquifer recharge pressures but
    taking into account uncertainty.

    Notes
    ------
    Saves the plot in file 'what_if_scenarios.png'
    '''
    t = np.arange(1980,2030,step = 0.1)
    MAS = [11.3]*len(t)
    t_pos = np.arange(1980,2020.01,step = 0.05)

    ps = posterior_pars()[0] #100 different parameter sets due to uncertainty


    t_forecast = np.arange(2019.99,2030,step=0.05)

    v=0.15
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    
    C = np.zeros(ps.shape[0])

    #looping through all parameters in parameter set with uncertainty, finding all 2030 concentration values,
    #then finding a confidence interval for the 2030 concentration values
    for pi in range(0,ps.shape[0]):
        C[pi] = LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0)[-1]
        ax.plot(t_forecast, LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0), 'm-', lw=0.3)
    confidence_int(C,"0.00 MPa")

    for pi in range(0,ps.shape[0]):
        C[pi] = LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.05)[-1]
        ax.plot(t_forecast, LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.05), 'g-', lw=0.3)
    confidence_int(C,"0.05 MPa")

    for pi in range(0,ps.shape[0]):
        C[pi] = LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.1)[-1]
        ax.plot(t_forecast, LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.1), 'b-', lw=0.3)
    confidence_int(C,"0.10 MPa")

    for pi in range(0,ps.shape[0]):
        C[pi] =  LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.15)[-1]
        ax.plot(t_forecast, LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.15), 'r-', lw=0.3)
    confidence_int(C,"0.15 MPa")

    for pi in range(0,ps.shape[0]):
        ax.plot(t_pos, LPM_Model_forecast(t_pos, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0), 'k-', lw=0.3)

    ax.errorbar(t0,c0,yerr=v,fmt='ro', label='Data', markersize=2.5)
    ax.plot(t, MAS, 'k--', alpha = 0.7, label='Maximum acceptable standard')

    ax.plot([], [], 'k-', label='posterior samples')
    ax.plot([], [], 'm-', label='$dP_{MAR}$ = 0.0 MPa')
    ax.plot([], [], 'g-', label='$dP_{MAR}$ = 0.05 MPa')
    ax.plot([], [], 'b-', label='$dP_{MAR}$ = 0.10 MPa')
    ax.plot([], [], 'r-', label='$dP_{MAR}$ = 0.15 MPa')
    ax.legend(loc=2)
    plt.title("Edendale Aquifer LPM: with and without ACS")
    ax.set_xlabel('Time [yrs]')
    ax.set_ylabel('Nitrate Concentration [mg/L]')
    
    plt.savefig('Plots'+ os.sep +"what_if_uncertainty.png")
    plt.close(fig)

def without_acs_uncertainty():
    '''
    Plots the comparision between nitrate concentrations with and without
    an active carbon sink but taking into account uncertainty

    Notes
    ------
    Saves the plot in file 'without_acs_uncertainty.png'
    '''
    t = np.arange(1980,2030,step = 0.1)
    MAS = [11.3]*len(t)
    t_pos = np.arange(1980,2020,step = 0.1)

    ps = posterior_pars()[0]

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    
    v = 0.15
    ax.errorbar(t0,c0,yerr=v,fmt='ro', label='Data', markersize=2.5)
    C = np.zeros(ps.shape[0])

    # Without ACS
    for pi in range(0,ps.shape[0]):
        C[pi] = LPM_Model_forecast(t_pos, ps[pi][0], 1, ps[pi][2], ps[pi][3],0)[-1]
        ax.plot(t_pos, LPM_Model_forecast(t_pos, ps[pi][0], 1, ps[pi][2], ps[pi][3],0), 'g-', lw=0.3)

    # With ACS
    for pi in range(0,ps.shape[0]):
        C[pi] -= LPM_Model_forecast(t_pos, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0)[-1]
        ax.plot(t_pos, LPM_Model_forecast(t_pos, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0), 'k-', lw=0.3)

    confidence_int(C,"2030 concentration values with vs without ACS")
   
    ax.plot([],[], 'k-', label = 'With ACS')
    ax.plot([],[], 'g-', label = 'Without ACS')
    ax.plot(t, MAS, 'k--', alpha = 0.7, label='Maximum acceptable standard')
    ax.legend()
    plt.title("Edendale Aquifer LPM: With and Without ACS")
    ax.set_xlabel('Time [yrs]')
    ax.set_ylabel('Nitrate Concentration [mg/L]')
    
    fig.savefig('Plots'+ os.sep +'without_acs_uncertainty.png', dpi = 200)
    plt.close(fig)

def misfit_plot(old):
    '''
    Produces a misfit plot

    Parameters
    ----------
    old: bool
        If intial is true will plot misfit of old model, otherwise will plot misfit
        of improved model
    '''

    t = np.arange(1980,2030,step = 0.1)
    t0, c0 = np.genfromtxt("Data"+ os.sep +'nl_n.csv', delimiter=',', skip_header=1).T

    if old:
        parameters = curve_fit(LPM_Model,t0,c0)[0]
    else:
        parameters = posterior_pars()[1]

    C = LPM_Model(t,*parameters)


    misfit=[] # stores misfit values

    # Interpolate concentration at data points for model
    C=np.interp(t0,t,C)

    # Calculate misfits
    for i in range(len(C)):
        misfit.append(c0[i]-C[i])

    # Plot
    f, ax2 = plt.subplots(1,1)
    ax2.plot(t0,misfit,'bx', label = 'Misfit')
    ax2.axhline(0., c='k', ls=':')
    ax2.set_title('Misfit of Best fit LPM model')
    plt.ylabel('Concentration misfit (mg/L)')
    plt.xlabel('Time (Years)')

    if old:
        file = 'misfit_plot_old_model'
    else:
        file = 'misfit_plot_new_model'

    #plt.show()
    f.savefig('Plots'+ os.sep + file, dpi = 200)
    plt.close(f)

def alpha_distribution():
    """
    Produces a histogram, with bounds for the 90 percent confidence interval for parameter alpha

    Notes: 
    ----------
    No parameters take, calls other functions. Prints confidence interval of alpha to screen
    """
    alpha_list = []
    #looping through 10 iterations of posterior_pars() uncertainty values
    for i in range(10):
        ps = posterior_pars()[0]
        
        for j in range(ps.shape[0]):
            alpha_list.append(ps[j][1])
        
    
    ci = confidence_int(alpha_list, "alpha") #confidence interval for alpha
    f, ax = plt.subplots(1,1)

    #plotting histogram and 90 percent conf int for alpha
    ax.hist(alpha_list, bins = 30, density = True)
    ax.axvline(x=ci[0], color = 'r')
    ax.axvline(x=ci[1], color = 'r', label = '90 percent confidence interval')
    plt.ylabel('Probability Density')
    plt.xlabel('Alpha')
    ax.legend()
    ax.set_title('Posterior Distribution of Alpha Parameter')

    #plt.show()
    f.savefig('Plots'+ os.sep + 'alpha_distribution.png', dpi = 300)
    plt.close(f)