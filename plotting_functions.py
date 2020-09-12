# This scipt contains the function that generate figures
import numpy as np
from matplotlib import pyplot as plt
from testing_functions import *

def plot_data():
    '''
    Plots the given nitrate concentration and cow data

    Notes
    ------
    Saves the plot in file 'nc_data.png'
    '''

    # reading in cow data
    data_cow = np.genfromtxt("nl_cows.txt", delimiter = ',', skip_header = 1)
    x_cow = data_cow[:,0]
    y_cow = data_cow[:,1]

    # reading in concentration data
    data_nconc = np.genfromtxt("nl_n.csv", delimiter = ',', skip_header = 1)
    x_nconc = data_nconc[:,0]
    y_nconc = data_nconc[:,1]

    # plotting cow and concentration data
    fig1 = plt.figure(figsize = (15,7.5))
    ax = fig1.add_subplot(111)
    ax.plot(x_cow, y_cow, color = 'red', label = 'Cows')
    ax.set_xlabel('Years', fontsize = 20)
    ax.set_ylabel('Number of Cows', fontsize = 20)
    ax2 = ax.twinx()
    ax2.plot(x_nconc, y_nconc, color = 'blue', label = 'Nitrate Concentration')
    ax2.set_ylabel('Nitrate Concentration(mg/L)', fontsize = 20)
    plt.title('Annual Cattle Numbers in Southland and Nitrate Concentrations of Edendale Farm', fontsize = 20)
    ax.plot([],[], 'b', label='Nitrate Concentration')
    ax.legend()
    #plt.show()
    fig1.savefig('nc_data.png', dpi = 200)


def initial_model():
    '''
    Plots the intial concentration model

    Notes
    ------
    Saves the plot in file 'initial_model.png'
    '''
    p = posterior_pars_old()
    t = np.arange(1980,2020,step=0.05)
    C = LPM_Model(t,*p)


    t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5)
    ax.plot(t, C, 'b-', label="Model")
    plt.title("Initial Model", fontsize=20)
    ax.legend()
    #plt.show()
    fig.savefig('initial_model.png', dpi = 200)

def improved_model():
    '''
    Plots the improved concentration model

    Notes
    ------
    Saves the plot in file 'improved_model.png'
    '''
    pos, p = posterior_pars()
    t = np.arange(1980,2030,step=0.05)
    C = LPM_Model(t,*p)


    t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5)
    ax.plot(t, C, 'b-', label="Model")
    plt.title("Improved Model", fontsize=20)
    ax.legend()
    #plt.show()
    fig.savefig('improved_model.png', dpi = 200)


def what_ifs():
    '''
    Plots the several different 'what-if' scenarios for different managed aquifer recharge pressures.

    Notes
    ------
    Saves the plot in file 'what_if_scenarios.png'
    '''
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
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5)
    # ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data', markersize=2.2)

    #ax.plot(t, LPM_Model_forecast(t, b1, alpha, bc, tau, 0), 'k-', label='Forecast best-fit', alpha=0.5)

    ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.0), 'm-')

    ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.05), 'g-')

    ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.1), 'b-')

    ax.plot(t_forecast, LPM_Model_forecast(t_forecast, b1, alpha, bc, tau, 0.15), 'r-')

    ax.plot(t, LPM_Model_forecast(t, b1, alpha, bc, tau, 0), 'k-', label='Best-fit model')
    ax.plot([], [], 'm-', label='$dP_{mar}$ = 0.0 MPa')
    ax.plot([], [], 'g-', label='$dP_{mar}$ = 0.05 MPa')
    ax.plot([], [], 'b-', label='$dP_{mar}$ = 0.10 MPa')
    ax.plot([], [], 'r-', label='$dP_{mar}$ = 0.15 MPa')
    plt.title("Future scenarios for potential values of $dP_{mar}", fontsize=20)
    ax.legend(loc=2)
    #plt.show()
    plt.savefig("what_if_scenarios.png")

def without_acs():
    '''
    Plots the comparision between nitrate concentrations with and without an active carbon sink.

    Notes
    ------
    Saves the plot in file 'without_acs.png'
    '''
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
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5)

    ax.plot(t, LPM_Model(t, *p), 'k-', label='Forecast with ACS', alpha=0.5)
    ax.plot(t, LPM_Model(t, b1, 1, bc, tau), 'g-', label='Forecast without ACS ', alpha=0.5)
    ax.legend(loc=2)
    plt.title("Improved model compared to model without ACS parameter", fontsize=20)
    #plt.show()
    fig.savefig('without_acs.png', dpi = 200)


def what_ifs_uncertainty():
    '''
    Plots the several different 'what-if' scenarios for different managed aquifer recharge pressures but
    taking into account uncertainty.

    Notes
    ------
    Saves the plot in file 'what_if_scenarios.png'
    '''
    
    t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
    t = np.arange(1980,2030,step = 0.1)

    t_pos = np.arange(1980,2020,step = 0.1)

    ps, p = posterior_pars()


    t_forecast = np.arange(2020,2030,step=0.05)

    v=0.3
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    #ax.plot(t0,c0, 'ro', label="data", markersize=2.5)


    for pi in range(0,ps.shape[0]):
        ax.plot(t_forecast, LPM_Model_forecast(t_forecast, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0), 'm-', lw=0.3)

    for pi in range(0,ps.shape[0]):
        ax.plot(t, LPM_Model_forecast(t, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.05), 'g-', lw=0.3)

    for pi in range(0,ps.shape[0]):
        ax.plot(t, LPM_Model_forecast(t, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.1), 'b-', lw=0.3)

    for pi in range(0,ps.shape[0]):
        ax.plot(t, LPM_Model_forecast(t, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0.15), 'r-', lw=0.3)

    for pi in range(0,ps.shape[0]):
        ax.plot(t_pos, LPM_Model_forecast(t_pos, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0), 'k-', lw=0.3)

    ax.errorbar(t0,c0,yerr=v,fmt='ro', label='Data', markersize=2.5)

    ax.plot([], [], 'k-', label='posterior samples')
    ax.plot([], [], 'm-', label='$dP_{mar}$ = 0.0 MPa')
    ax.plot([], [], 'g-', label='$dP_{mar}$ = 0.05 MPa')
    ax.plot([], [], 'b-', label='$dP_{mar}$ = 0.10 MPa')
    ax.plot([], [], 'r-', label='$dP_{mar}$ = 0.15 MPa')
    ax.legend(loc=2)
    plt.title("Potential effects caused by different values of dP_{mar}$", fontsize=20)
    #plt.show()
    plt.savefig("scenario_forecasts.png")

def without_acs_uncertainty():
    '''
    Plots the comparision between nitrate concentrations with and without
    an active carbon sink but taking into account uncertainty

    Notes
    ------
    Saves the plot in file 'without_acs_uncertainty.png'
    '''

    t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
    t = np.arange(1980,2020,step = 0.1)

    t_pos = np.arange(1980,2019.75,step = 0.1)

    ps, p = posterior_pars()

    t_forecast = np.arange(2020,2030,step=0.05)

    v=0.3
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    ax.plot(t0,c0, 'ro', label="Data", markersize=2.5)

    for pi in range(0,ps.shape[0]):
        ax.plot(t_pos, LPM_Model_forecast(t_pos, ps[pi][0], 1, ps[pi][2], ps[pi][3],0), 'g-', lw=0.3)


    for pi in range(0,ps.shape[0]):
        ax.plot(t_pos, LPM_Model_forecast(t_pos, ps[pi][0], ps[pi][1], ps[pi][2], ps[pi][3],0), 'k-', lw=0.3)

    ax.plot([],[], 'k-', label = 'With ACS')
    ax.plot([],[], 'g-', label = 'Without ACS')
    ax.legend()
    plt.title("Uncertainty caused by lack of alpha parameter", fontsize=20)
    #plt.show()
    fig.savefig('without_acs_uncertainty.png', dpi = 200)