import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from model_functions import *

def error_bar_plot():
    #reading data
    t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
    #t = np.arange(1980,2020,step = 0.5)

    #using curve_fit() to find best parameters to input
    pars = curve_fit(LMP_Model,t0,c0,[1,1,1,1,15])
    b = pars[0][0]
    b1 = pars[0][1]
    alpha = pars[0][2]
    bc = pars[0][3]
    tau = pars[0][4]
    C = LMP_Model(t0,b,b1,alpha,bc,tau)
    #setting variance
    v = 0.5
    #calculate the sum-of-squares objective function 
    S = np.sum(np.square(np.subtract(c0, C)/v))

    # plotting commands
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.plot(t0,C,'b-', label='model')
    ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data')
    ax.set_xlabel('time')
    ax.set_ylabel('pressure')
    ax.set_title('objective function: S={:3.2f}'.format(S))
    ax.legend()
    plt.show()

def grid_search(N, pars):
    b = np.linspace(pars[0][0]/2,pars[0][0]*1.5, N)
    # b1 = np.linspace(pars[0][1]/2,pars[0][1]*1.5, N)
    # alpha = np.linspace(pars[0][2]/2,pars[0][2]*1.5, N)
    bc = np.linspace(pars[0][3]/2,pars[0][3]*1.5, N)
    tau = np.linspace(pars[0][4]/2,pars[0][4]*1.5, N)

    # grid of parameter values: returns every possible combination of parameters in a and b
    A, B, C= np.meshgrid(b, bc, tau, indexing='ij')

    # empty 2D matrix for objective function
    S = np.zeros(A.shape)
    
    #reading data
    t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T

    v = 0.5

    # for i in range(len(b)):
    #     for j in range(len(b1)):
    #         # 2. compute the sum of squares objective function at each value 
    #         for k in range(len(alpha)):
    #             for l in range(len(bc)):
    #                 for m in range(len(tau)):
    #                     Cm = LMP_Model(t0, b[i], b1[j], alpha[k],bc[l], tau[m])
    #                     S[i,j,k,l,m] = np.sum(np.square(np.subtract(Cm, c0)/v))

    for i in range(len(b)):
        for j in range(len(bc)):
            for k in range(len(tau)):
                Cm = LMP_Model(t0, b[i], 0.5, 0.1, bc[j], tau[k])
                S[i,j,k] = np.sum(np.square(np.subtract(Cm, c0)/v))

    # 3. compute the posterior
    P = np.exp(-S/2)

    Pint = np.sum(P)*(b[1]-b[0])*(bc[1]-bc[0])*(tau[1]-tau[0])

    P = np.true_divide(P, Pint)
    return b, bc, tau, P

def get_samples(pars, P, N_samples):
    A, B, C = np.meshgrid(pars[0][0], pars[0][3], pars[0][4], indexing='ij')

    mean, covariance = fit_mvn([A,B,C], P)
    #mean = np.average(P)

    # 1. create samples using numpy function multivariate_normal (Google it)
    samples = np.random.multivariate_normal(mean, covariance, size = N_samples)
    
    return samples


def fit_mvn(parspace, dist):
    """Finds the parameters of a multivariate normal distribution that best fits the data

    Parameters:
    -----------
        parspace : array-like
            list of meshgrid arrays spanning parameter space
        dist : array-like 
            PDF over parameter space
    Returns:
    --------
        mean : array-like
            distribution mean
        cov : array-like
            covariance matrix		
    """
    
    # dimensionality of parameter space
    N = len(parspace)
    
    # flatten arrays
    parspace = [p.flatten() for p in parspace]
    dist = dist.flatten()
    
    # compute means
    mean = [np.sum(dist*par)/np.sum(dist) for par in parspace]
    
    # compute covariance matrix
        # empty matrix
    cov = np.zeros((N,N))
        # loop over rows
    for i in range(0,N):
            # loop over upper triangle
        for j in range(i,N):
                # compute covariance
            cov[i,j] = np.sum(dist*(parspace[i] - mean[i])*(parspace[j] - mean[j]))/np.sum(dist)
                # assign to lower triangle
            if i != j: cov[j,i] = cov[i,j]
            
    return np.array(mean), np.array(cov)

def ensemble(samples):
    t = np.arange(1980,2020,step = 0.5)
    v = 0.5

    # create a figure and axes (see TASK 1)
    fig = plt.figure(figsize=(20,10))
    ax = fig.add_subplot(111)
    
    # for each sample, solve and plot the model  (see TASK 1)
    for a,b,c in samples:
        pm = LMP_Model(t,a, 0.5, 0.1, b,c)
        ax.plot(t,pm,'k-', lw = 0.5, alpha = 0.3)

    # this command just adds a line to the legend
    #ax.plot([],[],'k-', lw=0.5,alpha=0.4, label='model ensemble')

    # get the data
    t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T	
    ax.axvline(1990, color='b', linestyle=':', label='calibration')
    
    # 4. plot Wairakei data as error bars
    # *hint* see TASK 1 for appropriate plotting commands
    ax.errorbar(t0,c0,yerr=v,fmt='ro', label='data')
    ax.set_xlabel('time')
    ax.set_ylabel('nitrate concentration')
    ax.legend()
    plt.show()
    #plt.savefig('task6_3.png')

#error_bar_plot()
t0, c0 = np.genfromtxt('nl_n.csv', delimiter=',', skip_header=1).T
#t = np.arange(1980,2020,step = 0.5)

#using curve_fit() to find best parameters to input
pars = curve_fit(LMP_Model,t0,c0,[1,1,1,1,15])
b, bc, tau, P = grid_search(3, pars)

samples = get_samples(pars, P, 3)

ensemble(samples)