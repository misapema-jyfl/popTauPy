#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 10:57:47 2020

@author: miha
"""

import numpy as np
import pandas as pd
from scipy.optimize import NonlinearConstraint
from scipy.optimize import minimize
import time
import numdifftools as nd
import concurrent.futures

abc_file_path = "./solution_abc_2020-06-29_w_noise_h=10e-6.csv"
voronov_file_path = "./voronov_k.csv"



#
# Get the abc file into a dataframe and convert column names (charge states)
# to integers.
#
df_abc = pd.read_csv(abc_file_path, index_col=(0), sep=",")
df_abc.columns = [int(name) for name in list(df_abc.columns)]

#
# Get the voronov coefficient file into a dataframe.
#
df_voronov = pd.read_csv(voronov_file_path, index_col="state")


cStates = list(df_abc.columns) # Charge states for which abc are available


def voronov_rate(chargeState, electronTemperature, bias_coefficients):
    '''
    Calculate the q -> q+1 electron-impact ionization rate coefficient
    according to the semi-empirical formula by Voronov.
    
    N.B. Declare the voronov_coefficient_file in the code before this function.
    
    Parameters
    ----------
    chargeState : int
        Charge state of ion in question.
    electronTemperature : float
        Electron temperature for the electron population responsible 
        for ionization.

    Returns
    -------
    rateCoefficient : float
        Rate coefficient of ionization [cm3/s].

    '''
    
    
    q = chargeState
    Te = electronTemperature
    
    i = q-cStates[0]
    
    dE = df_voronov["dE"][q].astype(np.float)
    K = df_voronov["K"][q].astype(np.float)
    A = df_voronov["A"][q].astype(np.float)
    X = df_voronov["X"][q].astype(np.float)
    P = df_voronov["P"][q].astype(np.float)

    U = dE / Te
    rateCoefficient = A*(U**K)*np.exp(-U)*((1+P*U**0.5)/(X + U))*bias_coefficients[i]
    
    return rateCoefficient



def calculate_confinement_time(q, Te, ne, bias_coefficients):

    # Get a,b,c
    a = df_abc[q]["a"]
    b = df_abc[q]["b"]
    c_l = df_abc[q-1]["c"]
    
    # Calculate the rate coefficients
    sv_l = voronov_rate(q-1, Te, bias_coefficients)
    sv = voronov_rate(q, Te, bias_coefficients)
    
    # Calculate (postdiction for) the confinement time
    tau = (b - ne*sv - (a*c_l/(ne*sv_l)))**(-1)

    return tau

def calculate_cx_rate(q, Te, ne, bias_coefficients):
    
    b = df_abc[q]["b"]

    tau = calculate_confinement_time(q, Te, ne, bias_coefficients)
    
    # Ionization rate.
    inz_rate = voronov_rate(q, Te, bias_coefficients)*ne
    
    cx_rate = b - inz_rate - 1/tau
    
    return cx_rate

def F(Te, ne, q, bias_coefficients):
    '''
    Penalty function to minimize.

    Parameters
    ----------
    Te : float
        Electron temperature.
    ne : float
        Electron density.
    q : int
        Charge state.

    Returns
    -------
    F : float
        Penalty function value.

    '''
    
    a_h = df_abc[q+1]["a"]    
    
    # Calculate the (biased) rate coefficient
    sv = voronov_rate(q, Te, bias_coefficients)    
    
    tau_ratio = (q/(q+1))*(a_h/(ne*sv))
    
    tau = calculate_confinement_time(q, Te, ne, bias_coefficients)
    tau_h = calculate_confinement_time(q+1, Te, ne, bias_coefficients)
    
    F = 100*np.abs(tau_ratio - tau/tau_h)/tau_ratio
        
    return F



def set_bias_coefficients(unc_lo, unc_hi, length):
    bias_coefficients = np.random.uniform(unc_lo,unc_hi,len(cStates)) 
    return bias_coefficients


def fill_solution_space(q, ne_lo, ne_hi, Te_lo, Te_hi, MC_unc_lo, MC_unc_hi, N):
    

    #
    # Create Voronov bias coefficients.
    #
    bCoeffs = set_bias_coefficients(unc_lo=MC_unc_lo, unc_hi=MC_unc_hi, length=cStates)
    
    #
    # Create the ne vector, and displace points on vector by small random amount.
    # Make the displacement s.t. each n-point may be displaced at most halfway 
    # toward the next n-point to its left or to its right.
    #
    ns = np.logspace(ne_lo, ne_hi, N)
    randoms = np.random.uniform(low=-.5, high=.5, size=N)
    
    for i in range(len(ns)):
        if i == 0:
            continue
        if i == len(ns) - 1:
            continue
        rand = randoms[i]
        n = ns
        if rand < 0:
            ns[i] = n[i] + rand*(n[i]-n[i-1]) # Displacement towards lower n-point
        if rand >= 0:
            ns[i] = n[i] + rand*(n[i+1]-n[i]) # Displacement towards higher n-point
    
    #
    # Set the constraints.
    #
    if q != max(cStates):
        cons = (NonlinearConstraint(lambda Te: calculate_confinement_time(q, Te, n, bCoeffs), 0, np.inf),
                NonlinearConstraint(lambda Te: calculate_confinement_time(q+1, Te, n, bCoeffs), 0, np.inf),
                NonlinearConstraint(lambda Te: calculate_cx_rate(q, Te, n, bCoeffs), 0, np.inf))
                
    else:
        cons = (
        NonlinearConstraint(lambda Te: calculate_confinement_time(q, Te, n, bCoeffs), 0, np.inf),
        NonlinearConstraint(lambda Te: calculate_cx_rate(q, Te, n, bCoeffs), 0, np.inf)
        )
    
    #
    # Set the bounds
    #
    bnds = [(Te_lo, Te_hi)]
    
    opt_nes = [] # To collect solutions
    opt_Tes = [] 
    opt_Fs = []
    opt_taus = []
    opt_inz_rates = []
    opt_cx_rates = []
    opt_eCs = [] # Energy contents (n*T)
    
    #
    # For each n-point, perform the minimization of the penalty function.
    #
    for n in ns:
        
        
        T0 = np.random.uniform(low=Te_lo, high=Te_hi) #Set initial guess T0
        
        result = minimize(fun=F, 
                          x0=T0,
                          args=(n, q, bCoeffs),
                          method="SLSQP", 
                          constraints=cons, 
                          bounds=bnds,
                          jac=nd.Jacobian(F),
                          options= {"maxiter":25,
                                    "ftol":1e-6})
        Te = result.x[0]
                
        # print(result)
        
        # If minimization was successful, save the result. Else don't.
        if result.success == True:
            
            opt_nes.append(n)
            opt_Tes.append(Te)
            opt_Fs.append(F(Te, n, q, bCoeffs))
            
            # Calculate tau and save
            tau = calculate_confinement_time(q, Te, n, bCoeffs)
            opt_taus.append(tau)
            
            # Calculate ionization rate and save
            i = q-cStates[0]
            inz_rate = voronov_rate(q, Te,bias_coefficients=bCoeffs)*n
            opt_inz_rates.append(inz_rate)
            
            # Calculate charge exchange rate and save
            cx_rate = calculate_cx_rate(q, Te, n, bias_coefficients=bCoeffs)
            opt_cx_rates.append(cx_rate)
            
            # Calculate energy contents and save
            eC = n*Te
            opt_eCs.append(eC)
            
            
    # Make a dictionary of the results
    result = {}
    result["n"] = opt_nes
    result["T"] = opt_Tes
    result["tau"] = opt_taus
    result["F"] = opt_Fs
    result["inz_rate"] = opt_inz_rates
    result["cx_rate"] = opt_cx_rates
    result["energy_content"] = opt_eCs
    
    
    return result # Return the result
    





##############################################################################
# Set the code parameters below.                                             #
##############################################################################
    
qs = [6]
ne_lo = 11
ne_hi = 12.41664051
Te_lo = 10
Te_hi = 10e3
MC_unc_lo = -.6
MC_unc_hi = .6
N = 10
number_of_MC_iters = 1

##############################################################################
# No need to touch anything below this line!                                 #
##############################################################################





def helper(arguments):
    return fun(arguments[0])


def fun(q):
    '''
    Wrapper function for multithreading.
    '''
    result = fill_solution_space(q, ne_lo, ne_hi, Te_lo, Te_hi, MC_unc_lo, MC_unc_hi, N)
    return result

def main():
    '''
    main() function for multithreading.
    '''
    
    for q in qs:        
               
        iterList = [(q, iteration) for iteration in range(number_of_MC_iters)]
        
        opt_nes = [] # To collect solutions
        opt_Tes = [] 
        opt_Fs = []
        opt_taus = []
        opt_inz_rates = []
        opt_cx_rates = []
        opt_eCs = [] # Energy contents (n*T)
        
        with concurrent.futures.ProcessPoolExecutor() as executor:
            
            iteration = 0
            
            t0 = time.time()
            for result in executor.map(helper, iterList):
                t_elapsed = time.time() - t0
                
                iteration += 1
                
                # Send results to collection lists.
                [opt_nes.append(n) for n in result["n"]]
                [opt_Tes.append(Te) for Te in result["T"]]
                [opt_Fs.append(F) for F in result["F"]]
                [opt_taus.append(tau) for tau in result["tau"]]
                [opt_inz_rates.append(inz_rate) for inz_rate in result["inz_rate"]]
                [opt_cx_rates.append(cx_rate) for cx_rate in result["cx_rate"]]
                [opt_eCs.append(eC) for eC in result["energy_content"]]
                    
                cumSols = len(opt_nes) # Calculate number of cumulative solutions
                potSols = N*(iteration) # Number of potential solutions up to this iteration
                
                # Print status.
                print("Iteration {} completed. Time elapsed:  {:.2f} s. Cumulative solutions: {} / {}".format(iteration,
                                                                                                         t_elapsed,
                                                                                                         cumSols,
                                                                                                         potSols))
        result = {} # Make a dictionary of the results
        result["n"] = opt_nes
        result["T"] = opt_Tes
        result["tau"] = opt_taus
        result["F"] = opt_Fs
        result["inz_rate"] = opt_inz_rates
        result["cx_rate"] = opt_cx_rates
        result["energy_content"] = opt_eCs
        
        # Save resultant dictionary to file
        df_result = pd.DataFrame(result)
        df_result.to_csv("./results/unc_lo={:.0f}%_unc_hi={:.0f}%_MC_iters={}_N={}_q={}.csv".format(MC_unc_lo*100,
                                                                                        MC_unc_hi*100,
                                                                                        number_of_MC_iters,
                                                                                        N,
                                                                                        q))


# Code execution.
if __name__ == '__main__':
    main()










