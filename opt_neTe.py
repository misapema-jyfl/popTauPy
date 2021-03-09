#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:58:53 2021

@author: miha
"""

import time
import numpy as np
import concurrent.futures
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from parameters import p
from scipy.optimize import NonlinearConstraint
# import numdifftools as nd

class Optimizer:
    
    def __init__(self, q):
        
        # Set current charge state 
        self.q = q
        
        # Retrieve parameters for runtime
        self.abc_file_path = p["abc_file_path"]
        self.voronov_file_path = p["voronov_file_path"]
        self.cStates = p["cStates"]
        self.ne_lo = p["ne_lo"]
        self.ne_hi = p["ne_hi"]
        self.Te_lo = p["Te_lo"]
        self.Te_hi = p["Te_hi"]
        self.MC_unc_lo = p["MC_unc_lo"]
        self.MC_unc_hi = p["MC_unc_hi"]
        self.N = p["N"]
        self.number_of_MC_iters = p["number_of_MC_iters"]

        # Set initial Voronov biases
        self.biases = {}
        self.biases[q-1] = 1
        self.biases[q] = 1
        self.biases[q+1] = 1
        
        
        # Get the abc file into a dataframe and convert column names 
        # (charge states) to integers.
        self.df_abc = pd.read_csv(self.abc_file_path, index_col=(0), sep=",")
        self.df_abc.columns = [int(name) for name in list(self.df_abc.columns)]


        # Get the voronov coefficient file into a dataframe.
        self.df_voronov = pd.read_csv(self.voronov_file_path, index_col="state")



    def voronov_rate(self, chargeState, electronTemperature):
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
    
        # Rename for brevity
        q = chargeState
        Te = electronTemperature
        
        # Retrieve formula parameters for charge state q         
        dE = self.df_voronov["dE"][q].astype(np.float)
        K = self.df_voronov["K"][q].astype(np.float)
        A = self.df_voronov["A"][q].astype(np.float)
        X = self.df_voronov["X"][q].astype(np.float)
        P = self.df_voronov["P"][q].astype(np.float)
    
        # Calculate the rate coefficient
        U = dE / Te
        
        rateCoefficient = A*(U**K)*np.exp(-U)*((1+P*U**0.5)/(X + U))
        
        return rateCoefficient


    def calculate_confinement_time(self, q, Te, ne):
        '''
        Calculate the confinement time according to 
        Eq. () in the manuscript.
        
        q = charge state
        Te = electron temperature (eV)
        ne = electron density (1/cm3)
        '''
        
        # Get a,b,c
        a = self.df_abc[q]["a"]
        b = self.df_abc[q]["b"]
        c_l = self.df_abc[q-1]["c"]
        
        # Calculate the rate coefficients
        sv_l = self.voronov_rate(q-1, Te)*self.biases[q-1]
        sv = self.voronov_rate(q, Te)*self.biases[q]
        
        # Calculate (postdiction for) the confinement time
        tau = (b - ne*sv - (a*c_l/(ne*sv_l)))**(-1)
    
        return tau


    def calculate_cx_rate(self, q, Te, ne):
        '''
        Calculate the charge exchange rate.
        
        q = charge state
        Te = electron temperature (eV)
        ne = electron density (1/cm3)
        '''
        b = self.df_abc[q]["b"]
    
        tau = self.calculate_confinement_time(q, Te, ne)
        
        inz_rate = self.biases[q]*self.voronov_rate(q, Te)*ne
        
        cx_rate = b - inz_rate - 1/tau
        
        return cx_rate




    def F(self, Te, ne, q):
        '''
        Penalty function to minimize.
        
        Te = electron temperature (eV)
        ne = electron density (1/cm3)
        q = charge state
        '''
        
        a_h = self.df_abc[q+1]["a"]    
        
        # Calculate the (biased) rate coefficient
        sv = self.voronov_rate(q, Te)*self.biases[q] 
        
        tau_ratio = (q/(q+1))*(a_h/(ne*sv))
        tau = self.calculate_confinement_time(q, Te, ne)
        tau_h = self.calculate_confinement_time(q+1, Te, ne)
        
        F = 100*np.abs(tau_ratio - tau/tau_h)/tau_ratio
            
        return F



    
        
    def make_constraints(self, n):
        '''
        Create constraints for the minimize_F algorithm
        '''
        
        q = self.q
        
        cons = (NonlinearConstraint(lambda Te: self.calculate_confinement_time(q, Te, n), 0, np.inf),
                    NonlinearConstraint(lambda Te: self.calculate_confinement_time(q+1, Te, n), 0, np.inf),
                    NonlinearConstraint(lambda Te: self.calculate_cx_rate(q, Te, n), 0, np.inf))
    
        return cons
    
    
    
    def minimize_F(self, n):
        '''
        Minimize the penalty function as a function of Te 
        at every given ne value.
        '''        
        
        # Make the return dictionary
        res = {}
        
        # Make the initial guess for Te
        initialGuessForTe = np.random.uniform(low=self.Te_lo, high=self.Te_hi)
        
        
        bnds = [(self.Te_lo, self.Te_hi)] # Bounded in Te
        
        
        # Minimize F
        result = minimize(fun=self.F,
                 x0=initialGuessForTe,
                 args=(n, self.q),
                 method="SLSQP",
                 bounds=bnds
                 ) 
                
        # Check whether minimization was successful
        if result.success == True:
            
            Te = result.x[0]
            
            # Calculate values of F, cx, inz, tau
            Fval = self.F(Te, n, self.q)
            tau = self.calculate_confinement_time(self.q, Te, n)
            inz_rate = self.biases[self.q]*self.voronov_rate(self.q, Te)*n
            cx_rate = self.calculate_cx_rate(self.q, Te, n)
            eC = n*Te
            tau_h = self.calculate_confinement_time(self.q + 1, Te, n)
            
            # Make sure that none of the values is unphysical
            if tau > 0 and inz_rate > 0 and cx_rate > 0 and eC > 0 and tau_h > 0:
            
                res["success"] = True
                res["tau"] = tau
                res["inz_rate"] = inz_rate
                res["cx_rate"] = cx_rate
                res["eC"] = eC
                res["F"] = Fval
                res["ne"] = n
                res["Te"] = Te
                
                return res
            else:
                res["success"] = False
        else:
            res["success"] = False
    
        return res
    

    

    def create_ne_array(self):
        '''
        Create the electron density array.
        
        By randomly displacing the elements created in a logspaced array,
        we can fill the ne axis over the course of the ~1000 
        iterations normally performed during the algorithm runtime.
        '''
        ns = np.logspace(self.ne_lo, self.ne_hi, num=self.N)
        randoms = np.random.uniform(low=-.5, high=.5, size=self.N)
        
        # Displace each element by a random amount towards 
        # the previous, or subsequent element in the vector.
        for i in range(len(ns)):
            
            if i == 0:
                # Don't displace the first element.
                continue 
            
            if i == len(ns) - 1:
                # Don't displace the last element.
                continue
            
            rand = randoms[i]
            n = ns
            
            if rand < 0:
                # Displacement towards previous element
                ns[i] = n[i] + rand*(n[i]-n[i-1]) 
            
            if rand >= 0:
                # Displacement towards next element
                ns[i] = n[i] + rand*(n[i+1]-n[i]) 
                
        n = ns
        
        return n




    def find_solution_set(self):
        '''
        Find the set of acceptable solutions.
        Call the 'minimize_F' method in each element
        of the array 'n' which spans the given n_e limits.
        '''
        
        # Create lists to hold the results
        optimized_nes = []
        optimized_Tes = []
        optimized_Fs = []
        optimized_taus = []
        optimized_inz_rates = []
        optimized_cx_rates = []
        optimized_energy_contents = []
        
        # Use multiprocessing
        if __name__ == "__main__":
            with concurrent.futures.ProcessPoolExecutor() as executor:
                
                # Create the ne array
                n = self.create_ne_array()
                
                # Run minimisation at each n
                for result in executor.map(self.minimize_F, n):
                    
                    # Check that there exists a solution
                    if result["success"] == True:
                        
                        # Get the values at (n, T)
                        n = result["ne"]
                        Te = result["Te"]
                        Fval = result["F"]
                        tau = result["tau"]
                        inz_rate = result["inz_rate"]
                        cx_rate = result["cx_rate"]
                        eC = result["eC"]
                        
                        # Append values to results
                        optimized_nes.append(n)
                        optimized_Tes.append(Te)
                        optimized_Fs.append(Fval)            
                        optimized_taus.append(tau)
                        optimized_inz_rates.append(inz_rate)
                        optimized_cx_rates.append(cx_rate)
                        optimized_energy_contents.append(eC)
        
        
        # Pack results to dictionary
        results = {}
        results["ne"] = optimized_nes
        results["Te"] = optimized_Tes
        results["F"] = optimized_Fs
        results["tau"] = optimized_taus
        results["inz_rate"] = optimized_inz_rates
        results["cx_rate"] = optimized_cx_rates
        results["eC"] = optimized_energy_contents
        
        return results
            




def make_voronov_biases(optimizeObject):
        '''
        Create a list of biases on the Voronov formula coefficients
        using a random, uniform selection for the biases.
        '''
        opt = optimizeObject
        
        lo = opt.MC_unc_lo
        hi = opt.MC_unc_hi
        Num = opt.number_of_MC_iters
        
        biases_l = np.random.uniform(low=lo, high=hi, size=Num) + 1
        biases= np.random.uniform(low=lo, high=hi, size=Num) + 1 
        biases_h = np.random.uniform(low=lo, high=hi, size=Num) + 1
        
        return [biases_l, biases, biases_h]







def run_algorithm(charge_state):
    '''
    Run the optimisation routine on a chosen charge state,
    using the runtime parameters designated in 
    the settings file 'parameters.py'
    '''
    
    print("Starting run for charge state {}...\n\n".format(str(charge_state)))
    
    q = charge_state # Rename for brevity
    
    # Instantiate the optimizer object for charge state q
    o = Optimizer(q)
    
    # Create the list of Voronov biases to use
    b = make_voronov_biases(optimizeObject=o)
    
    # Track number of iterations
    i=0
    
    # Find the solution set with each given set of Voronov biases
    # and pack the solution set to the output dataframe
    ne = []
    Te = []
    tau = []
    inz_rate = []
    cx_rate = []
    eC = []
    F = []
    used_biases_l = []
    used_biases = []
    used_biases_h = []
    
    absoluteStart = time.perf_counter() # For tracking elapsed time
    for _ in range(o.number_of_MC_iters):
        
        # Set the Voronov biases for this iteration    
        o.biases[q-1] = b[0][i]
        o.biases[q] = b[1][i]
        o.biases[q+1] = b[2][i]
        
        # Find solution set using above biases
        start = time.perf_counter()
        result = o.find_solution_set()
        finish = time.perf_counter()
        
        i += 1
        
        # Append the output lists
        [ne.append(el) for el in result["ne"]]
        [Te.append(el) for el in result["Te"]]
        [tau.append(el) for el in result["tau"]]
        [inz_rate.append(el) for el in result["inz_rate"]]
        [cx_rate.append(el) for el in result["cx_rate"]]
        [F.append(el) for el in result["F"]]
        [eC.append(el) for el in result["eC"]]
        [used_biases_l.append(o.biases[q-1]) for _ in range(len(result["F"]))]
        [used_biases.append(o.biases[q]) for _ in range(len(result["F"]))]
        [used_biases_h.append(o.biases[q+1]) for _ in range(len(result["F"]))]
        
        # Print runtime status...
        numSol = len(ne)
        totSol = i*o.N
        print("Finished iteration {} in {} s (total: {} s) | Cumulative accuracy: {} %"
              .format(i,
                      round(finish-start,1),
                      round(finish-absoluteStart,1),
                      round(100*numSol/totSol,1))
              )
        
    
    # Create the output dataframe
    df_out = pd.DataFrame()
    df_out["n"] = ne
    df_out["T"] = Te
    df_out["tau"] = tau
    df_out["inz_rate"] = inz_rate
    df_out["cx_rate"] = cx_rate
    df_out["eC"] = eC
    df_out["bias_l"] = used_biases_l
    df_out["bias"] = used_biases
    df_out["bias_h"] = used_biases_h
    df_out["F"] = F
    
    
    
    # Save results to file
    name = "unc_lo={:.0f}%_unc_hi={:.0f}%_MC_iters={}_N={}_q={}.csv".format(
        o.MC_unc_lo*100,
        o.MC_unc_hi*100,
        o.number_of_MC_iters,
        o.N,
        o.q)
    outputPath = p["output_directory"] + name
    df_out.to_csv( outputPath )
    
    print("Run completed!\n\n")
  


# Execute the routine
cStates = p["cStates"]
for q in cStates:
    run_algorithm(q)
