#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:58:53 2021

@author: miha
"""

# import time
import numpy as np
import concurrent.futures
import pandas as pd
from scipy.optimize import minimize
# from parameters import p
# from parameters import general_parameters as gp




class Optimizer:
    
    def __init__(self, q, params):
        
        g = params["general"]
        d = params["data"]
        o = params["optimizer"]
        p = params["parsing"]
        
        # Set the directory structure
        # making sure that the given paths point to directories
        # i.e. there is a "/" at the end of the string.
        self.workDir = g["working_directory"] if g["working_directory"][-1]=="/" else g["working_directory"] + "/"
        self.elementalDir = g["elemental_data_directory"] if g["elemental_data_directory"][-1] == "/" else g["elemental_data_directory"] + "/"
        
        
        # Set current charge state 
        self.q = q
        
        # Retrieve parameters for runtime
        self.abc_file_path = g["save_to_path"] + d["abc_file_name"]
        self.method = o["rate_coefficient_method"]
        self.species = d["species"].lower()
        
        if not d["available_charge_states"] == False:
            self.cStates = d["available_charge_states"][1:-1]
        else:
            self.cStates = p["available_charge_states"][1:-1]
            
        self.ne_lo = np.log10(o["ne_lo"])
        self.ne_hi = np.log10(o["ne_hi"])
        self.Ee_lo = o["Ee_lo"]
        self.Ee_hi = o["Ee_hi"]
        self.N = o["number_of_ne"]
        self.number_of_MC_iters = o["number_of_MC"]
        
        # Retrieve the elemental data parameters corresponding to 
        # the chosen method for evaluating the rate coefficient.
        # Retrieve the function by which to evaluate the rate coefficient.
        if self.method=="voronov":
            s = (self.elementalDir,self.method,"_",self.species,".csv")
            filePath = "".join(s)
            self.df_eldata = pd.read_csv(filePath, index_col="state")
        
        if self.method=="interpMB":
            s = (self.elementalDir,self.method,"_",self.species,".csv")
            filePath = "".join(s)
            self.df_eldata = pd.read_csv(filePath, index_col="state")
        
        # Set the function by which to evaluate the rate coefficients.
        self.rc = self.evaluate_rc()
            
        # Set initial uncertainty biases
        self.biases = {}
        self.biases[q-1] = 1
        self.biases[q] = 1
        self.biases[q+1] = 1
        
        # Retrieve the uncertainty bounds from the data file.
        self.MC_unc_lo = self.df_eldata["UNC"][q-1]/100
        self.MC_unc = self.df_eldata["UNC"][q]/100
        self.MC_unc_hi = self.df_eldata["UNC"][q+1]/100
        
        
        # Get the abc file into a dataframe and convert column names 
        # (charge states) to integers.
        self.df_abc = pd.read_csv(self.abc_file_path, index_col=(0), sep=",")
        self.df_abc.columns = [int(name) for name in list(self.df_abc.columns)]


        
    
    def evaluate_rc(self):
        '''
        Evaluate the rate coefficient using the chosen method.
        The method is specified in the 'parameters.py'
        
        Returns
        -------
        rc : callable
            Chosen function for rate coefficient evaluation.
            N.B. The callable must be a function of the 
            ion charge state and the average energy of the distribution!
        '''
        
        if self.method=="voronov":
            rc = self.voronov_rate
            
        if self.method=="interpMB":
            
            # Load the three necessary interpolation functions
            # into a dataframe object.
            # This is necessary for the actual rate coefficient function.
            self.df_interpFuns = pd.DataFrame()
            self.df_interpFuns[self.q-1] = [self.load_interp_fun(self.q-1, self.method)]
            self.df_interpFuns[self.q] = [self.load_interp_fun(self.q, self.method)]
            self.df_interpFuns[self.q+1] = [self.load_interp_fun(self.q+1, self.method)]
            
            rc = self.interpolated_MB
            
        return rc


    def load_interp_fun(self, chargeState, interpType):
        '''
        Returns the interpolation function of type interpType
        (if exists) for the ion of the given charge state.
        '''
        s = (self.elementalDir,interpType,"_",self.species,"_",
             str(chargeState),"+.npy")
        filename = "".join(s)
        interpFun = np.load(filename,allow_pickle=True).item()
        
        return interpFun


    def interpolated_MB(self, chargeState, averageEnergy):
        
        # Pick the correct interpolation function from the dataframe
        # and then calculate the rate coefficient at given <E>.
        f = self.df_interpFuns[chargeState].values[0]
        rateCoefficient = f(averageEnergy)*1e6 # Convert m3 -> cm3 
        
        return rateCoefficient



    def voronov_rate(self, chargeState, averageEnergy):
        '''
        Calculate the q -> q+1 electron-impact ionization rate coefficient
        according to the semi-empirical formula by Voronov.
        
        N.B. Declare the voronov_coefficient_file in the code before this function.
        
        Parameters
        ----------
        chargeState : int
            Charge state of ion in question.
        averageEnergy : float
            Average electron energy for the electron population responsible 
            for ionization.
    
        Returns
        -------
        rateCoefficient : float
            Rate coefficient of ionization [cm3/s].
    
        '''
    
        # Rename for brevity
        q = chargeState
        Te = averageEnergy*(2/3) # convert <Ee> to Te
        
        # Retrieve formula parameters for charge state q         
        dE = self.df_eldata["dE"][q].astype(np.float)
        K = self.df_eldata["K"][q].astype(np.float)
        A = self.df_eldata["A"][q].astype(np.float)
        X = self.df_eldata["X"][q].astype(np.float)
        P = self.df_eldata["P"][q].astype(np.float)
    
        # Calculate the rate coefficient
        U = dE / Te
        
        rateCoefficient = A*(U**K)*np.exp(-U)*((1+P*U**0.5)/(X + U))
        
        return rateCoefficient


    def calculate_confinement_time(self, q, Ee, ne):
        '''
        Calculate the confinement time according to 
        Eq. () in the manuscript.
        
        q = charge state
        Ee = Average electron energy (eV)
        ne = electron density (1/cm3)
        '''
        
        # Get a,b,c
        a = self.df_abc[q]["a"]
        b = self.df_abc[q]["b"]
        c_l = self.df_abc[q-1]["c"]
        
        # Calculate the rate coefficients
        sv_l = self.rc(q-1, Ee)*self.biases[q-1]
        sv = self.rc(q, Ee)*self.biases[q]
        
        # Calculate (postdiction for) the confinement time
        tau = (b - ne*sv - (a*c_l/(ne*sv_l)))**(-1)
    
        return tau


    def calculate_cx_rate(self, q, Ee, ne):
        '''
        Calculate the charge exchange rate.
        
        q = charge state
        Ee = Average electron energy (eV)
        ne = electron density (1/cm3)
        '''
        b = self.df_abc[q]["b"]
    
        tau = self.calculate_confinement_time(q, Ee, ne)
        
        inz_rate = self.biases[q]*self.rc(q, Ee)*ne
        
        cx_rate = b - inz_rate - 1/tau
        
        return cx_rate




    def F(self, Ee, ne, q):
        '''
        Penalty function to minimize.
        
        Ee = Average electron energy (eV)
        ne = electron density (1/cm3)
        q = charge state
        '''
        
        a_h = self.df_abc[q+1]["a"]    
        
        # Calculate the (biased) rate coefficient
        sv = self.rc(q, Ee)*self.biases[q] 
        
        tau_ratio = (q/(q+1))*(a_h/(ne*sv))
        tau = self.calculate_confinement_time(q, Ee, ne)
        tau_h = self.calculate_confinement_time(q+1, Ee, ne)
        
        F = 100*np.abs(tau_ratio - tau/tau_h)/tau_ratio
            
        return F

        
    # def make_constraints(self, n):
    #     '''
    #     Create constraints for the minimize_F algorithm
    #     '''
        
    #     q = self.q
        
    #     cons = (NonlinearConstraint(lambda Ee: self.calculate_confinement_time(q, Ee, n), 0, np.inf),
    #                 NonlinearConstraint(lambda Ee: self.calculate_confinement_time(q+1, Ee, n), 0, np.inf),
    #                 NonlinearConstraint(lambda Ee: self.calculate_cx_rate(q, Ee, n), 0, np.inf))
    
    #     return cons
    
    
    
    def minimize_F(self, n):
        '''
        Minimize the penalty function as a function of <Ee> 
        at every given ne value.
        '''        
        
        # Make the return dictionary
        res = {}
        
        # Make the initial guess for <Ee>
        initialGuessForEe = np.random.uniform(low=self.Ee_lo, high=self.Ee_hi)
        
        
        bnds = [(self.Ee_lo, self.Ee_hi)] # Bounded in <Ee>
        
        
        # Minimize F
        result = minimize(fun=self.F,
                 x0=initialGuessForEe,
                 args=(n, self.q),
                 method="SLSQP",
                 bounds=bnds
                 ) 
                
        # Check whether minimization was successful
        if result.success == True:
            
            Ee = result.x[0]
            
            # Calculate values of F, cx, inz, tau
            Fval = self.F(Ee, n, self.q)
            tau = self.calculate_confinement_time(self.q, Ee, n)
            inz_rate = self.biases[self.q]*self.rc(self.q, Ee)*n
            cx_rate = self.calculate_cx_rate(self.q, Ee, n)
            eC = n*Ee
            tau_h = self.calculate_confinement_time(self.q + 1, Ee, n)
            
            # Make sure that none of the values is unphysical
            if tau > 0 and inz_rate > 0 and cx_rate > 0 and eC > 0 and tau_h > 0:
            
                res["success"] = True
                res["tau"] = tau
                res["inz_rate"] = inz_rate
                res["cx_rate"] = cx_rate
                res["eC"] = eC
                res["F"] = Fval
                res["ne"] = n
                res["Ee"] = Ee
                
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
        optimized_Ees = []
        optimized_Fs = []
        optimized_taus = []
        optimized_inz_rates = []
        optimized_cx_rates = []
        optimized_energy_contents = []

        # # Use multiprocessing
        # if __name__ == "__main__":
        with concurrent.futures.ProcessPoolExecutor() as executor:
            
            # Create the ne array
            n = self.create_ne_array()
            
            # Run minimisation at each n
            for result in executor.map(self.minimize_F, n):
                
                # Check that there exists a solution
                if result["success"] == True:
                    
                    # Get the values at (n, T)
                    n = result["ne"]
                    Ee = result["Ee"]
                    Fval = result["F"]
                    tau = result["tau"]
                    inz_rate = result["inz_rate"]
                    cx_rate = result["cx_rate"]
                    eC = result["eC"]
                    
                    # Append values to results
                    optimized_nes.append(n)
                    optimized_Ees.append(Ee)
                    optimized_Fs.append(Fval)            
                    optimized_taus.append(tau)
                    optimized_inz_rates.append(inz_rate)
                    optimized_cx_rates.append(cx_rate)
                    optimized_energy_contents.append(eC)
        
        
        # Pack results to dictionary
        results = {}
        results["ne"] = optimized_nes
        results["Ee"] = optimized_Ees
        results["F"] = optimized_Fs
        results["tau"] = optimized_taus
        results["inz_rate"] = optimized_inz_rates
        results["cx_rate"] = optimized_cx_rates
        results["eC"] = optimized_energy_contents
        
        return results
            




def make_biases(optimizeObject):
        '''
        Create a list of biases on the Voronov formula coefficients
        using a random, uniform selection for the biases.
        '''
        opt = optimizeObject
        
        l = opt.MC_unc_lo
        m = opt.MC_unc
        h = opt.MC_unc_hi
        
        Num = opt.number_of_MC_iters
        
        biases_l = np.random.uniform(low=-l, high=l, size=Num) + 1
        biases = np.random.uniform(low=-m, high=m, size=Num) + 1 
        biases_h = np.random.uniform(low=-h, high=h, size=Num) + 1
        
        return [biases_l, biases, biases_h]







# def run_algorithm(charge_state):
#     '''
#     Run the optimisation routine on a chosen charge state,
#     using the runtime parameters designated in 
#     the settings file 'parameters.py'
#     '''
    
#     print("Starting run for charge state {}...\n\n".format(str(charge_state)))
    
#     q = charge_state # Rename for brevity
    
#     # Instantiate the optimizer object for charge state q
#     o = Optimizer(q)
    
#     # Create the list of Voronov biases to use
#     b = make_biases(optimizeObject=o)
    
#     # Track number of iterations
#     i=0
    
#     # Find the solution set with each given set of Voronov biases
#     # and pack the solution set to the output dataframe
#     ne = []
#     Ee = []
#     tau = []
#     inz_rate = []
#     cx_rate = []
#     eC = []
#     F = []
#     used_biases_l = []
#     used_biases = []
#     used_biases_h = []
    
#     absoluteStart = time.perf_counter() # For tracking elapsed time
#     for _ in range(o.number_of_MC_iters):
        
#         # Set the Voronov biases for this iteration    
#         o.biases[q-1] = b[0][i]
#         o.biases[q] = b[1][i]
#         o.biases[q+1] = b[2][i]
        
#         # Find solution set using above biases
#         start = time.perf_counter()
#         result = o.find_solution_set()
#         finish = time.perf_counter()
        
#         i += 1
        
#         # Append the output lists
#         [ne.append(el) for el in result["ne"]]
#         [Ee.append(el) for el in result["Ee"]]
#         [tau.append(el) for el in result["tau"]]
#         [inz_rate.append(el) for el in result["inz_rate"]]
#         [cx_rate.append(el) for el in result["cx_rate"]]
#         [F.append(el) for el in result["F"]]
#         [eC.append(el) for el in result["eC"]]
#         [used_biases_l.append(o.biases[q-1]) for _ in range(len(result["F"]))]
#         [used_biases.append(o.biases[q]) for _ in range(len(result["F"]))]
#         [used_biases_h.append(o.biases[q+1]) for _ in range(len(result["F"]))]
        
#         # Print runtime status...
#         numSol = len(ne)
#         totSol = i*o.N
#         print("Finished iteration {} in {} s (total: {} s) | Cumulative accuracy: {} %"
#               .format(i,
#                       round(finish-start,1),
#                       round(finish-absoluteStart,1),
#                       round(100*numSol/totSol,1))
#               )
        
    
#     # Create the output dataframe
#     df_out = pd.DataFrame()
#     df_out["n"] = ne
#     df_out["E"] = Ee
#     df_out["tau"] = tau
#     df_out["inz_rate"] = inz_rate
#     df_out["cx_rate"] = cx_rate
#     df_out["eC"] = eC
#     df_out["bias_l"] = used_biases_l
#     df_out["bias"] = used_biases
#     df_out["bias_h"] = used_biases_h
#     df_out["F"] = F
    
    
    
#     # Save results to file
#     name = "solset_MC_iters-{}_N-{}_q-{}.csv".format(
#         o.number_of_MC_iters,
#         o.N,
#         o.q)
#     outputPath = p["output_directory"] + name
#     df_out.to_csv( outputPath )
    
#     print("Run completed!\n\n")
  


# # Execute the routine
# cStates = p["cStates"]
# for q in cStates:
#     run_algorithm(q)
