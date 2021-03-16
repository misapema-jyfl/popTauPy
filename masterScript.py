#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 09:27:03 2021

Master script for running the optimization routines.
This script must be run from the command line, 
where the parameters file (.yaml format) must also be specified.

Example:
    >>> python3 masterScript.py params.yaml

If you wish not to execute a given script, 
set its execution option to false in the parameters file.

@author: miha
"""

import time
import sys
import yaml
from opt_abc import Fitter
from opt_neTe import Optimizer
from opt_neTe import make_biases
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

font = {"family":"normal",
        "weight":"normal",
        "size": 15}

matplotlib.rc("font", **font)


# The path to the parameters file is given in terminal.
paramFilePath = sys.argv[1]
params = yaml.safe_load(open(paramFilePath))

d = params["data"]
g = params["general"]
o = params["optimizer"]


def do_abc():
    # Perform the abc fitting...
    df_res = pd.DataFrame()
    for charge_state in d["available_charge_states"][1:-1]:
        
        print("Working on charge_state: " + str(charge_state) + "+...")
        
        F = Fitter(charge_state=charge_state, params=params)
        
        popt, pcov, chi2_reduced = F.do_fit()
        
        # Plot the final fit
        t,i = F.get_data(charge_state)
        t_i, t_f = F.determine_interpolation_limits()
        T = np.linspace(t_i, t_f, num=len(t))
        y = F.fit_rk4(a=popt[0], b=popt[1], c=popt[2])
        
        fig, ax = plt.subplots()
        
        # s = "{ " + str(charge_state) + "+}"
        
        ax.plot(t*1e3,i*1e3,c="k", label="Measured")
        ax.plot(T*1e3,y(T)*1e3,c="r", ls="--", label="Fitted")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Current (enA)")
        ax.legend()
        outName = "fig_fit_h{:.0e}_q{}".format(o["rk_time_step"], str(charge_state))     
        fig.tight_layout()    
        fig.savefig(g["save_to_path"] + outName + "+.eps", format="eps")
        fig.savefig(g["save_to_path"] + outName + "+.png", format="png", dpi=300)
        plt.close(fig)
    
        # Unpack the result
        [a, b, c] = popt
        [da, db, dc] = np.sqrt(np.diag(pcov))
        
        df_res[charge_state]=[a,b,c,da,db,dc,chi2_reduced]
        
    df_res.index=["a", "b", "c", "da", "db", "dc", "chi2"]
    
    print(df_res)
    
    df_res.to_csv(g["save_to_path"] + d["abc_file_name"])
    


def run_algorithm(charge_state):
    '''
    Run the optimisation routine on a chosen charge state,
    using the runtime parameters separately specified.
    '''
    
    print("Starting run for charge state {}...\n\n".format(str(charge_state)))
    
    q = charge_state # Rename for brevity
    
    # Instantiate the optimizer object for charge state q
    optimizer = Optimizer(q, params=params)
    
    # Create the list of Voronov biases to use
    b = make_biases(optimizeObject=optimizer)
    
    # Track number of iterations
    i=0
    
    # Find the solution set with each given set of Voronov biases
    # and pack the solution set to the output dataframe
    ne = []
    Ee = []
    tau = []
    inz_rate = []
    cx_rate = []
    eC = []
    F = []
    used_biases_l = []
    used_biases = []
    used_biases_h = []
    
    absoluteStart = time.perf_counter() # For tracking elapsed time
    for _ in range(optimizer.number_of_MC_iters):
        
        # Set the Voronov biases for this iteration    
        optimizer.biases[q-1] = b[0][i]
        optimizer.biases[q] = b[1][i]
        optimizer.biases[q+1] = b[2][i]
        
        # Find solution set using above biases
        start = time.perf_counter()
        result = optimizer.find_solution_set()
        finish = time.perf_counter()
        
        i += 1
        
        # Append the output lists
        [ne.append(el) for el in result["ne"]]
        [Ee.append(el) for el in result["Ee"]]
        [tau.append(el) for el in result["tau"]]
        [inz_rate.append(el) for el in result["inz_rate"]]
        [cx_rate.append(el) for el in result["cx_rate"]]
        [F.append(el) for el in result["F"]]
        [eC.append(el) for el in result["eC"]]
        [used_biases_l.append(optimizer.biases[q-1]) for _ in range(len(result["F"]))]
        [used_biases.append(optimizer.biases[q]) for _ in range(len(result["F"]))]
        [used_biases_h.append(optimizer.biases[q+1]) for _ in range(len(result["F"]))]
        
        # Print runtime status...
        numSol = len(ne)
        totSol = i*optimizer.N
        print("Finished iteration {} in {} s (total: {} s) | Cumulative accuracy: {} %"
              .format(i,
                      round(finish-start,1),
                      round(finish-absoluteStart,1),
                      round(100*numSol/totSol,1))
              )
        
    
    # Create the output dataframe
    df_out = pd.DataFrame()
    df_out["n"] = ne
    df_out["E"] = Ee
    df_out["tau"] = tau
    df_out["inz_rate"] = inz_rate
    df_out["cx_rate"] = cx_rate
    df_out["eC"] = eC
    df_out["bias_l"] = used_biases_l
    df_out["bias"] = used_biases
    df_out["bias_h"] = used_biases_h
    df_out["F"] = F
    
    
    
    # Save results to file
    name = "solset_MC_iters-{}_N-{}_q-{}.csv".format(
        optimizer.number_of_MC_iters,
        optimizer.N,
        optimizer.q)
    outputPath = g["save_to_path"] + name
    df_out.to_csv( outputPath )
    
    print("Run completed!\n\n")

def do_neTe():
    # Execute the routine
    cStates = d["available_charge_states"][2:-2]
    for q in cStates:
        run_algorithm(q)
        

if g["do_abc"]:
    do_abc()
if g["do_neTe"]:
    do_neTe()