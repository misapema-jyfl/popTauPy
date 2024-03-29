#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 09:27:03 2021

Master script for running the optimization routines.
This script must be run from the command line, 
where the parameters file (.yaml format) must also be specified.

Example:
    >>> python3 masterScript.py params.yaml

Use the GUI to create the parameters file.

@author: miha
"""

import time
import sys
import yaml
from parse_data import Parser
from opt_abc import Fitter
from opt_neTe import Optimizer
from opt_neTe import make_biases
from plot_results import SolSetPlotter
from plot_results import Plotter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from paramtest import ParamTest

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
p = params["parsing"]

if not d["available_charge_states"] == False:
    cStates = d["available_charge_states"]
else:
    cStates = p["available_charge_states"]


def do_abc():
    # Perform the abc fitting...
    df_res = pd.DataFrame()
    for charge_state in cStates[1:-1]:
        
        print("Working on charge_state: " + str(charge_state) + "+...")
        
        F = Fitter(charge_state=charge_state, params=params)
        
        popt, pcov, chi2_reduced = F.do_fit()
        
        
        t,i = F.get_data(charge_state)
        t_i, t_f = F.determine_interpolation_limits()
        T = np.linspace(t_i, t_f, num=len(t))
        y = F.fit_rk4(a=popt[0], b=popt[1], c=popt[2])
        
        # Save the fit in a DataFrame
        df = pd.DataFrame()
        df["t"] = T
        df["i"] = y(T)
        
        # Output the fit data to .csv
        outputFileName = "out_fit_h{:.0e}_q{}.csv".format(o["rk_time_step"], str(charge_state))
        df.to_csv(g["save_to_path"] + outputFileName)
        
        # Plot the final fit
        fig, ax = plt.subplots()
        
        # s = "{ " + str(charge_state) + "+}"
        
        ax.plot(t*1e3,i,c="k", label="Measured")
        ax.plot(T*1e3,y(T),c="r", ls="--", label="Fitted")
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel(r"Current ($e$A)")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
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
    Run the optimisation routine on chosen charge state,
    using the runtime parameters separately specified.
    '''
    
    print("Starting run for charge state {}...\n\n".format(str(charge_state)))
    
    q = charge_state # Rename for brevity
    
    # Instantiate the optimizer object for charge state q
    optimizer = Optimizer(q, params=params)
    
    # Create the list of uncertainty biases to use
    b = make_biases(optimizeObject=optimizer)
    
    # Track number of iterations
    i=0
    
    # Find the solution set with each given set of biases
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
    charge_states = cStates[2:-2]
    for q in charge_states:
        run_algorithm(q)
        


# Test given parameters file for errors.
pt = ParamTest(params)
errors, warnings = pt.do_tests()

# If no errors, run the script.
if errors == 0:
    if g["do_parse"]:
        print("\nParsing data...")
        P = Parser(params)
        P.do_parse()
    if g["do_abc"]:
        print("\nDoing a,b,c fitting...")
        do_abc()
    if g["do_neTe"]:
        print("\nDetermining solution set for given a,b,c values...")
        do_neTe()
    if g["do_plots"]:
        
        print("\nPlotting results...")
        
        if params["plotting"]["plot_solution_sets"]["do_plot"]:
            SSP = SolSetPlotter(params)
            SSP.do_solution_set_plots()
        
        P = Plotter(params)
        P.do_result_plots()
        
    
print("\nProcess finished!\n")