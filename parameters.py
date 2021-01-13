# -*- coding: utf-8 -*-

import numpy as np

'''
Parameters for running parse_data.py
'''

parsing_parameters = {
    
    "charge_states":[1,2,3,4,5,6,7,8,9,10,11,12],
    
    "1+_control_signals" : ["/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K1+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K2+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K3+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K4+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K5+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K6+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K7+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K8+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K9+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K10+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K11+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/1+_K12+.csv"],
    
    "n+_signals" : ["/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K1+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K2+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K3+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K4+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K5+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K6+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K7+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K8+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K9+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K10+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K11+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/footer_removed_data/K12+.csv"],
    
    "save_to_path" : "/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data",
    
    "elemental_symbol" : "K"
    }


'''
Parameters for running the optimisation routine 'opt_abc.py'.
'''
# charge_states = list of available charge states (in numerical order)
# 1+ control signals = list of data files corresponding to the 1+ control signals
# parsed_data_files = list of the parsed data files (t=0 set the same, b.g. reduced, in numerical order)
# h = Runge-Kutta algorithm time step (units of s)
# output_file_name = file name under which to save the obtained abc parameters  
d = {
     
     "charge_states" : [1,2,3,4,5,6,7,8,9,10,11,12],
     
     "1+_control_signals" : ["/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K1+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K2+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K3+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K4+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K5+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K6+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K7+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K8+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K9+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K10+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K11+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/1+_K12+.csv"],
     
     "parsed_data_files": ["/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K1+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K2+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K3+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K4+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K5+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K6+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K7+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K8+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K9+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K10+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K11+.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/data/2020-11-10_vary_pulse_duration/5ms/parsed_data/K12+.csv"],
     
     "h" : 1000e-6,
     
     "output_directory" :  "/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/",
     
     "output_file_name" : "abc_pulse=5ms_h=1e-3.csv"
     
     }






'''
Parameters for running the optimisation routine 'opt_neTe.py'
'''

p = {
     "abc_file_path" 		: "./solution_abc_2020-06-29_w_noise_h=10e-6.csv", 	# Path to the file containing the abc parameters
     "voronov_file_path" 	: "./voronov_k.csv", 					# Relative path to the voronov coefficient file
     "cStates" 		: [5], 						# List of charge states on which to perform optimisation
     "ne_lo" 			: 11, 							# log10 of the lower limit of electron density (1/cm3)
     "ne_hi" 			: 12.41664051, 					# log10 of the upper limit of electron density (1/cm3)
     "Te_lo" 			: 10, 							# Lower limit of Te (eV)
     "Te_hi" 			: 10e3, 						# Upper limit of Te (eV)
     "MC_unc_lo" 		: -.6, 						# Lower limit of Voronov uncertainty 
     "MC_unc_hi" 		: .6, 							# Upper limit of Voronov uncertainty
     "N" 			: 100, 						# Number of elements in the ne-vector
     "number_of_MC_iters" 	: 10 							# Number of Monte Carlo iterations to perform
     }






'''
Parameters for plotting.
'''

# RK4 fits
# charge_states = list of charge states for which to perform fitting
# h = RK4 time step to use in fitting
plotting_fits = {
    "charge_states" : [4,5,6,7,8,9,10,11],
    "h" : 10e-6,
    "output_file_name" : "fit_h=10e-6", # Will automatically append "_q=X+"
    "output_file_type" : ".eps"
    }






# solution_set_files = list of files containing the solution sets obtained from opt_neTe.py
# charge_states = list of charge states corresponding to above files
# Te_lo = lower limit for Te to accept in plots (eV)
# Te_hi = upper -"-
# ne_lo = lower limit for ne to accept in plots (1/cm3)
# ne_hi = upper -"-
# F = upper limit of the penalty function value.
plotting_results = {
    
    
"solution_set_files" :
[
"/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/unc_lo=-60%_unc_hi=60%_MC_iters=100_N=1000_q=4.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/unc_lo=-60%_unc_hi=60%_MC_iters=100_N=1000_q=5.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/unc_lo=-60%_unc_hi=60%_MC_iters=100_N=1000_q=6.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/unc_lo=-60%_unc_hi=60%_MC_iters=100_N=1000_q=7.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/unc_lo=-60%_unc_hi=60%_MC_iters=100_N=1000_q=8.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/unc_lo=-60%_unc_hi=60%_MC_iters=100_N=1000_q=9.csv",
"/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/unc_lo=-60%_unc_hi=60%_MC_iters=100_N=1000_q=10.csv"
],
    
"charge_states" : [4,5,6,7,8,9,10],


#
# General plotting parameters
#
"Te_lo" : 0, # To accept all data points, set lower limit = 0, and upper = np.inf
"Te_hi" : np.inf,
"ne_lo" : 0,
"ne_hi" : np.inf,
"F" : 1e-4,
"conf" : 0.341,

#
# Plot number of solutions against F
#
"plotting_vs_F" : 
{
"plot_num_of_solutions_vs_F" : False, # Set True / False
"list_of_Fs" : [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6],
"y_lo" : 0,
"y_hi" : None,
},

#
# For plotting characteristic times vs F
#
"plotting_time_vs_F" : 
{
 "plot_or_not" : True, # Set True / False
 "list_of_Fs"  : [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6],
 "plot_tau" : False,
 "plot_inz" : False,
 "plot_cx"  : False,
 "tau_marker" : ".",
 "tau_color"  : "r",
 "inz_marker" : "s",
 "inz_color"  : "k",
 "cx_marker"  : "^",
 "cx_color"   : "b"
},

#
# For plotting characteristic times vs charge state
#
"plotting_time_vs_q" : 
{
 "plot_or_not" : True, # Set True / False
 "plot_tau" : False,
 "plot_inz" : False,
 "plot_cx"  : False,
 "tau_yscale" : "linear",
 "tau_marker" : ".",
 "tau_color"  : "r",
 "inz_yscale" : "linear",
 "inz_marker" : "s",
 "inz_color"  : "k",
 "cx_yscale"  : "linear",
 "cx_marker"  : "^",
 "cx_color"   : "b"
},


#
# For plotting energy content vs charge state
#
"plotting_Ee_vs_q" :
{
 "plot_or_not" : True, # Set True / False
 "marker" : "s",
 "color" : "k",
 "y_hi" : 1e17,
 "y_lo" : 1e14
},


#
# For plotting triple products vs charge state
#
"plotting_triple_vs_q" :
{
 "plot_or_not" : True, # Set True / False
 "marker" : "s",
 "color" : "m",
 "y_hi" : None,
 "y_lo" : None
},

#
# Select whether or not you wish to output data as .csv
#
"output_tau_vs_q" : False,
"output_inz_vs_q" : False,
"output_cx_vs_q" : False,
"output_Ee_vs_q" : False,
"output_triple_vs_q" : True,



"output_directory" : "/home/miha/Work/research/ppp-2/pulse_duration/results/5ms/"

}



