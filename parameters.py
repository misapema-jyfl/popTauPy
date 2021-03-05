# -*- coding: utf-8 -*-
'''

Settings file for running the optimisation code popTauPy.













IMPORTANT INFORMATION BELOW!






N.B. When designating an output directory, include the '/' at the end 
of the address! E.g. "./results/" (NOT "./results") 






IMPORTANT INFORMATION ABOVE!

'''












import numpy as np

'''
Parameters for running parse_data.py

charge_states : Charge states for which transient data is available.
1+_control_signals : The 1+ control signals corresponding to n+ transients.
n+_signals : The extracted n+ transients.
save_to_path : Directory under which to save the parsed data files.
elemental_symbol : E.g. 'K' for potassium. Used for naming convention.
'''

parsing_parameters = {
    
"charge_states":[1,2,3,4,5,6,7,8,9,10,11,12],

"1+_control_signals" : 
[
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K1+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K2+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K3+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K4+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K5+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K6+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K7+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K8+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K9+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K10+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K11+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/1+_K12+.csv"
],
    
"n+_signals" : 
[
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K1+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K2+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K3+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K4+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K5+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K6+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K7+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K8+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K9+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K10+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K11+.csv",
"/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/footer_removed_data/K12+.csv"
],
    
"save_to_path" : "/home/miha/Work/research/ppp-2/experimental_data/gas_dosing/data/141e-9/parsed_data/",

"elemental_symbol" : "K"
}


    
    
    
    
    
    
    
    
    
    
    
'''
Parameters for running the optimisation routine 'opt_abc.py'.

charge_states : list of available charge states (in numerical order)
1+_control_signals : list of data files corresponding to the 1+ control signals
parsed_data_files : list of the parsed data files (t=0 set the same, b.g. reduced, in numerical order)
h : Runge-Kutta algorithm time step (units of s)
output_directory : Directory under which to save the obtained abc parameters.
output_file_name : File name under which to save the obtained abc parameters. 

'''

# You may use the clever generator below to generate the list of 
# 1+ control signals and n+ extraction signals (parsed!).
# Simply specify the directory where the data files are located,
# in the variable 'fileLocDir', and the charge states in list 'cStates',
# and place the variable containing the list into the 
# dictionary. The loop assumes that the files are named using 
# the logic "1+_KX+.csv" for the 1+ control signal corresponding 
# to transient of Kx+, and "KX+.csv" for the transient of Kx+.
# You may also specify the list of files manually, if you wish.
fileLocDir = "/n/work01/misapema/data/gas/81e-9/"
cStates = [1,2,3,4,5,6,7,8,9,10,11,12]

# Make the list of control signals
ctrl_files = []
for q in cStates:
    s = ("1+_K",str(q),"+.csv")
    s = "".join(s)
    f = fileLocDir + s
    ctrl_files.append(f)

# Make the list of N+ signals
n_files = []
for q in cStates:
    s = ("K",str(q),"+.csv")
    s = "".join(s)
    f = fileLocDir + s
    n_files.append(f)
    
d = {
     
"charge_states" : cStates,

"1+_control_signals" : ctrl_files,
     
"parsed_data_files": n_files,
     
"h" : 10e-6,

"output_directory" :  "/n/work01/misapema/results/gas/81e-9/",

"output_file_name" : "abc_pulse=81e-9_h10e-6s.csv"
     
}




    
    
    
    
    
    
    
    


'''
Parameters for running the optimisation routine 'opt_neTe.py'
'''

# By default, the path to the abc file is the one specified 
# in dictionary 'd' for the output file name. 
# You may change this manually, if you wish to optimize 
# using another set of abc parameters.
abcLoc = d["output_directory"] + d["output_file_name"]

# By default, the available charge states are 
# determined based on the charge states given in 
# dictionary 'd'. You may also wish to specify a different list.
availableStates = d["charge_states"][1:-1]

# By default, the output directory is the one specified 
# in the dictionary 'd', i.e. same as for the abc file.
outDir = d["output_directory"]

p = {
     "abc_file_path" 		: abcLoc, 	# Path to the file containing the abc parameters
     "voronov_file_path" 	: "/n/work01/misapema/popTauPy-v1.1.1/voronov_k.csv", 					# Relative path to the voronov coefficient file
     "cStates" 		: availableStates, 						# List of charge states on which to perform optimisation
     "ne_lo" 			: 11, 							# log10 of the lower limit of electron density (1/cm3)
     "ne_hi" 			: 12.41664051, 					# log10 of the upper limit of electron density (1/cm3)
     "Te_lo" 			: 10, 							# Lower limit of Te (eV)
     "Te_hi" 			: 10e3, 						# Upper limit of Te (eV)
     "MC_unc_lo" 		: -.6, 						# Lower limit of Voronov uncertainty 
     "MC_unc_hi" 		: .6, 							# Upper limit of Voronov uncertainty
     "N" 			: 1000, 						# Number of elements in the ne-vector
     "number_of_MC_iters" 	: 1000, 							# Number of Monte Carlo iterations to perform
     "output_directory" : outDir # Directory to save resultant solution sets.
     }










'''
Parameters for plotting.
'''

# RK4 fits
# charge_states = list of charge states for which to perform fitting
# h = RK4 time step to use in fitting
plotting_fits = {
    "abc_file_path" : "/home/miha/Work/research/ppp-2/analysis/uW_sweep/results/150W/abc_pulse=150W_h10e-6s.csv",
    "charge_states" : [5,6,7,8,9,10,11],
    "h" : 10e-6,
    "output_directory" : "/home/miha/Work/research/ppp-2/analysis/uW_sweep/results/150W/"
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
"/home/miha/Work/research/ppp-2/analysis/uW_sweep/results/250W/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=6.csv",
"/home/miha/Work/research/ppp-2/analysis/uW_sweep/results/250W/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=7.csv",
"/home/miha/Work/research/ppp-2/analysis/uW_sweep/results/250W/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=8.csv",
"/home/miha/Work/research/ppp-2/analysis/uW_sweep/results/250W/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=9.csv",
"/home/miha/Work/research/ppp-2/analysis/uW_sweep/results/250W/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=10.csv"
],
    
"charge_states" : [6,7,8,9,10],


#
# General plotting parameters
#
"Te_lo" : 0, # To accept all data points, set lower limit = 0, and upper = np.inf
"Te_hi" : np.inf,
"ne_lo" : 0,
"ne_hi" : np.inf,
"F" : 1e-7,
"conf" : 0.341,
"bias_max" : 0.6, # Maximum length of a component of the bias vector

#
# Plot number of solutions against F
#
"plotting_vs_F" : 
{
"plot_num_of_solutions_vs_F" : True, # Set True / False
"list_of_Fs" : [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9],
"y_lo" : 1,
"y_hi" : None,
},


#
# Plot number of solutions against bias_max
#
"plotting_vs_bias_max":
{
"plot_num_of_solutions_vs_bias_max" : False, # Set True/False
"list_of_bias_max" : [0.6, 0.5, 0.4, 0.3, 0.2, 0.1],
"y_lo" : 1,
"y_hi" : None,
},


#
# For plotting characteristic times vs F
#
"plotting_time_vs_F" : 
{
 "plot_or_not" : False, # Set True / False
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
# For plotting characteristic times vs bias max
#
"plotting_time_vs_bMax" : 
{
 "plot_or_not" : False, # Set True / False
 "list_of_Bs"  : [0.6, 0.5, 0.4, 0.3, 0.2, 0.1],
 "plot_tau" : True,
 "plot_inz" : True,
 "plot_cx"  : True,
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
 "plot_tau" : True,
 "plot_inz" : True,
 "plot_cx"  : True,
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
 "y_hi" : None,
 "y_lo" : None
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
"output_tau_vs_q" : True,
"output_inz_vs_q" : True,
"output_cx_vs_q" : True,
"output_Ee_vs_q" : True,
"output_triple_vs_q" : True,



"output_directory" : "/home/miha/Work/research/ppp-2/analysis/uW_sweep/results/250W/"

}



