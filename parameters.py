# -*- coding: utf-8 -*-


'''
Parameters for running the optimisation routine 'opt_abc.py'
'''
# charge_states = list of available charge states (in numerical order)
# 1+ control signals = list of data files corresponding to the 1+ control signals
# parsed_data_files = list of the parsed data files (t=0 set the same, b.g. reduced, in numerical order)
# h = Runge-Kutta algorithm time step (units of s)
# output_file_name = file name under which to save the obtained abc parameters  
d = {
     
     "charge_states" : [1,2,3,4,5,6,7,8,9,10,11,12],
     
     "1+_control_signals" : ["/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K1+.csv",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K2+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K3+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K4+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K5+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K6+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K7+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K8+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K9+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K10+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K11+.txt",
"/home/miha/Work/research/ppp/data/2020 06 29 Measurements/1+_ctrl/1+_K12+.txt"],
     
     "parsed_data_files": ["/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K1+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K2+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K3+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K4+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K5+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K6+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K7+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K8+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K9+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K10+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K11+.csv",
"/home/miha/Work/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K12+.csv"],
     
     "h" : 1000e-6,
     
     "output_file_name" : "abc_2020-06-29_w_noise_h=1e-3.csv"
     
     }






'''
Parameters for running the optimisation routine 'opt_neTe.py'
'''

p = {
     "abc_file_path" 		: "./solution_abc_2020-06-29_w_noise_h=10e-6.csv", 	# Relative path to the file containing the abc parameters
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







