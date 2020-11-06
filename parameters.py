# -*- coding: utf-8 -*-

'''
Parameters for running the optimisation routine 'opt_neTe.py'
'''

p = {
     "abc_file_path" : "./solution_abc_2020-06-29_w_noise_h=10e-6.csv", # Relative path to the file containing the abc parameters
     "voronov_file_path" : "./voronov_k.csv", # Relative path to the voronov coefficient file
     "cStates" : [5], # List of charge states on which to perform optimisation
     "ne_lo" : 11, # log10 of the lower limit of electron density (1/cm3)
     "ne_hi" : 12.41664051, # log10 of the upper limit of electron density (1/cm3)
     "Te_lo" : 10, # Lower limit of Te (eV)
     "Te_hi" : 10e3, # Upper limit of Te (eV)
     "MC_unc_lo" : -.6, # Lower limit of Voronov uncertainty 
     "MC_unc_hi" : .6, # Upper limit of Voronov uncertainty
     "N" : 100, # Number of elements in the ne-vector
     "number_of_MC_iters" : 10 # Number of Monte Carlo iterations to perform
     }