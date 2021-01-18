# popTauPy

A numerical code for determining plasma characteristic values from 
short pulse mode 1+ injection induced extraction current transients
in a CB-ECRIS.

## Instructions for use

* These instructions are incomplete will be updated at a later date.

### a, b, c fitting

* The sample data (see link in manuscript) has already been parsed to be acceptable for the script
'opt_abc.py'. 

* The file solution_abc_2020-06-29_w_noise_h=10e-6.csv contains a,b,c parameters, their uncertainties and the reduced chi**2 values corresponding to the fits made to the sample data. 

### Postidctions from a,b,c parameters

The script in 'opt_neTe.py' runs the optimisation routine. You may change 
the code runtime parameters (e.g. limits of electron density, and temperature),
by modifying the dictionary within 'parameters.py'. 
