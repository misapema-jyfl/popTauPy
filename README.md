# popTauPy

A numerical code for determining plasma characteristic values from 
short pulse mode 1+ injection induced extraction current transients
in a CB-ECRIS.

## Instructions for use

* These instructions will be updated at a later date.

### a, b, c fitting

* The sample data has already been parsed to be acceptable for the script
'opt_abc.py'

### Postidctions from a,b,c parameters

The script in 'opt_neTe.py' runs the optimisation routine. You may change 
the code runtime parameters (e.g. limits of electron density, and temperature),
by modifying the dictionary within 'parameters.py'. 

* The abc_file_path must contain a path to the file, which contains the 
abc parameters.

* The voronov_file_path must contain a path to the file, which contains
the fitting parameters for the Voronov semi-empirical formula for the
rate coefficient.

* The parameter 'cStates' must be a list of integers, which correspond to the
charge states of the ions, for which data is available in the abc file.

* The electron density limits (ne_lo & ne_hi) must be given as the base-10 
logarithm of the electron number density (1/cm^3). 

* The electron temperature (Te_lo & Te_hi) must be given in units of eV.

* The parameter N defines the number of elements in the ne-vector. The higher
this is, the denser will be the obtained solution set of (ne, Te)-pairs.

* The parameter number_of_MC_iterations defines the number of times the 
optimisation routine will be run. In each iteration a different set of 
random biases will be applied on the rate coefficients provided by the 
Voronov formula (within uncertainty limits, see below). The higher this value
is, the more combinations of biases will be tested.

* The Voronov formula uncertainty lower and upper limits are given in the 
parameters 'unc_MC_lo' and 'unc_MC_hi', respectively. Voronov reports in 
his paper () that the uncertainty of the literature values used in his 
fitted formula, is 60 %. This corresponds to unc_MC_lo=-.6, and unc_MC_hi=.6.

The resultant data file (.csv) will be output to the directory "./results/".
The the results folder will be created automatically by the code, if
it does not exist already.
