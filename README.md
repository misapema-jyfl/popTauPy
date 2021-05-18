# popTauPy

A numerical code for determining plasma characteristic values from 
short pulse mode 1+ injection induced extraction current transients
in a CB-ECRIS.

## Major changes since v1.1.0
* A GUI has been developed to make it easier to create 
the parameters.yaml files.

* The code is now run from the terminal by calling the 
masterScript.py with the parameters (.yaml) file 
given as its argument. E.g:
> python3 masterScript.py params.yaml

* Changed parameters input from parameters.py to a .yaml file.
As an example, see the file 'params.yaml'.

* Added possibility to use interpolation functions 
to evaluate the rate coefficients. Currently,
only the Maxwell-Boltzmann distribution is available.
This interpolation function can be chosen by specifying 
"method" = "interpMB" in the parameters.py file.

* Solutions from 'opt_neTe.py' are now output in terms of 
average electron energy. This makes the plotting routines of 
version 1.2.0 incompatible with those of previous versions.

* Merged solution set plotting with 'plot_results.py'
and consequently removed 'plot_solution_sets.py'

* .png plots (dpi=300) will now also be plotted automatically 
when plotting results (along with the .eps figures).

* Included fit plotting in 'opt_abc.py'. Plots will be saved in the 
same directory with the abc solution file.

* The uncertainties for the rate coefficient function are now 
to be determined in the corresponding elemental data file.
E.g. in the case where the rate coefficient is evaluated using the 
Voronov formula, one would specify the charge state specific uncertainties
in the file "./elemental_data/voronov_k.csv"
