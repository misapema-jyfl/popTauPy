# Developer's note

This repository is no longer maintained. The new version of the code
can be found from the [ct-analyzer](https://www.github.com/misapema-jyfl/ct-analyzer) repository.

# popTauPy

A numerical code for determining plasma characteristic values from 
short pulse mode 1+ injection induced extraction current transients
in a CB-ECRIS.

## Compatibility

This code has been developed and tested on Linux Ubuntu (20.04). 
It is known that the code is currently incompatible with Windows OS.
Sorry for the inconvenience. Use Linux.

## Instructions for use

* Use the GUI to generate a parameters file.

* Execute the masterScript.py from terminal, passing the parameters file
as its argument. E.g: In the working directory, run
> python3 masterScript.py params.yaml

The code expects to find the 1+ injection pulse 
corresponding to the measured n+ current as separate files, 
and it expects the time series data to be stored in two columns,
with time in the first column, and the signal in the second column.

A conversion factor must be passed as a parameter to the code
to convert the signal to units of Amperes. 

## Disclaimer

This code is by no means foolproof, and there are many ways to
use it in a way that produces undesirable behavior. 
There is also a lack of flexibility and features.
Then again, there are currently only two people in the world
who are actively using the code, and the other guy will 
(hopefully) forgive me.

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
