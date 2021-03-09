# popTauPy

A numerical code for determining plasma characteristic values from 
short pulse mode 1+ injection induced extraction current transients
in a CB-ECRIS.

## Major changes since v1.1.1

* Solutions from 'opt_neTe.py' are now output in terms of 
average electron energy. This makes the plotting routines of 
version 1.2.0 incompatible with those of previous versions.

* Merged solution set plotting with 'plot_results.py'
and consequently removed 'plot_solution_sets.py'

* .png plots (dpi=300) will now also be plotted automatically 
when plotting results (along with the .eps figures).

* Included fit plotting in 'opt_abc.py'. Plots will be saved in the 
same directory with the abc solution file.

* The lower and upper limits of electron density are now specified 
simply as (1/cm3) and no longer as the corresponding base-10 logs.

* Removed option to plot as a function of the maximum component 
of the bias vector.

* The uncertainties for the rate coefficient function are now 
to be determined in the corresponding elemental data file.
E.g. in the case where the rate coefficient is evaluated using the 
Voronov formula, one would specify the charge state specific uncertainties
in the file "./elemental_data/voronov_k.csv"

## Minor changes

### Updated parameters.py
* Simplified specification of data file locations.
* Added default locations for some parameters.
* Added key "elemental_data_dir" to dictionary 'p'
* Added key "method" to dictionary 'p'
* Added key "species" to dictionary 'p'
* Removed keys "MC_unc_lo" and "MC_unc_hi" from dictionary 'p'


### Updated opt_neTe.py
* Commented out import of numdifftools as unnecessary.
* Added 'if __name__ == "__main__":' condition
in the context of multiprocessing to fix compatibility issue with Windows
operating systems.
* voronov_rate is now called as a function of the average energy of
the MB distribution, instead of the temperature.
