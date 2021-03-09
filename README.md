# popTauPy

A numerical code for determining plasma characteristic values from 
short pulse mode 1+ injection induced extraction current transients
in a CB-ECRIS.

## Changes since v1.1.1

* Updated parameters.py
** Simplified specification of data file locations.
** Added default locations for some parameters.

* Updated opt_neTe.py
** Commented out import of numdifftools as unnecessary.
** Added 'if __name__ == "__main__":' condition
in the context of multiprocessing to fix compatibility issue with Windows
operating systems.