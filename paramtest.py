#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 11:05:47 2021

@author: miha
"""

import os

class ParamTest:
    def __init__(self, parameters):
        self.params = parameters
        self.availableMethods = ["voronov", "interpMB"]
        
        
    def do_tests(self):
        print("\nRunning tests...")
        errors = 0
        warnings = 0
        errors += self.checkDirs()
        e, w = self.checkOptimizerValues()
        errors += e
        warnings += w
        
        print("Process excited with {} errors and {} warnings."\
              .format(errors,warnings))
        print("Commencing...\n")
        
        return errors, warnings
        
    def checkDirs(self):
        '''
        Check that each directory path is given and that they point 
        to an actual directory.
        '''
        
        errors = 0
        # Always check the absolutely necessary directories
        if not os.path.isdir(self.params["general"]["working_directory"]):
            print("Error: Invalid working directory!")
            errors+=1
        if not os.path.isdir(self.params["general"]["elemental_data_directory"]):
            print("Error: Invalid elemental data directory!")
            errors+=1
        if not os.path.isdir(self.params["general"]["save_to_path"]):
            print("Error: Invalid save to path!")
            errors+=1
        
        # Only check parsing paths if do_parse == True
        if self.params["general"]["do_parse"]:
            if not os.path.isdir(self.params["parsing"]["raw_data_path"]):
                print("Error: Invalid path to raw data directory!")
                errors+=1
            if not os.path.isdir(self.params["parsing"]["parsed_data_path"]):
                print("Error: Invalid path to parsed data directory!")
                errors+=1
                
        # Only check data path if do_abc == True
        if self.params["general"]["do_abc"]:
            if not os.path.isdir(self.params["data"]["parsed_data_path"]):
                print("Error: Invalid path to parsed data directory!")
                errors+=1
                
        # # Only check raw data path if do_parse == True
        # if self.params["general"]["do_parse"]:
        #     if not os.path.isdir(self.params["parsing"]["raw_data_path"]):
        #         errors+=1
        #     if not os.path.isdir(self.params["parsing"]["parsed_data_path"]):
        #         errors+=1
        
        return errors
    
    def checkOptimizerValues(self):
        '''
        Check that given optimizer values are OK.
        '''
        
        warnings = 0
        errors = 0
        
        # Only check the parameters if do_neTe == True
        if self.params["general"]["do_neTe"]:
            if not self.params["optimizer"]["rk_time_step"] <= 1.0e-3:
                print("Warning: Runge-Kutta time step greater than 1ms.")
                warnings+=1
            if not self.params["optimizer"]["rate_coefficient_method"] in self.availableMethods:
                print("Error: Chosen method not available!")
                print("Available methods: ", self.availableMethods)
                errors+=1
            if self.params["optimizer"]["Ee_hi"] > 10.0e+3:
                print("Warning: Cross section data may be unreliable beyond 10 keV!")
                warnings += 1
                
        return errors, warnings