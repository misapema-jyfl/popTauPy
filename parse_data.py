#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 10:42:06 2021

@author: miha
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import sys
font = {"family": "normal",
        "weight": "normal",
        "size": 15}

matplotlib.rc("font",**font)


class Parser:
    
    def __init__(self, params):
        
        self.element = params["parsing"]["species"].upper()
        
        self.fileDir = params["parsing"]["raw_data_path"]
        self.outDir = params["parsing"]["parsed_data_path"]
        
        self.header = params["parsing"]["header"]
        if self.header == 0:
            self.header = []
        elif self.header > 0:
            self.header = np.linspace(0,self.header,self.header)
            
        self.footer = params["parsing"]["footer"]
        self.sep = params["parsing"]["separator"]
        
        # Specify the naming conventions. 
        # Use placeholder {} for the (variable) name specifier.
        # self.nPlusInConvention = params["parsing"]["n_plus_naming_convention"]
        # self.onePlusInConvention = params["parsing"]["one_plus_naming_convention"]
        # self.names = params["parsing"]["names"]
        
        # Retrieve 1+ & n+ filenames
        self.onePlusFiles = params["parsing"]["one_plus_file_names"]
        self.nPlusFiles = params["parsing"]["n_plus_file_names"]        
        self.names = [(one, n) for one, n in zip(self.onePlusFiles, self.nPlusFiles)]
        
        self.cStates = params["parsing"]["available_charge_states"]
        
        # Specify the output file naming conventions for 
        # the one plus data and the n+ data.
        # Use the placeholder {} in the convention.
        # The charge state will be substituted for the placeholder.
        self.onePlusOutConvention = "1+_{}{}+.csv"
        self.nPlusOutConvention = "{}{}+.csv"
    
        
    def nameFile(self, convention, name):
        '''
        Creates a filename based on a naming convention and 
        a variable name.
        
        Give the naming convention as a string with a placeholder.
        E.g. convention = "2021 02 24 - {} courant N+.txt"
        
        Give the name as a string.
        E.g. name = "16 59 53"
        
        The function will then place the name into the placeholder 
        within the convention.
        '''
        filename = convention.format(name)
        
        return filename
    
    
    def loadData(self, fileName):
        '''
        Read and return the t, i data from the given file as a DataFrame.
        '''
        
        filePath = self.fileDir + fileName
        
        try:
            df = pd.read_csv(filePath, sep=self.sep, engine="python",
                             skipfooter=self.footer, header=None, 
                             skiprows=self.header,
                             names=["t","i"])
        except:
            print("Error: file {} not found!".format(fileName))
            print("Is this the correct path to the file?\n")
            print(filePath)
            print("\nHope that helps! Exiting...")
            sys.exit()
            
        return df
    

    def parseData(self, cState, name):
        
        # Get the filenames for the one plus and n+ transient signals 
        # onePlusFileName = self.nameFile(self.onePlusInConvention, name)
        # nPlusFileName = self.nameFile(self.nPlusInConvention, name)
        onePlusFileName = name[0]
        nPlusFileName = name[1]
        
        # Open the data files
        dfOne = self.loadData(onePlusFileName)
        dfN = self.loadData(nPlusFileName)
        
        # Get the time series data
        t_1, i_1 = dfOne["t"], dfOne["i"]
        t_n, i_n = dfN["t"], dfN["i"]
        
        # Remove time offset, i.e. set 1+ rise onset as t=0
        t_off = t_1[ i_1 > 0.1*max(i_1) ].values[0]
        t_1 = t_1 - t_off
        t_n = t_n - t_off
        
        # Remove background from n+ signal
        bg = np.average(i_n[t_n<0])
        i_n = i_n-bg
        
        # Remove background from 1+ signal
        bg = np.average(i_1[t_1<0])
        i_1 = i_1-bg
        
        # Normalize 1+ command signal
        i_1 = max(i_n)*i_1/max(i_1)

        # Pack parsed data to DataFrame and save to .csv
        df_1 = pd.DataFrame([t_1,i_1]).transpose()
        df_n = pd.DataFrame([t_n,i_n]).transpose()
        
        s = (self.outDir, self.onePlusOutConvention.format(self.element, cState))
        onePlusOutName = "".join(s)
        
        s = (self.outDir, self.nPlusOutConvention.format(self.element, cState))
        nPlusOutName = "".join(s)
        
        df_1.to_csv(onePlusOutName, index=None)
        df_n.to_csv(nPlusOutName, index=None)
        
        # Return the data frames for possible plotting
        return df_1, df_n
    
    
    def do_parse(self):
        
        # Check that figures have a directory to save to
        s = (self.outDir, "figures/")
        figDir = "".join(s)
        if not os.path.isdir(figDir):
            os.mkdir(figDir)
                
        for cState, name in zip(self.cStates, self.names):
        
            # Make plots
            fig, ax = plt.subplots()
            
            # Do parsing            
            df_1, df_n = self.parseData(cState, name)
            
            # Plot for sanity check
            ax.plot(df_1["t"], df_1["i"])
            ax.plot(df_n["t"], df_n["i"])
            
            # Determine save to path and save
            s = (figDir, self.element, str(cState), ".png")
            savePath = "".join(s)
            fig.savefig(savePath)
        