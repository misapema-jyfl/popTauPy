#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 12:56:36 2021

A Graphical User Interface for generating parameter .yaml files
to run the Consecutive Transients Analyzer.

@author: miha
"""

import tkinter as tk

# Insantiate root window
root = tk.Tk()

## Configure the root window
root.title("Parameters generator GUI")
root.geometry("500x500")
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)


# Execute mainloop
root.mainloop()