#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 08:26:57 2021

Utilities class for the GUI.

@author: miha
"""

import tkinter as tk

class Utils:
    def __init__(self):
        self.padx, self.pady = 5, 5
        return None
    
    def activateCheck(self, variable, button):
        '''Check the state of the variable,
        and set the state of the button.'''
        if variable:
            button.config(state=tk.NORMAL)
        else:
            button.config(state=tk.DISABLED)
        
    def createCheckBox(self, boxText,
                       row, 
                       column, 
                       weight=(1,1),
                       sticky="",
                       command=None,
                       window=None):
        '''
        Creates a check box.

        Parameters
        ----------
        boxText : str
            Label text for the check box.
        row : int
            Row to place in.
        column : int
            Column to place in.
        sticky : str, optional
            Sticky. E.g. "we" The default is "".
        weight : array, optional
            Relative weight of the grid cell. Default is (1,1).
        command : callable, optional
            Function to call upon check. Default is None.
        window : window object
            Window into which to broadcast the check box.
            
        Returns
        -------
        variable : IntVar()
            Variable containing the value of the check box. 1 = checked, 
            0 = not checked.

        '''
        
        if window==None:
            window=self.root
        
        window.columnconfigure(index=column, weight=weight[1])
        window.rowconfigure(index=row, weight=weight[0])
        
        variable = tk.IntVar()
        
        cbox = tk.Checkbutton(window, 
                              text=boxText, 
                              variable=variable,
                              command=command)
        

            
        cbox.grid(row=row, 
                  column=column,
                  sticky=sticky,
                  ipadx=self.padx,
                  ipady=self.pady)
        
        variable.trace_add("write", self.dummy)
        
        return variable
        
    
    def createLabel(self, labelText,
                    row,
                    column,
                    sticky="",
                    weight=(1,1),
                    window=None):
        '''
        Creates a label

        Parameters
        ----------
        labelText : str
            Text to show.
        row : int
            Row to place in.
        column : int
            Column to place in.
        sticky : str, optional
            Sticky. E.g. "we" The default is "".
        weight : array, optional
            Relative weight of the grid cell. Default is (1,1).
        window : window object
            Window into which to broadcast the Label.
            
        Returns
        -------
        tk.Label object.

        '''
        
        if window==None:
            window=self.root
        
        window.columnconfigure(index=column, weight=weight[1])
        window.rowconfigure(index=row, weight=weight[0])
        
        label = tk.Label(window, 
                         text=labelText)
        label.grid(column=column,
                   row=row,
                   sticky=sticky,
                   ipadx=self.padx, 
                   ipady=self.pady)
        
        return label
    
    def createButton(self, buttonText,
                     row,
                     column,
                     command,
                     sticky="",
                     weight=(1,1),
                     state=tk.NORMAL,
                     window=None,
                     padx=None,
                     pady=None):
        '''
        Creates a button

        Parameters
        ----------
        buttonText : str
            Text to show on button.
        row : int
            Row to place in.
        column : int
            Column to place in.
        command : callable
            Function to call upon button press.
        sticky : str, optional
            Sticky. E.g. "we" The default is "".
        weight : array, optional
            Relative weight of the grid cell. Default is (1,1).
        state : tk.NORMAL or tk.DISABLED
            State of the Button. Default is tk.NORMAL.
        window : window object
            Window into which to broadcast the Button.

        Returns
        -------
        button : tk.Button object
            The button object.

        '''
        
        if window==None:
            window=self.root
        
        window.columnconfigure(index=column, weight=weight[1])
        window.rowconfigure(index=row, weight=weight[0])
        
        if padx == None:
            padx = self.padx
        if pady == None:
            pady = self.pady
        
        button = tk.Button(window,
                           text=buttonText,
                           command=command,
                           state=state)
        button.grid(column=column, 
                    row=row,
                    sticky=sticky,
                    ipadx=padx, 
                    ipady=pady)
        
        return button

    def createTextBox(self, row,
                      column,
                      columnspan=1,
                      sticky="",
                      weight=(1,1),
                      window=None,
                      command=None):
        '''
        Creates a text box whose input can be stored.

        Parameters
        ----------
        row : int
            Row to grid text box into.
        column : int
            Column to grid text box into.
        columnspan : int
            Number of columns for box to span.
        sticky : string
            Where to stick to in the grid cell.
        weight : array, optional
            Relative weight of the grid cell. Default is (1,1).
        window : window object
            Window into which to broadcast the box.
        command : callable, optional
            Set the callback function for the trace of the inputVar.
            Default callback is Utils.dummy.
            
        Returns
        -------
        inputVar : Stringvar
            Text input variable.
        textBox : tk.Entry object
            TextBox object.

        '''
        
        # if window==None:
        #     window=self.root
        
        if command == None:
            command = self.dummy
            
        window.columnconfigure(index=column, weight=weight[1])
        window.rowconfigure(index=row, weight=weight[0])
        
        inputVar = tk.StringVar()
        textbox = tk.Entry(window,
                           textvariable=inputVar)
        textbox.grid(row=row, 
                     column=column,
                     columnspan=columnspan,
                     sticky=sticky,
                     ipadx=self.padx, 
                     ipady=self.pady)
        
        
            
        inputVar.trace_add("write", command)
        
        return inputVar, textbox

    def dummy(self, *_):
        '''Dummy function for tests'''
        return 0