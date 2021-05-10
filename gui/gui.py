#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 12:56:36 2021

A Graphical User Interface for generating parameter .yaml files
to run the Consecutive Transients Analyzer.

@author: miha
"""

import tkinter as tk


class Excel(tk.Frame):
    def __init__(self, master, rows, columns, width):
        super().__init__(master)

        for i in range(columns):
            self.make_entry(0, i+1, width, f'C{i}', False) 

        for row in range(rows):
            self.make_entry(row+1, 0, 5, f'R{row}', False)
                
            for column in range(columns):
                self.make_entry(row+1, column+1, width, '', True)

    def make_entry(self, row, column, width, text, state):
        e = tk.Entry(self, width=width)
        if text: e.insert(0, text)
        e['state'] = tk.NORMAL if state else tk.DISABLED
        e.coords = (row-1, column-1)
        e.grid(row=row, column=column)

class Window:
    
    def __init__(self):
        
        self.root = tk.Tk()
        
        self.root.title("Parameters YAML file generator")
        self.root.minsize(500,220)
        self.padx, self.pady = 5, 5
        self.createLayout()
        
    def createLayout(self):
        '''Creates the layout'''
        
        # Specify working directory
        self.createLabel(labelText="Working directory:",
                         row=0,
                         column=0,
                         weight=(0,0),
                         sticky="n")
        self.workingDirInp = self.createTextBox(row=0,
                                                column=1,
                                                columnspan=1,
                                                weight=(0,1),
                                                sticky="nwe")
        
        # Specify save to path
        self.createLabel(labelText="Save to path:",
                         row=1,
                         column=0,
                         weight=(0,0),
                         sticky="nwe")
        self.saveToPathInp = self.createTextBox(row=1,
                                                column=1,
                                                columnspan=1,
                                                weight=(0,1),
                                                sticky="nwe")
    
        # Create check boxes for routines to run
        self.createLabel(labelText="Routines to run:",
                         row=2,
                         column=0,
                         sticky="wn",
                         weight=(0,0))
        self.parseVar = self.createCheckBox(boxText=" Parse data",
                            row=3,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            command=lambda: self.activateCheck(self.parseVar.get(),
                                                               self.parseButton))
        self.fitVar = self.createCheckBox(boxText=" abc fitting",
                            row=4,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            command=lambda: self.activateCheck(self.fitVar.get(),
                                                               self.fitButton))
        self.optVar = self.createCheckBox(boxText=" (ne, Ee) optimization",
                            row=5,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            command=lambda: self.activateCheck(self.optVar.get(),
                                                               self.optButton))
        self.plotVar = self.createCheckBox(boxText=" Plot data",
                            row=6,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            command=lambda: self.activateCheck(self.plotVar.get(),
                                                               self.plotButton))
        # Create open dialog buttons
        self.parseButton = self.createButton(buttonText="Parsing",
                          row=3,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=self.createParsingWindow,
                          state=tk.DISABLED)
        self.fitButton = self.createButton(buttonText="abc fitting",
                          row=4,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=self.dummy,
                          state=tk.DISABLED)
        self.optButton = self.createButton(buttonText="(ne, Ee)-optimization",
                          row=5,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=self.dummy,
                          state=tk.DISABLED)
        self.plotButton = self.createButton(buttonText="Plotting",
                          row=6,
                          column=1,
                          weight=(0,1),
                          sticky="news",
                          command=self.dummy,
                          state=tk.DISABLED)
        
    
        

    def createParsingWindow(self):
        '''Create the parsing dialog window.'''
        parsingWindow = tk.Toplevel(self.root)
        parsingWindow.title("Parsing parameters")
        parsingWindow.minsize(500, 220)
        self.parsingWindow = parsingWindow
        # Ask for the 1+ naming convention
        self.createLabel("1+ naming convention:",
                         row=0,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        self.onePlusConvention, self.onePlusTextBox = self.createTextBox(row=0,
                           column=1,
                           sticky="nwe",
                           weight=(0,1),
                           window=parsingWindow)
        # Ask for the N+ naming convention
        self.createLabel("N+ naming convention:",
                         row=1,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        self.nPlusConvention, self.nPlusTextBox = self.createTextBox(row=1,
                           column=1,
                           sticky="nwe",
                           weight=(0,1),
                           window=parsingWindow)
        # Ask for the filename variables
        self.createLabel("Filename variables:",
                         row=2,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        self.filenameVars, self.filenameTextBox = self.createTextBox(row=2,
                           column=1,
                           sticky="nwe",
                           weight=(0,1),
                           window=parsingWindow)
        self.createButton(buttonText="Accept",
                          row=2,
                          column=2,
                          sticky="nwe",
                          weight=(0,0),
                          window=parsingWindow,
                          command=lambda: self.parsingAcceptFilenames(
                                              self.filenameVars.get()
                                              ),
                          padx=0,
                          pady=0)
        self.createButton(buttonText="Clear",
                          row=1,
                          column=2,
                          sticky="nwe",
                          weight=(0,0),
                          window=parsingWindow,
                          command=lambda: self.parsingClearFilenames(),
                          padx=0,
                          pady=0)
        
    def parsingPrintFilenames(self):
        '''Prints out the specified filenames for user to check.'''
        
        for i, name in enumerate(self.onePlusFilenames):
            
            self.createLabel(labelText=name,
                             row = i+3,
                             column = 0,
                             sticky = "wn",
                             weight = (0,0),
                             window = self.parsingWindow)
        
    def parsingClearFilenames(self):
        '''Clear input'''
        
        self.onePlusFilenames = []
        self.nPlusFilenames = []
        self.nPlusTextBox.delete(0,"end")
        self.onePlusTextBox.delete(0,"end")
        self.filenameTextBox.delete(0,"end")
        self.parsingPrintFilenames()
    
    def parsingAcceptFilenames(self, textVar):
        '''Broadcasts the specified filenames to lists,
        and stores them as the class variables:
            self.onePlusFilenames and self.nPlusFilenames.'''
        
        onePlusFilenames = []
        nPlusFilenames = []
        
        splittedTextVar = textVar.split("\n")
        
        c = True
        i = 0
        while c:
            var = splittedTextVar[i]
            
            onePlusFilename = self.onePlusConvention.get().format(var)
            nPlusFilename = self.nPlusConvention.get().format(var)
            onePlusFilenames.append(onePlusFilename)
            nPlusFilenames.append(nPlusFilename)
            
            i+=1
            if splittedTextVar[i] == "":
                c = False
            # If this starts to look like an infinite loop, terminate.
            elif i >= 100: # When this hard-coded value becomes obsolete,
                           # I'll congratulate the Nobel-laureate.
                c = False
        self.onePlusFilenames = onePlusFilenames
        self.nPlusFilenames = nPlusFilenames
        
        self.parsingPrintFilenames()
        
        
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
        None.

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
                      window=None):
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
            
        Returns
        -------
        inputVar : Stringvar
            Text input variable.

        '''
        
        if window==None:
            window=self.root
            
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
        
        inputVar.trace_add("write", self.dummy)
        
        return inputVar, textbox

    def dummy(self, *_):
        '''Dummy function for tests'''
        return 0

# Instantiate the window and run its mainloop
window = Window()
window.root.mainloop()