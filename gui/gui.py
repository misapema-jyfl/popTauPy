#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 12:42:59 2021

@author: miha
"""

import tkinter as tk
import numpy as np 

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
    
class CreateToolTip(object):
    """
    create a tooltip for a given widget
    """
    def __init__(self, widget, text='widget info'):
        self.waittime = 500     #miliseconds
        self.wraplength = 180   #pixels
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.widget.bind("<ButtonPress>", self.leave)
        self.id = None
        self.tw = None

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.showtip)

    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)

    def showtip(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background="#ffffff", relief='solid', borderwidth=1,
                       wraplength = self.wraplength)
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()
            

    
    
class MainWindow:
    
    def __init__(self):
        
        self.root = tk.Tk()
        self.root.title("Parameters YAML file generator")
        self.root.minsize(500,220)
        self.parameters = {
            "workingDir":"",
            "saveToPath":"",
            "species":"",
            "cStates":[],
            "onePlusFilenames":[],
            "nPlusFilenames":[],
            "header":0,
            "footer":0,
            "separator":",",
            "pathToRawData":"",
            "abcFilename":"",
            "rkStepsize":1.0e-3,
            "ne_lo":0,
            "ne_hi":0,
            "Ee_lo":0,
            "Ee_hi":0,
            "N_ne":1000,
            "N_MC":1000,
            "method":""}
        self.workingDir = tk.StringVar()
        self.saveToPath = tk.StringVar()
        self.createLayout()
        
    def createLayout(self):
        '''Creates the MainWindow layout'''
        u = Utils()
        
        # Specify working directory
        workLabel = u.createLabel(labelText="Working directory:",
                            row=0,
                            column=0,
                            weight=(0,0),
                            sticky="n",
                            window=self.root)
        CreateToolTip(workLabel, "Directory containing the code files.\
                      E.g. './misapema/popTauPy-v1.2.3/'\
                          (without the \'\')")
        self.workingDir, wDirBox = u.createTextBox(row=0,
                                                column=1,
                                                columnspan=1,
                                                weight=(0,1),
                                                sticky="nwe",
                                                window=self.root)
        
        # Specify save to path
        saveToLabel = u.createLabel(labelText="Save to path:",
                         row=1,
                         column=0,
                         weight=(0,0),
                         sticky="nwe",
                         window=self.root)
        CreateToolTip(saveToLabel, "All results will be saved to this directory.")
        self.saveToPath, saveToBox = u.createTextBox(row=1,
                                                column=1,
                                                columnspan=1,
                                                weight=(0,1),
                                                sticky="nwe",
                                                window=self.root)
    
        # Create check boxes for routines to run
        u.createLabel(labelText="Routines to run:",
                         row=2,
                         column=0,
                         sticky="wn",
                         weight=(0,0),
                         window=self.root)
        self.parseVar = u.createCheckBox(boxText=" Parse data",
                            row=3,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            command=lambda: u.activateCheck(self.parseVar.get(),
                                                               self.parseButton),
                            window=self.root)
        self.fitVar = u.createCheckBox(boxText=" abc fitting",
                            row=4,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            command=lambda: u.activateCheck(self.fitVar.get(),
                                                               self.fitButton),
                            window=self.root)
        self.optVar = u.createCheckBox(boxText=" (ne, Ee) optimization",
                            row=5,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            command=lambda: u.activateCheck(self.optVar.get(),
                                                               self.optButton),
                            window=self.root)
        self.plotVar = u.createCheckBox(boxText=" Plot data",
                            row=6,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            command=lambda: u.activateCheck(self.plotVar.get(),
                                                               self.plotButton),
                            window=self.root)
        # Create open dialog buttons
        p = ParsingWindow(self)
        self.parseButton = u.createButton(buttonText="Parsing",
                          row=3,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=p.createParsingWindow,
                          state=tk.DISABLED,
                          window=self.root)
        a = abcWindow(self)
        self.fitButton = u.createButton(buttonText="abc fitting",
                          row=4,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=a.createABCWindow,
                          state=tk.DISABLED,
                          window=self.root)
        o = optimizeWindow(self)
        self.optButton = u.createButton(buttonText="(ne, Ee)-optimization", # TODO!
                          row=5,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=o.createOptWindow,
                          state=tk.DISABLED,
                          window=self.root)
        self.plotButton = u.createButton(buttonText="Plotting", # TODO!
                          row=6,
                          column=1,
                          weight=(0,1),
                          sticky="news",
                          command=u.dummy,
                          state=tk.DISABLED,
                          window=self.root)
        self.button = u.createButton(buttonText="Button", #TODO! Placeholder button.
                                     row=1000,
                                     column=1,
                                     weight=(0,1),
                                     sticky="en",
                                     command=self.broadCastParameters,
                                     state=tk.NORMAL,
                                     window=self.root)
        
    def broadCastParameters(self, *_):
        self.parameters["workingDir"] = self.workingDir.get()
        self.parameters["saveToPath"] = self.saveToPath.get()
        print("Parameters:")
        for key in self.parameters:
            print(key, self.parameters[key])
            
class ParsingWindow:
    
    def __init__(self, rootWindow):
        '''Create the parsing dialog window attached to the rootWindow.'''
        self.u = Utils()
        self.dataDir = tk.StringVar()
        self.header = tk.StringVar()
        self.footer = tk.StringVar()
        self.separator = tk.StringVar()
        self.parsingFilenameLabels = {"1+":[],"n+":[]}
        self.species = tk.StringVar()
        self.cState_i = tk.StringVar()
        self.cState_f = tk.StringVar()
        self.rootWindow=rootWindow
        
    
    def createParsingWindow(self):
        '''Generates the parsing window.'''
        u = self.u
        
        parsingWindow = tk.Toplevel(self.rootWindow.root)
        parsingWindow.title("Parsing parameters")
        parsingWindow.minsize(500, 350)
        self.parsingWindow=parsingWindow
        
        # Ask for the data storage directory
        dataDirLabel = u.createLabel(labelText="Path to raw data directory",
                                     row=0,
                                     column=0,
                                     sticky="nw",
                                     weight=(0,0),
                                     window=parsingWindow)
        CreateToolTip(dataDirLabel,"Specify path to directory\
                      containing the raw data files.")
        self.dataDir, self.dataDirTextBox = u.createTextBox(row=0,
                           column=1,
                           sticky="nwe",
                           weight=(0,1),
                           window=parsingWindow,
                           command=self.parsingSaveActive)
        
        # Ask for the 1+ naming convention
        onePlusLabel = u.createLabel("1+ naming convention:",
                         row=1,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        CreateToolTip(onePlusLabel,"Use {} as placeholder.\
                                             File ending is mandatory.\
                                                 E.g. 1+_{}.txt")
        self.onePlusConvention, self.onePlusTextBox = u.createTextBox(row=1,
                           column=1,
                           sticky="nwe",
                           weight=(0,1),
                           window=parsingWindow,
                           command=self.parsingSaveActive)
        # Ask for the N+ naming convention
        nPlusLabel = u.createLabel("N+ naming convention:",
                         row=2,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        CreateToolTip(nPlusLabel,"Use {} as placeholder.\
                                             File ending is mandatory.\
                                                 E.g. n+_{}.txt")
        self.nPlusConvention, self.nPlusTextBox = u.createTextBox(row=2,
                           column=1,
                           sticky="nwe",
                           weight=(0,1),
                           window=parsingWindow,
                           command=self.parsingSaveActive)
        # Ask for the filename variables
        variablesLabel = u.createLabel("Filename variables:",
                         row=3,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        CreateToolTip(variablesLabel, "Copy and paste a column from Excel.\
                      Given values will replace the placeholders ({})\
                          in naming conventions specified above.")
        self.filenameVars, self.filenameTextBox = u.createTextBox(row=3,
                           column=1,
                           sticky="nwe",
                           weight=(0,1),
                           window=parsingWindow)
        
        
        parsingClearButton = u.createButton(buttonText="Clear",
                          row=2,
                          column=2,
                          sticky="nwe",
                          weight=(0,0),
                          window=parsingWindow,
                          command=lambda: self.parsingClearFilenames(),
                          padx=0,
                          pady=0)
        CreateToolTip(parsingClearButton, "Clear stored filenames.")
        
        parsingAcceptButton = u.createButton(buttonText="Accept",
                          row=3,
                          column=2,
                          sticky="nwe",
                          weight=(0,0),
                          window=parsingWindow,
                          command=lambda: self.parsingAcceptFilenames(
                                              self.filenameVars.get()),
                          padx=0,
                          pady=0)
        CreateToolTip(parsingAcceptButton,"Accept and store filenames.")
        
        # Specify data file properties
        # Header
        parsingHeaderLabel=u.createLabel(labelText="Header:", 
                         row=4,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        CreateToolTip(parsingHeaderLabel,\
                      "If file has no header, set 0.\
                          Else, set length of header in rows.")
        self.header, headerBox = u.createTextBox(row=4,
                        column=1,
                        sticky="nw",
                        weight=(0,0),
                        window=parsingWindow,
                        command=self.parsingSaveActive)
        headerBox.config(width=5)
        # Footer
        parsingFooterLabel = u.createLabel(labelText="Footer:", 
                         row=5,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        CreateToolTip(parsingFooterLabel,\
                      "If file has no footer, set 0.\
                          Else, set length of footer in rows.")
        self.footer, footerBox = u.createTextBox(row=5,
                        column=1,
                        sticky="nw",
                        weight=(0,1),
                        window=parsingWindow,
                        command=self.parsingSaveActive)
        footerBox.config(width=5)
        # Data separator
        parsingSepLabel = u.createLabel(labelText="Separator:", 
                         row=6,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        CreateToolTip(parsingSepLabel,"Give data separator.\
                      E.g. '\\t' or ','\
                      without the \'\'")
        self.separator, separatorBox = u.createTextBox(row=6,
                        column=1,
                        sticky="nw",
                        weight=(0,1),
                        window=parsingWindow,
                        command=self.parsingSaveActive)
        separatorBox.config(width=5)
        
        # Ion species and charge states
        speciesLabel = u.createLabel(labelText="Ion species: ", 
                      row=7, 
                      column=0,
                      sticky="nw",
                      weight=(0,0),
                      window=parsingWindow)
        CreateToolTip(speciesLabel, "Name of element, e.g. \'K\'.")
        self.species, speciesBox = u.createTextBox(row=7,
                        column=1,
                        sticky="nw",
                        weight=(0,0),
                        command=self.parsingSaveActive,
                        window=parsingWindow)
        speciesBox.config(width=5)
        
        availableLabel = u.createLabel(labelText="Available charge states",
                      row=8, 
                      column=0,
                      window=parsingWindow,
                      sticky="nw",
                      weight=(0,0))
        CreateToolTip(availableLabel, "Transients for at least 5 charge states are required.\
                      the transients must correspond to consecutive charge states.\
                          Charge states must be in numerically ascending order.\
                              The order must be the same as for the filenames.")
        
        cStateLabel_i = u.createLabel(labelText="i: ", 
                      row=9, 
                      column=0,
                      sticky="ne",
                      weight=(0,0),
                      window=parsingWindow)
        CreateToolTip(cStateLabel_i, "First measured charge state.")
        
        self.cState_i, cStateBox_i = u.createTextBox(row=9,
                        column=1,
                        sticky="nw",
                        weight=(0,0),
                        command=self.parsingSaveActive,
                        window=parsingWindow)
        cStateBox_i.config(width=5)
        
        cStateLabel_f = u.createLabel(labelText="f: ", 
                      row=10, 
                      column=0,
                      sticky="ne",
                      weight=(0,0),
                      window=parsingWindow)
        CreateToolTip(cStateLabel_f, "Last measured charge state.")
        
        self.cState_f, cStateBox_f = u.createTextBox(row=10,
                        column=1,
                        sticky="nw",
                        weight=(0,0),
                        command=self.parsingSaveActive,
                        window=parsingWindow)
        cStateBox_f.config(width=5)
        
        # Create labels for 1+ and n+ filenames
        u.createLabel("1+ filenames:",
                         row=11,
                         column=0,
                         sticky="nwe",
                         weight=(0,0),
                         window=parsingWindow)
        u.createLabel("N+ filenames:",
                         row=11,
                         column=1,
                         sticky="nwe",
                         weight=(0,1),
                         window=parsingWindow)
        
        # Create Save button
        # TODO! Would be nice if this stuck to the lower edge of the window
        # even after the window expands.
        # TODO! When the window expands, the filenames may stretch beyond 
        # the screen. At least theoretically. Not that anyone will ever
        # be able to capture that many consecutive transients...
        self.saveButton = u.createButton(buttonText="Save",
                       row=1000,
                       column=2,
                       command=self.parsingSaveParameters,
                       sticky="ews",
                       weight=(0,0),
                       window=parsingWindow,
                       state=tk.DISABLED)
    
    def parsingSaveActive(self, *_):
        '''If all the fields have been filled (correctly OR incorrectly)
        change the state of the Save button to NORMAL. If the fields
        are cleared, set state back to DISABLED.'''
        lenDir = len(self.dataDir.get())
        lenOnePlus = len(self.parsingFilenameLabels["1+"])
        lenNPlus = len(self.parsingFilenameLabels["n+"])
        lenH = len(self.header.get())
        lenF = len(self.footer.get())
        lenSep = len(self.separator.get())
        lenSpecies = len(self.species.get())
        lenStates = len(self.cState_i.get())*len(self.cState_f.get())
        if (lenDir*lenOnePlus*lenNPlus\
            *lenH*lenF*lenSep*lenSpecies*lenStates > 0):
            self.saveButton.config(state=tk.NORMAL)
        else:
            self.saveButton.config(state=tk.DISABLED)
            
    def parsingPrintFilenames(self):
        '''Prints out the specified filenames for user to check.'''
        u = self.u
        self.parsingFilenameLabels={"1+":[], "n+":[]}
        
        # Print the 1+ filenames
        for i, name in enumerate(self.onePlusFilenames):
            
            label = u.createLabel(labelText=name,
                             row = i+12,
                             column = 0,
                             sticky = "wne",
                             weight = (0,0),
                             window = self.parsingWindow)
            self.parsingFilenameLabels["1+"].append(label)
            
        # Print the N+ filenames
        for i, name in enumerate(self.nPlusFilenames):
            
            label = u.createLabel(labelText=name,
                             row = i+12,
                             column = 1,
                             sticky = "wne",
                             weight = (0,1),
                             window = self.parsingWindow)
            self.parsingFilenameLabels["n+"].append(label)
        
    def parsingClearFilenames(self):
        '''Clear input'''
        
        # Destroy the filename labels
        for label in self.parsingFilenameLabels["1+"]:
            label.destroy()
        for label in self.parsingFilenameLabels["n+"]:
            label.destroy()
            
        self.onePlusFilenames = []
        self.nPlusFilenames = []
        # self.nPlusTextBox.delete(0,"end")
        # self.onePlusTextBox.delete(0,"end")
        # self.filenameTextBox.delete(0,"end")
        self.parsingPrintFilenames()
        self.parsingSaveActive()
    
    def parsingAcceptFilenames(self, textVar):
        '''Broadcasts the specified filenames to lists,
        and stores them as the class variables:
            self.onePlusFilenames and self.nPlusFilenames.'''
        
        onePlusFilenames = []
        nPlusFilenames = []
        
        splittedTextVar = textVar.split("\n")
        
        c = True
        i = 0
        # Check that an input is given
        if splittedTextVar[i] == "":
            c = False
        
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
        self.parsingSaveActive()
        
    def parsingSaveParameters(self):
        '''Save the parameters given.
        TODO! The validity of given parameters need to be checked!
        '''
        params = self.rootWindow.parameters
        params["onePlusFilenames"]=self.onePlusFilenames
        params["nPlusFilenames"]=self.nPlusFilenames
        params["header"]=self.header.get()
        params["footer"]=self.footer.get()
        params["separator"]=self.separator.get()
        params["pathToRawData"]=self.dataDir.get()
        params["species"]=self.species.get()
        # Create the list of charge states
        cState_i = int(self.cState_i.get())
        cState_f = int(self.cState_f.get())
        listOfStates = [i for i in range(cState_i, cState_f+1)]
        params["cStates"]=listOfStates
        self.parsingWindow.destroy()

class abcWindow:
    
    def __init__(self, rootWindow):
        '''Create the parsing dialog window attached to the rootWindow.'''
        self.u = Utils()
        self.abcFilename = tk.StringVar()
        self.species = tk.StringVar()
        self.cState_i = tk.StringVar()
        self.cState_f = tk.StringVar()
        self.h = tk.StringVar()
        self.rootWindow=rootWindow
        
    def createABCWindow(self):
        '''Create the layout'''
        abcWindow = tk.Toplevel(self.rootWindow.root)
        abcWindow.title("fitting parameters")
        abcWindow.minsize(500, 195)
        self.abcWindow=abcWindow
        
        
        u = self.u
        
        # Ion species
        u.createLabel(labelText="Ion species: ",
                                     row=0,
                                     column=0,
                                     sticky="ne",
                                     weight=(0,0),
                                     window=abcWindow)
        self.species, speciesBox = u.createTextBox(row=0,
                                                   column=1,
                                                   sticky="nw",
                                                   weight=(0,0),
                                                   window=abcWindow,
                                                   command=self.abcSaveActive) 
        speciesBox.config(width=5)
        
        # If the species was already given, use the given value.
        c = len(self.rootWindow.parameters["species"]) > 0
        if c: 
            speciesBox.config(state=tk.DISABLED)
            speciesBox.config(textvariable=self.species.set(
                value=self.rootWindow.parameters["species"])
                )
        
        # Available charge states
        availableLabel = u.createLabel(labelText="Available charge states",
                      row=1, 
                      column=0,
                      window=abcWindow,
                      sticky="ne",
                      weight=(0,1))
        CreateToolTip(availableLabel, "Transients for at least 5 charge states are required.\
                      the transients must correspond to consecutive charge states.\
                          Charge states must be in numerically ascending order.\
                              The order must be the same as for the filenames.")
        
        cStateLabel_i = u.createLabel(labelText="i: ", 
                      row=2, 
                      column=0,
                      sticky="ne",
                      weight=(0,0),
                      window=abcWindow)
        CreateToolTip(cStateLabel_i, "First measured charge state.")
        
        self.cState_i, cStateBox_i = u.createTextBox(row=2,
                        column=1,
                        sticky="nw",
                        weight=(0,1),
                        command=self.abcSaveActive, 
                        window=abcWindow)
        cStateBox_i.config(width=5)
        # If the cState was already given, use the given value.
        c = len(self.rootWindow.parameters["cStates"]) > 0
        if c: 
            cStateBox_i.config(state=tk.DISABLED)
            cStateBox_i.config(textvariable=self.cState_i.set(
                value=self.rootWindow.parameters["cStates"][0])
                )
        
        
        cStateLabel_f = u.createLabel(labelText="f: ", 
                      row=3, 
                      column=0,
                      sticky="ne",
                      weight=(0,0),
                      window=abcWindow)
        CreateToolTip(cStateLabel_f, "Last measured charge state.")
        
        self.cState_f, cStateBox_f = u.createTextBox(row=3,
                        column=1,
                        sticky="nw",
                        weight=(0,1),
                        command=self.abcSaveActive, 
                        window=abcWindow)
        cStateBox_f.config(width=5)
        # If the cState was already given, use the given value.
        c = len(self.rootWindow.parameters["cStates"]) > 0
        if c: 
            cStateBox_f.config(state=tk.DISABLED)
            cStateBox_f.config(textvariable=self.cState_f.set(
                value=self.rootWindow.parameters["cStates"][-1])
                )
            
        
        # Runge-Kutta method stepsize
        hLabel = u.createLabel(labelText="RK4-method stepsize (s): ",
                                      row=4,
                                      column=0,
                                      sticky="ne",
                                      weight=(0,0),
                                      window=abcWindow)
        CreateToolTip(hLabel, "Specify the stepsize for the RK4 algorithm,\
                      used for fitting the differential current balance equation\
                          to the measurement data.")
        self.h, hBox = u.createTextBox(row=4,
                                    column=1,
                                    columnspan=1,
                                    sticky="nw",
                                    weight=(0,1),
                                    window=abcWindow,
                                    command=self.abcSaveActive) 
        hBox.config(width=10)
        
        # filename for abc.csv
        filenameLabel = u.createLabel(labelText="abc filename: ",
                                      row=5,
                                      column=0,
                                      sticky="ne",
                                      weight=(0,0),
                                      window=abcWindow)
        CreateToolTip(filenameLabel, "Specify name under which to save abc file.\
                      No spaces. Ending must be .csv.\
                          E.g. abc_h1e-3.csv")
        self.abcFilename, filenameBox = u.createTextBox(row=5,
                                                    column=1,
                                                    columnspan=1,
                                                    sticky="new",
                                                    weight=(0,1),
                                                    window=abcWindow,
                                                    command=self.abcSaveActive) 
    
        # Save button
        self.saveButton = u.createButton(buttonText="Save", 
                                 row=1000, 
                                 column=2, 
                                 window=abcWindow,
                                 state=tk.DISABLED,
                                 sticky="se",
                                 weight=(0,0),
                                 command=self.abcSaveParameters)
        
    def abcSaveActive(self, *_):
        '''If all the fields have been filled (correctly OR incorrectly)
        change the state of the Save button to NORMAL. If the fields
        are cleared, set state back to DISABLED.'''
        lenSpecies = len(self.species.get())
        lenStates = len(self.cState_i.get())*len(self.cState_f.get())
        lenH = len(self.h.get())
        lenF = len(self.abcFilename.get())
        
        if (lenSpecies*lenStates*lenH*lenF > 0):
            self.saveButton.config(state=tk.NORMAL)
        else:
            self.saveButton.config(state=tk.DISABLED)
    
    def abcSaveParameters(self):
        '''Save the parameters given.
        TODO! The validity of given parameters need to be checked!
        '''
        params = self.rootWindow.parameters
        params["species"]=self.species.get()
        # Create the list of charge states
        cState_i = int(self.cState_i.get())
        cState_f = int(self.cState_f.get())
        listOfStates = [i for i in range(cState_i, cState_f+1)]
        params["cStates"]=listOfStates
        params["abcFilename"]=self.abcFilename.get()
        params["rkStepsize"]=float(self.h.get())
        self.abcWindow.destroy()

class optimizeWindow:
    
    def __init__(self, rootWindow):
        '''Create the optimization dialog window attached to the rootWindow.'''
        self.u = Utils()
        self.rootWindow=rootWindow
        
    def createOptWindow(self):
        '''Create the layout'''
        optWindow = tk.Toplevel(self.rootWindow.root)
        optWindow.title("Optimization parameters")
        optWindow.minsize(500, 300)
        u = self.u
        self.ne_lo = tk.StringVar()
        self.ne_hi = tk.StringVar()
        self.Ee_lo = tk.StringVar()
        self.Ee_hi = tk.StringVar()
        self.freq = tk.StringVar()
        self.N = tk.StringVar()
        self.MC = tk.StringVar()
        self.rc_method = tk.StringVar()
        self.abcFilename = tk.StringVar()
        self.optWindow=optWindow  
        
        # Save button
        self.saveButton = u.createButton(buttonText="Save", 
                                 row=1000, 
                                 column=3, 
                                 window=optWindow,
                                 state=tk.NORMAL,
                                 sticky="se",
                                 weight=(0,0),
                                 command=self.optSaveParameters)
        # abc file
        abcLabel = u.createLabel(labelText="abc file to use: ",
                                 row=0,
                                 column=0,
                                 sticky="ne",
                                 weight=(0,0),
                                 window=optWindow)
        self.abcFilename, abcBox = u.createTextBox(row=0,
                        column=1,
                        sticky="nw",
                        weight=(0,1),
                        command=self.optSaveActive, 
                        window=optWindow)
        abcBox.config(width=20)
        
        # If the cState was already given, use the given value.
        c = len(self.rootWindow.parameters["abcFilename"]) > 0
        if c: 
            abcBox.config(state=tk.DISABLED)
            abcBox.config(textvariable=self.abcFilename.set(
                value=self.rootWindow.parameters["abcFilename"])
                )
        
        # ne lower and upper limits
        neLabel = u.createLabel(labelText="Electron density limits",
                                   row=1,
                                   column=0,
                                   sticky="nw",
                                   weight=(0,1),
                                   window=optWindow)
        CreateToolTip(neLabel, "Units of 1/cm3")
        neLabel_lo = u.createLabel(labelText="Lower",
                                   row=2,
                                   column=0,
                                   sticky="e",
                                   weight=(0,0),
                                   window=optWindow)
        CreateToolTip(neLabel_lo, "Can be estimated using the 1+ stopping method\
                                    (doi: 10.1103/PhysRevAccelBeams.19.053402).\
                                    Default value corresponds to 500 W @ 14.5 GHz.")
        self.ne_lo, neBox_lo = u.createTextBox(row=2,
                                               column=1,
                                               sticky="w",
                                               weight=(0,1),
                                               window=optWindow,
                                               command=self.optSaveActive) 
        neBox_lo.config(width=10)
        self.ne_lo.set(value=f"{1e+11:.3e}")
        neLabel_hi = u.createLabel(labelText="Upper",
                                   row=3,
                                   column=0,
                                   sticky="e",
                                   weight=(0,0),
                                   window=optWindow)
        CreateToolTip(neLabel_hi, "Default is the cut-off frequency.")
        self.ne_hi, neBox_hi = u.createTextBox(row=3,
                                               column=1,
                                               sticky="we",
                                               weight=(0,0),
                                               window=optWindow,
                                               command=self.optSaveActive) 
        neBox_hi.config(width=10)
        
        # Cut-off density calculation
        cutoffLabel = u.createLabel(labelText="Calculate cut-off density",
                                         row=2, column=2, sticky="we",
                                         weight=(0,1), window=optWindow)
        CreateToolTip(cutoffLabel, "Input the heating frequency (Hz) below\
                      to automatically set the cut-off density \
                          as the upper limit.")
        u.createLabel(labelText="uW frequency: ",
                                  row=3,column=2,sticky="e",
                                  weight=(0,0), window=optWindow)
        self.freq, freqBox = u.createTextBox(row=3, column=3, sticky="we",
                                        weight=(0,1), window=optWindow, 
                                        command=self.calculateCutoff) 
        self.freq.set(value=14.5e9)
        # <Ee> lower and upper limits
        EeLabel = u.createLabel(labelText="Electron energy <Ee> limits",
                                   row=4,
                                   column=0,
                                   sticky="nw",
                                   weight=(0,1),
                                   window=optWindow)
        CreateToolTip(EeLabel, "Units of eV")
        EeLabel_lo = u.createLabel(labelText="Lower",
                                   row=5,
                                   column=0,
                                   sticky="e",
                                   weight=(0,0),
                                   window=optWindow)
        CreateToolTip(EeLabel_lo, "Default is 10 eV, i.e.\
                      order of magnitude of the plasma potential.")
        self.Ee_lo, EeBox_lo = u.createTextBox(row=5,
                                               column=1,
                                               sticky="w",
                                               weight=(0,1),
                                               window=optWindow,
                                               command=self.optSaveActive) 
        EeBox_lo.config(width=10)
        self.Ee_lo.set(value=10)
        EeLabel_hi = u.createLabel(labelText="Upper",
                                   row=6,
                                   column=0,
                                   sticky="e",
                                   weight=(0,0),
                                   window=optWindow)
        CreateToolTip(EeLabel_hi, "Default is 10 keV.\
                                    N.B. recommended cross section data \
                                    only goes up to 10 keV.")
        self.Ee_hi, EeBox_hi = u.createTextBox(row=6,
                                               column=1,
                                               sticky="we",
                                               weight=(0,0),
                                               window=optWindow,
                                               command=self.optSaveActive) 
        self.Ee_hi.set(value=9.999e3) 
        EeBox_hi.config(width=10)
        
        # Solution set settings
        u.createLabel(labelText="Solution set settings", row=7, column=0,
                      sticky="w", weight=(0,0), window=optWindow)
        # Number of elements in ne vector
        NLabel = u.createLabel(labelText="N",
                               row=8, column=0, sticky="e",
                               weight=(0,0), window=optWindow)
        CreateToolTip(NLabel, "Number of elements in ne-array.\
                      Determines the density of the solution space.\
                          Default is 1000.")
        self.N, NBox = u.createTextBox(row=8, column=1, sticky="we",
                                       weight=(0,0),window=optWindow,
                                       command=self.optSaveActive) 
        self.N.set(value=1000)
        
        # Number of Monte Carlo iterations
        MCLabel = u.createLabel(labelText="MC iterations",
                               row=9, column=0, sticky="e",
                               weight=(0,0), window=optWindow)
        CreateToolTip(MCLabel, "Number of Monte Carlo Iterations to perform.\
                      The bigger the number, the better the uncertainty \
                          of the rate coefficients is sampled.\
                          Default is 1000.")
        self.MC, MCBox = u.createTextBox(row=9, column=1, sticky="we",
                                       weight=(0,0),window=optWindow,
                                       command=self.optSaveActive) 
        self.MC.set(value=1000)
        
        

        
        
        # Rate coefficient evaluation method
        choices = {"Voronov", "Maxwell-Boltzmann"}
        self.rc_method.set("Maxwell-Boltzmann")
        popupLabel = tk.Label(optWindow, text="Rate coefficient evaluation: ")
        popupLabel.grid(row=10, column=0)
        popupMenu = tk.OptionMenu(optWindow, self.rc_method, *choices)
        popupMenu.grid(row=10, column=1)
        popupMenu.config(width=20)
        self.rc_method.trace("w", self.optSaveActive)
        
    def calculateCutoff(self, *_):
        '''Calculates the cut-off density for input into ne_hi.'''
        eps0 = 8.8542e-12
        e = 1.6021773e-19
        me = 9.10938356e-31
        f = float(self.freq.get())
        cutoff = eps0*me*(2*np.pi*f)**2/e**2
        cutoff = f"{cutoff:.3e}"
        self.ne_hi.set(value=cutoff)
        
    def optSaveActive(self, *_):
        '''If all the fields have been filled (correctly OR incorrectly)
        change the state of the Save button to NORMAL. If the fields
        are cleared, set state back to DISABLED.'''
        lenLo = len(self.ne_lo.get())*len(self.Ee_lo.get())
        lenHi = len(self.ne_hi.get())*len(self.Ee_hi.get())
        lenN = len(self.N.get())
        lenMC = len(self.MC.get())
        lenRC = len(self.rc_method.get())
        lenABC = len(self.abcFilename.get())
        if (lenLo*lenHi*lenN*lenMC*lenRC*lenABC > 0):
            self.saveButton.config(state=tk.NORMAL)
        else:
            self.saveButton.config(state=tk.DISABLED)
    
    def optSaveParameters(self):
        '''Save the parameters given.
        TODO! The validity of given parameters need to be checked!
        '''
        params = self.rootWindow.parameters
        params["abcFilename"]=self.abcFilename.get()
        params["ne_lo"]=float(self.ne_lo.get())
        params["ne_hi"]=float(self.ne_hi.get())
        params["Ee_lo"]=float(self.Ee_lo.get())
        params["Ee_hi"]=float(self.Ee_hi.get())
        # TODO! Set the method
        if self.rc_method.get() == "Maxwell-Boltzmann":
            params["method"]="interpMB"
        elif self.rc_method.get() == "Voronov":
            params["method"]="voronov"
        params["N"]=int(self.N.get())
        params["N_MC"]=int(self.MC.get())
        self.optWindow.destroy()
        
# Instantiate the window and run its mainloop
window = MainWindow()
window.root.mainloop()