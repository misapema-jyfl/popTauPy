#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 12:42:59 2021

@author: miha
"""

import tkinter as tk
from utilities import Utils

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
        self.parameters = {"species":"",
                                  "onePlusFilenames":[],
                                  "nPlusFilenames":[],
                                  "cStates":[],
                                  "header":0,
                                  "footer":0,
                                  "separator":",",
                                  "pathToRawData":""}
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
        self.workingDirInp = u.createTextBox(row=0,
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
        self.saveToPathInp = u.createTextBox(row=1,
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
        self.fitButton = u.createButton(buttonText="abc fitting", # TODO!
                          row=4,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=a.createABCWindow,
                          state=tk.DISABLED,
                          window=self.root)
        self.optButton = u.createButton(buttonText="(ne, Ee)-optimization", # TODO!
                          row=5,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=u.dummy,
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
        print("Parsing parameters:")
        print(self.parameters)
            
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
        self.rootWindow=rootWindow
        self.species = tk.StringVar()
        
    def createABCWindow(self):
        '''Create the layout'''
        abcWindow = tk.Toplevel(self.rootWindow.root)
        abcWindow.title("fitting parameters")
        abcWindow.minsize(500, 195)
        self.abcWindow=abcWindow
        
        
        u = self.u
        
        # Ion species
        speciesLabel = u.createLabel(labelText="Ion species: ",
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
                                                   command=u.dummy) # TODO!
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
                        command=u.dummy, # TODO!
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
                        command=u.dummy, # TODO!
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
                                                        command=u.dummy) # TODO!
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
                                                        command=u.dummy) # TODO!
        
        # Save button
        saveABC = u.createButton(buttonText="Save", 
                                 row=1000, 
                                 column=2, 
                                 window=abcWindow,
                                 state=tk.DISABLED,
                                 sticky="se",
                                 weight=(0,0),
                                 command=u.dummy) # TODO!
        
        
# Instantiate the window and run its mainloop
window = MainWindow()
window.root.mainloop()