#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 12:42:59 2021

@author: miha
"""

import tkinter as tk

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
    
    
class MainWindow:
    
    def __init__(self):
        
        self.root = tk.Tk()
        self.root.title("Parameters YAML file generator")
        self.root.minsize(500,220)
        self.createLayout()
    
    def createLayout(self):
        '''Creates the MainWindow layout'''
        u = Utils()
        
        # Specify working directory
        u.createLabel(labelText="Working directory:",
                            row=0,
                            column=0,
                            weight=(0,0),
                            sticky="n",
                            window=self.root)
        self.workingDirInp = u.createTextBox(row=0,
                                                column=1,
                                                columnspan=1,
                                                weight=(0,1),
                                                sticky="nwe",
                                                window=self.root)
        
        # Specify save to path
        u.createLabel(labelText="Save to path:",
                         row=1,
                         column=0,
                         weight=(0,0),
                         sticky="nwe",
                         window=self.root)
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
        p = ParsingWindow(self.root)
        self.parseButton = u.createButton(buttonText="Parsing",
                          row=3,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=p.createParsingWindow,
                          state=tk.DISABLED,
                          window=self.root)
        self.fitButton = u.createButton(buttonText="abc fitting",
                          row=4,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=u.dummy,
                          state=tk.DISABLED,
                          window=self.root)
        self.optButton = u.createButton(buttonText="(ne, Ee)-optimization",
                          row=5,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=u.dummy,
                          state=tk.DISABLED,
                          window=self.root)
        self.plotButton = u.createButton(buttonText="Plotting",
                          row=6,
                          column=1,
                          weight=(0,1),
                          sticky="news",
                          command=u.dummy,
                          state=tk.DISABLED,
                          window=self.root)
        
class ParsingWindow:
    
    def __init__(self, rootWindow):
        '''Create the parsing dialog window attached to the rootWindow.'''
        self.u = Utils()
        self.dataDir = tk.StringVar()
        self.header = tk.StringVar()
        self.footer = tk.StringVar()
        self.separator = tk.StringVar()
        self.parsingFilenameLabels = {"1+":[],"n+":[]}
        self.rootWindow=rootWindow
        
    
    def createParsingWindow(self):
        
        u = Utils()
        
        # To track the state of the Save button
        saveIsActive = tk.IntVar(0)
        
        parsingWindow = tk.Toplevel(self.rootWindow)
        parsingWindow.title("Parsing parameters")
        parsingWindow.minsize(500, 220)
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
        
        # Create labels for 1+ and n+ filenames
        u.createLabel("1+ filenames:",
                         row=7,
                         column=0,
                         sticky="nwe",
                         weight=(0,0),
                         window=parsingWindow)
        u.createLabel("N+ filenames:",
                         row=7,
                         column=1,
                         sticky="nwe",
                         weight=(0,1),
                         window=parsingWindow)
        # Create Save button
        self.saveButton = u.createButton(buttonText="Save",
                       row=1000,
                       column=2,
                       command=u.dummy,
                       sticky="news",
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
        if (lenDir*lenOnePlus*lenNPlus*lenH*lenF*lenSep > 0):
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
                             row = i+8,
                             column = 0,
                             sticky = "wne",
                             weight = (0,0),
                             window = self.parsingWindow)
            self.parsingFilenameLabels["1+"].append(label)
            
        # Print the N+ filenames
        for i, name in enumerate(self.nPlusFilenames):
            
            label = u.createLabel(labelText=name,
                             row = i+8,
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
        TODO! The validity of given parameters need to be checked!'''
        parsingParams = {}
        parsingParams["onePlusFilenames"]=self.onePlusFilenames
        parsingParams["nPlusFilenames"]=self.nPlusFilenames
        parsingParams["header"]=self.header
        parsingParams["footer"]=self.footer
        parsingParams["separator"]=self.separator
        parsingParams["pathToRawData"]=self.dataDir
        self.rootWindow.parsingParameters=parsingParams
        
# Instantiate the window and run its mainloop
window = MainWindow()
window.root.mainloop()