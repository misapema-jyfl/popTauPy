#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 12:42:59 2021

Dialog for creating parameters .YAML files 
for the consecutive transients analyzer code.

@author: miha
"""

import tkinter as tk
import tkinter.filedialog as fd
import numpy as np 
import yaml 

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
            "one_plus_naming_convention":"",
            "n_plus_naming_convention":"",
            "onePlusFilenames":[],
            "nPlusFilenames":[],
            "header":0,
            "footer":0,
            "multiplying_factor":0,
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
            "method":"",
            "solution_set_files":[],
            "plotting":{
                "available_charge_states":[],
                "ne_lo":0,
                "ne_hi":0,
                "Ee_lo":0,
                "Ee_hi":0,
                "F_hi":0,
                "confidence":100,
                "data_point":"median",
                "plot_solution_sets":False,
                "plot_solutions_v_F":False,
                "plot_char_v_q":False,
                "plot_ec_v_q":False,
                "plot_triple_v_q":False
                }}
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
        self.workingDir.set(value="./")
        
        
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
        self.saveToPath.set(value="")
        u.createButton(buttonText="Set", 
                                row=1, column=2, 
                                command=self.setParameter,
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
        self.optButton = u.createButton(buttonText="(ne, Ee)-optimization", 
                          row=5,
                          column=1,
                          weight=(0,1),
                          sticky="new",
                          command=o.createOptWindow,
                          state=tk.DISABLED,
                          window=self.root)
        pl = PlottingWindow(self)
        self.plotButton = u.createButton(buttonText="Plotting",
                          row=6,
                          column=1,
                          weight=(0,1),
                          sticky="news",
                          command=pl.createPlottingWindow,
                          state=tk.DISABLED,
                          window=self.root)
        
        # Field for inputting yaml file filename
        self.yamlName = tk.StringVar(value="../params.yaml")
        yamlEntry = tk.Entry(self.root,
                             textvariable=self.yamlName)
        yamlEntry.grid(row=1000,column=1,sticky="e",
                       ipadx=5, ipady=5)
        self.yamlName.trace_add(mode="write", callback=u.dummy)
        
        self.button = u.createButton(buttonText="Make YAML", 
                                     row=1000,
                                     column=2,
                                     weight=(0,1),
                                     sticky="en",
                                     command=self.broadCastParameters,
                                     state=tk.NORMAL,
                                     window=self.root)
    
    def setParameter(self, *_):
        '''Set values of workingDir and SaveToPath'''
        self.parameters["workingDir"] = self.workingDir.get()
        self.parameters["saveToPath"] = self.saveToPath.get()
        
    def broadCastParameters(self, *_):
        
        
        plotting = self.parameters["plotting"] # plotting parameters 
        
        # Create the YAML file
        data = dict(
            general = dict(
                working_directory = self.parameters["workingDir"],
                elemental_data_directory = self.parameters["workingDir"] + "elemental_data/",
                save_to_path = self.parameters["saveToPath"],
                do_parse = self.parseVar.get(),
                do_abc = self.fitVar.get(),
                do_neTe = self.optVar.get(),
                do_plots = self.plotVar.get()
                ),
            parsing = dict(
                raw_data_path = self.parameters["pathToRawData"],
                parsed_data_path = self.parameters["saveToPath"],
                species = self.parameters["species"],
                available_charge_states = self.parameters["cStates"],
                header = self.parameters["header"],
                footer = self.parameters["footer"],
                multiplying_factor = self.parameters["multiplying_factor"],
                separator = self.parameters["separator"],
                one_plus_file_names = self.parameters["onePlusFilenames"],
                n_plus_file_names = self.parameters["nPlusFilenames"]
                ),
            data = dict(
                parsed_data_path = self.parameters["saveToPath"],
                available_charge_states = self.parameters["cStates"],
                species = self.parameters["species"],
                abc_file_name = self.parameters["abcFilename"]
                ),
            optimizer = dict(
                rk_time_step = self.parameters["rkStepsize"],
                rate_coefficient_method = self.parameters["method"],
                ne_lo = self.parameters["ne_lo"],
                ne_hi = self.parameters["ne_hi"],
                Ee_lo = self.parameters["Ee_lo"],
                Ee_hi = self.parameters["Ee_hi"],
                number_of_ne = self.parameters["N_ne"],
                number_of_MC = self.parameters["N_MC"]
                ),
            plotting = dict(
                solution_set_files = self.parameters["solution_set_files"],
                available_charge_states = plotting["available_charge_states"],
                ne_lo = plotting["ne_lo"],
                ne_hi = plotting["ne_hi"],
                Ee_lo = plotting["Ee_lo"],
                Ee_hi = plotting["Ee_hi"],
                F_hi = plotting["F_hi"],
                confidence = plotting["confidence"]/100, # Conversion from % -> 1
                plot_solution_sets = dict(
                    do_plot = plotting["plot_solution_sets"],
                    sigma = 10, # TODO!
                    x_scale = "log",
                    y_scale = "log"
                    ),
                plot_number_of_solutions_vs_F = dict(
                    do_plot = plotting["plot_solutions_v_F"],
                    list_of_F_values = [1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9], # TODO!
                    y_lo = 1,
                    y_hi = False
                    ),
                plot_characteristic_time_vs_charge = dict(
                    do_plot = plotting["plot_char_v_q"],
                    output_data = True,
                    conf_yscale = "linear", # TODO!
                    conf_marker = "o",
                    conf_color = "r",
                    inz_yscale = "linear",
                    inz_marker = "s",
                    inz_color = "k",
                    cx_yscale = "linear",
                    cx_marker = "^",
                    cx_color = "b"
                    ),
                plot_energy_content_vs_charge = dict(
                    do_plot = plotting["plot_ec_v_q"],
                    output_data = True, # TODO!
                    marker = "s",
                    color = "k",
                    y_hi = False,
                    y_lo = False,
                    ),
                plot_triple_product_vs_charge = dict(
                    do_plot = plotting["plot_triple_v_q"],
                    output_data = True, # TODO!
                    marker = "s",
                    color = "m",
                    y_hi = False,
                    y_lo = False
                    )
                )
            )
        
        filename = str(self.yamlName.get())
        with open(filename, "w") as outfile:
            yaml.dump(data, outfile, default_flow_style=False)
            
class ParsingWindow:
    
    def __init__(self, rootWindow):
        '''Create the parsing dialog window attached to the rootWindow.'''
        self.u = Utils()
        self.dataDir = tk.StringVar()
        self.header = tk.StringVar()
        self.footer = tk.StringVar()
        self.separator = tk.StringVar()
        self.factor = tk.StringVar()
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
        # Multiplying factor
        parsingFactorLabel = u.createLabel(labelText="Conversion factor:", 
                         row=7,
                         column=0,
                         sticky="nw",
                         weight=(0,0),
                         window=parsingWindow)
        CreateToolTip(parsingFactorLabel, "Conversion factor required\
                      to convert the transient signal to units of Amperes.")
        self.factor, factorBox = u.createTextBox(row=7,
                        column=1,
                        sticky="nw",
                        weight=(0,1),
                        window=parsingWindow,
                        command=self.parsingSaveActive)
        factorBox.config(width=10)
        # Ion species and charge states
        speciesLabel = u.createLabel(labelText="Ion species: ", 
                      row=8, 
                      column=0,
                      sticky="nw",
                      weight=(0,0),
                      window=parsingWindow)
        CreateToolTip(speciesLabel, "Name of element, e.g. \'K\'.")
        self.species, speciesBox = u.createTextBox(row=8,
                        column=1,
                        sticky="nw",
                        weight=(0,0),
                        command=self.parsingSaveActive,
                        window=parsingWindow)
        speciesBox.config(width=5)
        
        availableLabel = u.createLabel(labelText="Available charge states",
                      row=9, 
                      column=0,
                      window=parsingWindow,
                      sticky="nw",
                      weight=(0,0))
        CreateToolTip(availableLabel, "Transients for at least 5 charge states are required.\
                      the transients must correspond to consecutive charge states.\
                          Charge states must be in numerically ascending order.\
                              The order must be the same as for the filenames.")
        
        cStateLabel_i = u.createLabel(labelText="i: ", 
                      row=10, 
                      column=0,
                      sticky="ne",
                      weight=(0,0),
                      window=parsingWindow)
        CreateToolTip(cStateLabel_i, "First measured charge state.")
        
        self.cState_i, cStateBox_i = u.createTextBox(row=10,
                        column=1,
                        sticky="nw",
                        weight=(0,0),
                        command=self.parsingSaveActive,
                        window=parsingWindow)
        cStateBox_i.config(width=5)
        
        cStateLabel_f = u.createLabel(labelText="f: ", 
                      row=11, 
                      column=0,
                      sticky="ne",
                      weight=(0,0),
                      window=parsingWindow)
        CreateToolTip(cStateLabel_f, "Last measured charge state.")
        
        self.cState_f, cStateBox_f = u.createTextBox(row=11,
                        column=1,
                        sticky="nw",
                        weight=(0,0),
                        command=self.parsingSaveActive,
                        window=parsingWindow)
        cStateBox_f.config(width=5)
        
        # Create labels for 1+ and n+ filenames
        u.createLabel("1+ filenames:",
                         row=12,
                         column=0,
                         sticky="nwe",
                         weight=(0,0),
                         window=parsingWindow)
        u.createLabel("N+ filenames:",
                         row=12,
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
        lenFactor = len(self.factor.get())
        if (lenDir*lenOnePlus*lenNPlus\
            *lenH*lenF*lenSep*lenSpecies*lenStates*lenFactor > 0):
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
                             row = i+13,
                             column = 0,
                             sticky = "wne",
                             weight = (0,0),
                             window = self.parsingWindow)
            self.parsingFilenameLabels["1+"].append(label)
            
        # Print the N+ filenames
        for i, name in enumerate(self.nPlusFilenames):
            
            label = u.createLabel(labelText=name,
                             row = i+13,
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
        params["one_plus_naming_convention"]=self.onePlusConvention.get()
        params["n_plus_naming_convention"]=self.nPlusConvention.get()
        params["header"]=int(self.header.get())
        params["footer"]=int(self.footer.get())
        params["multiplying_factor"]=float(self.factor.get())
        params["separator"]=str(self.separator.get())
        params["pathToRawData"]=str(self.dataDir.get())
        params["species"]=str(self.species.get())
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
        params["species"]=str(self.species.get())
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
        CreateToolTip(abcLabel, "The file is expected to be located in the\
                      save to path.")
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
        
        
        
        
class PlottingWindow:
    
    def __init__(self, rootWindow):
        '''Create the parsing dialog window attached to the rootWindow.'''
        self.u = Utils()
        self.ne_lo = tk.StringVar()
        self.ne_hi = tk.StringVar()
        self.Ee_lo = tk.StringVar()
        self.Ee_hi = tk.StringVar()
        self.F_hi = tk.StringVar()
        self.confidence = tk.StringVar()
        self.dataPoint = tk.StringVar()
        
        self.rootWindow=rootWindow
        
        
    
    def createPlottingWindow(self):
        '''Generates the plotting window.'''
        u = self.u
        
        window = tk.Toplevel(self.rootWindow.root)
        window.title("Plotting parameters")
        window.minsize(500, 350)
        self.window=window
        
        # Save button
        self.saveButton = u.createButton(buttonText="Save", 
                                 row=1000, 
                                 column=3, 
                                 window=window,
                                 state=tk.NORMAL,
                                 sticky="se",
                                 weight=(0,0),
                                 command=self.saveParameters)
        
        
        # Available charge states
        availableLabel = u.createLabel(labelText="Designate by charge states",
                      row=0, 
                      column=0,
                      window=window,
                      sticky="ne",
                      weight=(0,1))
        CreateToolTip(availableLabel, "Specify charge states for which \
                                      solution sets are available.")
        
        # If the charge states are already known,
        # designate the solution sets according to the 
        # naming convention.
        # TODO! Note that the optimization dialog needs to be
        # already filled by this point.
        if len(self.rootWindow.parameters["cStates"])>0:
            self.designate_solution_sets()
        
        
        cStateLabel_i = u.createLabel(labelText="i: ", 
                      row=1, 
                      column=0,
                      sticky="ne",
                      weight=(0,0),
                      window=window)
        CreateToolTip(cStateLabel_i, "First measured charge state.")
        
        self.cState_i, cStateBox_i = u.createTextBox(row=1,
                        column=1,
                        sticky="nw",
                        weight=(0,1),
                        command=self.saveActive, 
                        window=window)
        cStateBox_i.config(width=5)
        # If the cState was already given, use the given value.
        c = len(self.rootWindow.parameters["cStates"]) > 0
        if c: 
            cStateBox_i.config(state=tk.DISABLED)
            cStateBox_i.config(textvariable=self.cState_i.set(
                value=self.rootWindow.parameters["cStates"][0]+2)
                )
        
        
        cStateLabel_f = u.createLabel(labelText="f: ", 
                      row=2, 
                      column=0,
                      sticky="ne",
                      weight=(0,0),
                      window=window)
        CreateToolTip(cStateLabel_f, "Last measured charge state.")
        
        self.cState_f, cStateBox_f = u.createTextBox(row=2,
                        column=1,
                        sticky="nw",
                        weight=(0,1),
                        command=self.saveActive, 
                        window=window)
        cStateBox_f.config(width=5)
        # If the cState was already given, use the given value.
        c = len(self.rootWindow.parameters["cStates"]) > 0
        if c: 
            cStateBox_f.config(state=tk.DISABLED)
            cStateBox_f.config(textvariable=self.cState_f.set(
                value=self.rootWindow.parameters["cStates"][-1]-2)
                )
        
        # Designate manually 
        openButton = tk.Button(window, text="Manually select files to use.", 
                               command=self.selectFiles)
        openButton.grid(row=2, column=2, sticky="e")
        
        # If the cStates were already given, disable the manual choice.
        c = len(self.rootWindow.parameters["cStates"]) > 0
        if c: 
            openButton.config(state=tk.DISABLED)
        
        
        u.createLabel(labelText="Set limits on results to consider in plots: ",
                      row=3, column=0, sticky="w",
                      weight=(0,0), window=window)
        
        
        neLabel_lo = u.createLabel(labelText="ne lower",
                                   row=4,
                                   column=0,
                                   sticky="e",
                                   weight=(0,0),
                                   window=window)
        CreateToolTip(neLabel_lo, "Can be estimated using the 1+ stopping method\
                                    (doi: 10.1103/PhysRevAccelBeams.19.053402).\
                                    Default value corresponds to 500 W @ 14.5 GHz.")
        self.ne_lo, neBox_lo = u.createTextBox(row=4,
                                               column=1,
                                               sticky="w",
                                               weight=(0,1),
                                               window=window,
                                               command=self.saveActive) 
        neBox_lo.config(width=10)
        self.ne_lo.set(value=f"{1e+11:.3e}")
        neLabel_hi = u.createLabel(labelText="ne upper",
                                   row=5,
                                   column=0,
                                   sticky="e",
                                   weight=(0,0),
                                   window=window)
        CreateToolTip(neLabel_hi, "Default is the cut-off frequency.")
        self.ne_hi, neBox_hi = u.createTextBox(row=5,
                                               column=1,
                                               sticky="we",
                                               weight=(0,0),
                                               window=window,
                                               command=self.saveActive) 
        neBox_hi.config(width=10)
        
        # Cut-off density calculation
        cutoffLabel = u.createLabel(labelText="Calculate cut-off density",
                                         row=4, column=2, sticky="we",
                                         weight=(0,1), window=window)
        CreateToolTip(cutoffLabel, "Input the heating frequency (Hz) below\
                                  to automatically set the cut-off density \
                                  as the upper limit.")
        u.createLabel(labelText="uW frequency: ",
                                  row=5,column=2,sticky="e",
                                  weight=(0,0), window=window)
        self.freq, freqBox = u.createTextBox(row=5, column=3, sticky="we",
                                        weight=(0,1), window=window, 
                                        command=self.calculateCutoff) 
        self.freq.set(value=14.5e9)
        # <Ee> lower and upper limits
        # EeLabel = u.createLabel(labelText="Electron energy <Ee> limits",
        #                            row=4,
        #                            column=0,
        #                            sticky="nw",
        #                            weight=(0,1),
        #                            window=window)
        # CreateToolTip(EeLabel, "Units of eV")
        EeLabel_lo = u.createLabel(labelText="<Ee> lower",
                                   row=6,
                                   column=0,
                                   sticky="e",
                                   weight=(0,0),
                                   window=window)
        CreateToolTip(EeLabel_lo, "Default is 10 eV, i.e.\
                      order of magnitude of the plasma potential.")
        self.Ee_lo, EeBox_lo = u.createTextBox(row=6,
                                               column=1,
                                               sticky="w",
                                               weight=(0,1),
                                               window=window,
                                               command=self.saveActive) 
        EeBox_lo.config(width=10)
        self.Ee_lo.set(value=10)
        EeLabel_hi = u.createLabel(labelText="<Ee> upper",
                                   row=7,
                                   column=0,
                                   sticky="e",
                                   weight=(0,0),
                                   window=window)
        CreateToolTip(EeLabel_hi, "Default is 10 keV.\
                                    N.B. recommended cross section data \
                                    only goes up to 10 keV.")
        self.Ee_hi, EeBox_hi = u.createTextBox(row=7,
                                               column=1,
                                               sticky="we",
                                               weight=(0,0),
                                               window=window,
                                               command=self.saveActive) 
        self.Ee_hi.set(value=9.999e3) 
        EeBox_hi.config(width=10)
        
        # Penalty function upper limit
        fLabel = u.createLabel(labelText="Upper limit of penalty function: ", 
                               row=8, column=0, sticky="e", weight=(0,0),
                               window=window)
        # CreateToolTip(fLabel, "") # TODO!
        self.F_hi, FBox_hi = u.createTextBox(row=8,
                                               column=1,
                                               sticky="we",
                                               weight=(0,0),
                                               window=window,
                                               command=self.saveActive) 
        self.F_hi.set(value=1e-6) 
        FBox_hi.config(width=10)
        
        u.createLabel(labelText="Uncertainty bar estimation: ",
                      row=9, column=0, sticky="w",
                      weight=(0,0), window=window)
        
        confidenceLabel = u.createLabel(labelText="Confidence %: ",
                                        row=10, column=0, sticky="e",
                                        weight=(0,0), window=window)
        CreateToolTip(confidenceLabel, "Percentage of results that the \
                      uncertainty bars will enclose below and above \
                          the data point. (34.1 % corresponds to 1 sigma)")
        self.confidence, confBox = u.createTextBox(row=10,
                                               column=1,
                                               sticky="we",
                                               weight=(0,0),
                                               window=window,
                                               command=self.saveActive) 
        self.confidence.set(value=34.1) 
        confBox.config(width=10)
        
        # Data point selection
        choices = {"Median"} # TODO! Most probable value...
        self.dataPoint.set("Median")
        popupLabel = tk.Label(window, text="Data point: ")
        popupLabel.grid(row=11, column=0, sticky="e")
        popupMenu = tk.OptionMenu(window, self.dataPoint, *choices)
        popupMenu.grid(row=11, column=1)
        popupMenu.config(width=20)
        self.dataPoint.trace("w", self.saveActive)
        
        # Create check boxes for things to plot
        u.createLabel(labelText="Things to plot:",
                         row=12,
                         column=0,
                         sticky="w",
                         weight=(0,0),
                         window=window)
        self.plotSolSets = u.createCheckBox(boxText=" Solution sets",
                            row=13,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            window=window)
        self.plotSolSets.set(value=1)
        self.plotVsF = u.createCheckBox(boxText=" Number of solutions vs. F",
                            row=14,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            window=window)
        self.plotVsF.set(value=1)
        self.plotChar = u.createCheckBox(boxText=" Characteristic times vs. q",
                            row=15,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            window=window)
        self.plotChar.set(value=1)
        self.plotEC = u.createCheckBox(boxText=" Energy content",
                            row=16,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            window=window)
        self.plotEC.set(value=1)
        self.plotTriple = u.createCheckBox(boxText="Triple product",
                            row=17,
                            column=0,
                            weight=(0,0),
                            sticky="nw",
                            window=window)
        self.plotTriple.set(value=1)
        
        filesToUseLabel = tk.Label(window, text="Solution set files to use:")
        filesToUseLabel.grid(row=18, column=0, columnspan=4, sticky="we")
        # If the cStates were already given
        c = len(self.rootWindow.parameters["cStates"]) > 0
        if c: 
            l = tk.Label(window, text="AUTOMATICALLY CHOSEN")
            l.grid(row=19, column=0, columnspan=4, sticky="we")
            
            
    def designate_solution_sets(self):
        '''If the charge states have already been specified
        automatically designates the solution set file filenames
        according to the solution set naming convention.'''
        params = self.rootWindow.parameters
        ci = params["cStates"][0]
        cf = params["cStates"][-1]
        N = str(params["N_ne"])
        MC = str(params["N_MC"])
        self.solution_set_files = []
        for i in range(ci+2, cf-2+1):
            s = (params["saveToPath"],"solset_MC_iters-",MC,"_N-",N,"_q-",str(i),".csv")
            fname = "".join(s)
            self.solution_set_files.append(fname)
        
    
    def selectFiles(self, *_):
        '''Manually select solution set files to use in plotting.'''
        filez = fd.askopenfilenames(parent=self.window, title='Choose a file')
        self.solution_set_files = []
        
        for i, file in enumerate(filez):
            label = tk.Label(self.window, text=file)
            label.grid(row=19+i, column=0, columnspan=4, sticky="we")
            self.solution_set_files.append(file)
    
    def calculateCutoff(self, *_):
        '''Calculates the cut-off density for input into ne_hi.'''
        eps0 = 8.8542e-12
        e = 1.6021773e-19
        me = 9.10938356e-31
        f = float(self.freq.get())
        cutoff = eps0*me*(2*np.pi*f)**2/e**2
        cutoff = f"{cutoff:.3e}"
        self.ne_hi.set(value=cutoff)
    
    def saveActive(self, *_):
        '''If all the fields have been filled (correctly OR incorrectly)
        change the state of the Save button to NORMAL. If the fields
        are cleared, set state back to DISABLED.'''
        lenLo = len(self.ne_lo.get())*len(self.Ee_lo.get())
        lenHi = len(self.ne_hi.get())*len(self.Ee_hi.get())
        
        if (lenLo*lenHi > 0):
            self.saveButton.config(state=tk.NORMAL)
        else:
            self.saveButton.config(state=tk.DISABLED)
            
    def saveParameters(self):
        '''Save the parameters given.
        TODO! The validity of given parameters need to be checked!
        '''
        ci = int(self.cState_i.get())
        cf = int(self.cState_f.get())
        self.solset_cStates = [i for i in range(ci,cf+1)]
        params = self.rootWindow.parameters
        params["solution_set_files"] = self.solution_set_files
        params["plotting"]["ne_lo"] = float(self.ne_lo.get())
        params["plotting"]["ne_hi"] = float(self.ne_hi.get())
        params["plotting"]["Ee_lo"] = float(self.Ee_lo.get())
        params["plotting"]["Ee_hi"] = float(self.Ee_hi.get())
        params["plotting"]["F_hi"] = float(self.F_hi.get())
        params["plotting"]["confidence"] = float(self.confidence.get())
        params["plotting"]["data_point"] = str(self.dataPoint.get())
        params["plotting"]["plot_solution_sets"] = int(self.plotSolSets.get())
        params["plotting"]["plot_solutions_v_F"] = int(self.plotVsF.get())
        params["plotting"]["plot_char_v_q"] = int(self.plotChar.get())
        params["plotting"]["plot_ec_v_q"] = int(self.plotEC.get())
        params["plotting"]["plot_triple_v_q"] = int(self.plotTriple.get())
        params["plotting"]["available_charge_states"] = self.solset_cStates
        self.window.destroy()
        
# Instantiate the window and run its mainloop
window = MainWindow()
window.root.mainloop()