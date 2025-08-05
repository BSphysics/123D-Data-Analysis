# -*- coding: utf-8 -*-
"""
Created on Wed May 25 15:42:25 2022

@author: Ben
"""

from tkinter import *
import PySimpleGUI as sg

def pSHGGUI(initialDir, defaultDir, defaultThreshold):
    sg.theme("DarkBlue15")
    
    
    file_list_column = [
        [
            sg.Text("Select data folder"),
    
            sg.In(size=(25, 5), enable_events=True, key="folder"),
    
            sg.FolderBrowse(initial_folder = initialDir),
        ]  

    ]
    
    config_column = [
        [sg.Text("Select options for the pSHG analysis")], 
        [sg.T("         ")],
        [sg.T("         "), sg.Checkbox('Whole image histograms?', default=False, key="input1")],
        [sg.T("         ")],
        [sg.T("         "), sg.Checkbox('Show single pixel fits?', default=False, key="input2")],
        [sg.T("         ")],
        [sg.T("         "), sg.Checkbox('Polar histogram plot?',  default=False, key="input3")],
        [sg.T("         ")],
        [sg.T("         "), sg.Checkbox('Show Arrow plot?',  default=True, key="input4")],
        [sg.T("         ")],
        [sg.T("         "), sg.Checkbox('Use lasso ROI?',  default=False, key="input5")],
        [sg.T("         ")],
        [sg.T("         "), sg.Checkbox('Use slider intensity threshold?',  default=True, key="input6")],
        [sg.T("         ")],
        [sg.Text('Intensity threshold', size =(15, 1)), sg.InputText(size =(10, 1),enable_events=True, key="thresh")],
        [sg.T("         ")],
        [sg.Button("Run")]
    ]
    
    layout = [[
        sg.Column(file_list_column),
        sg.VSeperator(),
        sg.Column(config_column)]]

    
    # Create the window
    window = sg.Window("pSHG setup", layout )
    
    # Create an event loop
    while True:
        event, values = window.read()
        # End program if user closes window or
        # presses the OK button
        
        
        if event == "Run" or event == sg.WIN_CLOSED:
            break
    
    window.close()
    
    data_path = values["folder"]
    if not data_path:
        data_path = defaultDir
        print('\n Using default directory (see below)')
        print('\n' + data_path)
        
    else:
        print('\n Using user defined directory (see below)')
        print('\n' + data_path)
        
    thresh = values["thresh"]
    if not thresh and values["input6"] == False:
        threshold = defaultThreshold
        print('\n Using default threshold i.e. ' + str(threshold))
    elif not thresh and values["input6"] == True:
        threshold = 0.0
        print('\n Using slider threshold')
    else:
        threshold = float(thresh)
        print('\n Using user defined threshold i.e. ' + str(threshold))
        
         
    if values["input1"] == True:
        print('\n Plotting whole image histograms')
    else:
        print('')
        
    if values["input2"] == True:
        print('\n Use the pSHG fit viewer')
    else:
        print(' ')
        
    if values["input3"] == True:
        print('\n Plotting the polar histogram')
    else:
        print(' ')        
       
    return data_path, values["input1"], values["input2"], values["input3"], values["input4"], values["input5"], values["input6"], threshold