# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 12:17:04 2022

@author: Ben
"""

import os
from tkinter import *
import PySimpleGUI as sg
from pathlib import Path

def imageViewerGUI(initialDir, defaultDir):
    sg.theme("DarkGreen1")
    
    
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
        [sg.T("         "), sg.Checkbox('Show scale bar?', default=True, key="input1")],
        [sg.T("         ")],
        [sg.T("         "), sg.Checkbox('Empty', default=False, key="input2")],
        [sg.T("         ")],
        [sg.Text('Scale bar size (um)', size =(15, 1)), sg.InputText(size =(10, 1),enable_events=True, key="scaleBarSize")],
        [sg.T("         ")],
        [sg.Button("Run")]
    ]
    
    layout = [[
        sg.Column(file_list_column),
        sg.VSeperator(),
        sg.Column(config_column)]]

    
    # Create the window
    window = sg.Window("Image Viewer", layout )
    
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
              
         
    if values["input1"] == True:
        print('\n')
    else:
        print('')
        
    if values["input2"] == True:
        print('\n')
    else:
        print(' ')
        
        
    scaleBarSize = values["scaleBarSize"]
    if not scaleBarSize:
         scaleBarSize = 20
         print('\n Using default scaleBarSize ' + str(scaleBarSize) + ' um')
    else:
        scaleBarSize = float(scaleBarSize)
        print('\n Using user defined scaleBarSize ' + str(scaleBarSize) + ' um')
      
      
    save_path = data_path + '\\SHG tiffs'
    shg_path = Path(save_path)

    if not os.path.exists(save_path):
        os.mkdir(save_path)
      
    save_path = data_path + '\\SHG pngs'
    shg_png_path = Path(save_path)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

        
    save_path = data_path + '\\TPEF tiffs'
    tpf_path = Path(save_path)

    if not os.path.exists(save_path):
        os.mkdir(save_path)
        
    save_path = data_path + '\\TPEF pngs'
    tpf_png_path = Path(save_path)

    if not os.path.exists(save_path):
        os.mkdir(save_path)
        
    save_path = data_path + '\\Z stack'
    zstackPath = Path(save_path)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    save_path = data_path + '\\Z stack' + '\\SHG pngs'
    shg_png_zstackPath = Path(save_path)   
    if not os.path.exists(save_path):
        os.mkdir(save_path)
            
    save_path = data_path + '\\Z stack' + '\\TPEF pngs'
    tpf_png_zstackPath = Path(save_path)   
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    
    return data_path, values["input1"], values["input2"], scaleBarSize , shg_path, shg_png_path,tpf_path,tpf_png_path, zstackPath, shg_png_zstackPath,tpf_png_zstackPath 