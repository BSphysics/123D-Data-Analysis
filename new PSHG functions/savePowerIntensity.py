# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 11:38:52 2022

@author: Ben
"""

import numpy as np

def savePowerIntensity(pshg,ptpf,data_path, filenames):
    
    IntensityvsAngle = np.mean(np.mean(pshg,1),1)
    meanIntensityvsAngle = IntensityvsAngle / np.mean(IntensityvsAngle)
    np.save(data_path + '/meanIntensityvsAngle' + filenames[0][0:3], meanIntensityvsAngle )
        
    pTPFvsAngle = np.mean(np.mean(ptpf,1),1)
    meanpTPFvsAngle = pTPFvsAngle /np.mean(pTPFvsAngle)
    np.save(data_path + '/meanpTPFvsAngle' + filenames[0][0:3], meanpTPFvsAngle )    
    