# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 11:02:52 2022

@author: Ben
"""
import numpy as np


def saveDepthIntensity(pshg, ptpf, data_path):

    IntensityvsAngle = np.mean(np.mean(pshg,1),1)
    meanIntensityvsAngle = IntensityvsAngle / np.mean(IntensityvsAngle)
    np.save(data_path + '/meanIntensityvsAngle' + data_path[-3:], meanIntensityvsAngle )
    
    
    pTPFvsAngle = np.mean(np.mean(ptpf,1),1)
    meanpTPFvsAngle = pTPFvsAngle /np.mean(pTPFvsAngle)
    np.save(data_path + '/meanpTPFvsAngle' + data_path[-3:], meanpTPFvsAngle )