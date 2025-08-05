# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:23:46 2022

@author: Ben
"""
from rectROISelector import region_selector
from pathlib import Path
import os
import numpy as np
from pSHG_histograms import pSHGhistograms


def selectROI(phi2, I2, normSumDownscaleSHG, data_path, pshg, filenames):
   
    sPath = data_path + '\\ROIs'
    save_path = Path(sPath)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
        print("Directory " , save_path ,  " Created ")
    else:    
        print("Directory " , save_path ,  " already exists")
            
 
    roi = np.asarray(region_selector(np.sum(pshg,0), 0, np.max(np.sum(pshg,0)), 'SELECT ROI', sPath, filenames[0][:-10]))
    binned_roi = np.round(roi/4).astype('int')

    pSHGhistograms(phi2[binned_roi[2]:binned_roi[3], binned_roi[0]:binned_roi[1]], normSumDownscaleSHG[binned_roi[2]:binned_roi[3], binned_roi[0]:binned_roi[1]], 'Phi2', sPath)
    pSHGhistograms(I2[binned_roi[2]:binned_roi[3], binned_roi[0]:binned_roi[1]], normSumDownscaleSHG[binned_roi[2]:binned_roi[3], binned_roi[0]:binned_roi[1]], 'I2', sPath)


    return binned_roi