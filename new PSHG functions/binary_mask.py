# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 22:00:13 2018

@author: Ben
"""

import numpy as np

def binary_mask(threshLow, threshHigh ,im):
    
    img = np.copy(im)
    threshLow_idx = img < threshLow
   
    threshMid_idx = img >= threshLow
    
    threshHigh_idx = img > threshHigh
   
    img[threshLow_idx] = 0
    img[threshMid_idx] = 1
    img[threshHigh_idx] = 0
    

    
    return img