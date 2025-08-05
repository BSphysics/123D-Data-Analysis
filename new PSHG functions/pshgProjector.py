# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:05:17 2022

@author: Ben
"""
from skimage.morphology import disk
from skimage.filters import median
import numpy as np
from scipy.signal import convolve2d

roi = np.ones((3,3))
def pshgProjector(data_path, pshg, angles):
    
    anglesArr = np.zeros((len(pshg),pshg.shape[1],pshg.shape[2]))

    for idx in range(len(angles)):
        anglesArr[idx,:,:] = np.ones((pshg.shape[1],pshg.shape[2]))*angles[idx]

    cos2 = np.cos(2*anglesArr)
    cos4 = np.cos(4*anglesArr)
    sin2 = np.sin(2*anglesArr)
    sin4 = np.sin(4*anglesArr)

    a0 = 1*np.mean(pshg,0) + 1e-10
    a2 = 2*np.mean(pshg*cos2,0)
    a4 = 2*np.mean(pshg*cos4,0)
    b2 = 2*np.mean(pshg*sin2,0)
    b4 = 2*np.mean(pshg*sin4,0)
    
    # a0c = median(a0, disk(3)) 
    # a2c = median(a2, disk(3)) / a0c
    # a4c = median(a4, disk(3)) / a0c
    # b2c = median(b2, disk(3)) / a0c
    # b4c = median(b4, disk(3)) / a0c
    
    a0c = convolve2d(a0,roi, mode = 'same')  
    a2c = convolve2d(a2,roi, mode = 'same') / a0c
    a4c = convolve2d(a4,roi, mode = 'same') / a0c
    b2c = convolve2d(b2,roi, mode = 'same') / a0c
    b4c = convolve2d(b4,roi, mode = 'same') / a0c


    I2 = np.sqrt(a2c**2 + b2c**2)
    Phi = 0.5*np.arctan2(b2c,a2c)*180/np.pi
    Phi2 = (Phi>=0)*Phi + (Phi<=0)*(Phi+180)
    
    I4 = np.sqrt(a4c**2 + b4c**2)
    Phi = 0.25*np.arctan2(b4c,a4c)*180/np.pi
    Phi4 = (Phi>=0)*Phi + (Phi<=0)*(Phi+180)
    
    I4a = I4*np.sin(4*(Phi4 - Phi2)*(np.pi/180))    #Note - this is how I4a and I4s are defined in the MATLAB code,
    I4s = I4*np.cos(4*(Phi4 - Phi2)*(np.pi/180))    # but in the 2018 supplementary materials, the sin and cos are the otherway around
    
    return I2, Phi2, I4, Phi4, I4a, I4s