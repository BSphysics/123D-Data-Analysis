# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 10:55:35 2022

@author: Ben
"""
from pathlib import Path
import os
from rectangleFinder import rectangle_finder
from pSHG_histograms_ANA import pSHGhistogramsANA
import matplotlib.pyplot as plt

def anaROI(data_path, phi2, I2, dsImgs, imgs):
    
    sPath = data_path + '\\ROIs'
    save_path = Path(sPath)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
        print("Directory " , save_path ,  " Created ")
    else:    
        print("Directory " , save_path ,  " already exists")
    
 
    ROISize = 50
    downscaledROISize = 12
    print(data_path)
    coords, ROIfilename = rectangle_finder(ROISize, dsImgs, imgs)
    plt.savefig(data_path + '//' + ROIfilename[:-4] +' ROIs.png', bbox_inches='tight',pad_inches=0)
    
    # print(coords)
    
    for idx in range(len(coords)):  
        n = pSHGhistogramsANA(phi2[coords[idx][1]:coords[idx][1]+downscaledROISize,coords[idx][0]:coords[idx][0]+downscaledROISize],'Phi2', sPath, str(idx))
    
    print('\n\n')
    
    for idx in range(len(coords)):  
        n = pSHGhistogramsANA(I2[coords[idx][1]:coords[idx][1]+downscaledROISize,coords[idx][0]:coords[idx][0]+downscaledROISize],'I2', sPath, str(idx))
        
    for idx in range(len(coords)):  
        plt.figure()
        plt.imshow(dsImgs[coords[idx][1]:coords[idx][1]+downscaledROISize,coords[idx][0]:coords[idx][0]+downscaledROISize], cmap = 'gray', vmin = 0, vmax = 1)
        
        plt.savefig(data_path + '//' + ROIfilename[:-4] +' ROI' + str(idx) + '.png', bbox_inches='tight',pad_inches=0)
        