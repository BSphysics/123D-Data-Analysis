# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 22:57:00 2018

@author: Ben
"""

import os
from skimage import io
from natsort import natsorted
import numpy as np


def tiff_loader(data_path):
    
    files = os.listdir(data_path)
    
    files_tiff = [i for i in files if i.endswith('.tif')]
    
    files_tiff = natsorted(files_tiff)
        
    imgs = []
    
    for idx in range(len(files_tiff)):
    
        filename = os.fsdecode(files_tiff[idx])    
        
        filename = r'/' + filename
        
        im = io.imread(data_path + filename)
        
        imgs.append(im)
        
    imgs = np.array(imgs, dtype='int64')
       
    
    return files_tiff, imgs
