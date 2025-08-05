# -*- coding: utf-8 -*-
"""
Created on Thu May  5 22:46:40 2022

@author: Ben
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import LassoSelector
from pathlib import Path as Path1
from matplotlib.path import Path as Path2
import os

def lassoROI(im, phi2, I2, mask, time2wait, data_path):
    
    sPath = data_path + '\\lasso ROI results'
      
    save_path = Path1(sPath)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    numberOfExistingFolders = len(next(os.walk(sPath))[1]) # Neat way to count the number of oflders already in a directory
    lassoSavePath = sPath + '\\ROI' + str(numberOfExistingFolders + 1)
    os.mkdir(Path1(lassoSavePath)) # Make a new folder each time this part of the code is run => prevents overwriting accidentally
    
    fig, ax = plt.subplots()
    ax.imshow(im*mask,cmap = 'gray', vmin = 0, vmax = np.max(im)*0.5)
    
    def onSelect(verts):
        global ind, p
        p = Path2(verts)
        grid = [(i,j) for j in range(int(im.shape[0])) for i in range(int(im.shape[0]))]
        ii = np.nonzero([p.contains_point(xy) for xy in grid])
        ind = [grid[i] for i in ii[0]]

    def onPress(event):
        print('\n  \n')
        
    def onRelease(event):
        print( '\n \n')    
    lineProps = {'color': 'red', 'linewidth': 2, 'alpha': 0.8}

    lsso = LassoSelector( ax = ax , onselect = onSelect, lineprops = lineProps, button = 1)

    fig.canvas.mpl_connect('button_press_event' , onPress)
    fig.canvas.mpl_connect('button_release_event' , onRelease)

    plt.show()
    plt.title(str(time2wait) + ' seconds to draw an ROI')
    plt.pause(time2wait)

    width, height = len(im), len(im)
    x, y = np.mgrid[:height, :width]
    coors=np.hstack((x.reshape(-1, 1), y.reshape(-1,1))) # coors.shape is (4000000,2)

    m = p.contains_points(coors)
    n = np.reshape(m, (width, height))
    n = n.transpose(1,0)*1.0
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
    ax1.imshow(n*im, cmap = 'gray')
    ax1.axis('off')
    ax1.set_title('Intensity ROI')
    
    w = np.ma.masked_where(n*phi2*mask == 0, n*phi2)
    
        
    mycm = plt.cm.get_cmap("hsv").copy()
    mycm.set_bad(color='black') 
    
    ax2.imshow(w, cmap = mycm, vmin = 0, vmax = 180)
    ax2.axis('off')
    ax2.set_title('Phi2 ROI')
    
    
    p3 = ax3.imshow(n*I2*mask, vmin = 0, vmax = 0.3, cmap = 'hot')
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax3. axis('off')
    ax3.set_title('I2 ROI')
    plt.colorbar(p3, cax=cax)
    
    filename = 'ROI mask.png'
    plt.savefig(lassoSavePath + '\\' + filename)
    
    filename = 'ROI mask'
    np.save(lassoSavePath + '\\' + filename , mask)
    
    return n  , str(lassoSavePath)