# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 22:37:56 2022

@author: Ben
"""

from tkinter import filedialog
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
from matplotlib.patches import Rectangle
import os

def rectangle_finder(roiPix, dsImgs, imgs):

    data_path = r'C:\Users\Ben\OneDrive - University of Exeter\!Work\Work.2022\Data\123D\Ana\areas for pSHG analysis'
    f = filedialog.askopenfilename(initialdir=data_path , filetypes=[("Tiffs", ".tif")])
    im = io.imread(f)
    filename = os.path.split(f)[1]
    #%%
    
    print('size of ROI labelled tiff image = ' + str(im.shape))
    label_img = label(im[:,:,0])
    regions = regionprops(label_img)
    
    fig, (ax1,ax2,ax3 , ax4) = plt.subplots(1, 4, figsize = (10,6))
    
    ax1.imshow(np.sum(im,2))
    ax1.axis('off')
    ax2.imshow(np.sum(im,2))
    ax2.axis('off')
    ax3.imshow(np.sum(imgs,0), cmap = 'gray')
    ax3.axis('off')
    ax4.imshow(dsImgs, cmap = 'gray')
    ax4.axis('off')
    idx=0
    colours = ['r', 'b', 'y']
    coords = []
    for props in regions:
      
        y0, x0 = props.centroid
        minr, minc, maxr, maxc = props.bbox
        # bx = (minc, maxc, maxc, minc, minc)
        # by = (minr, minr, maxr, maxr, minr)
           
        ax2.add_patch(Rectangle([minc,minr], roiPix, roiPix, color = colours[idx], alpha = 0.4))
        ax3.add_patch(Rectangle([minc,minr], roiPix, roiPix, color = colours[idx], alpha = 0.4))
        
        coords.append([int(np.round(minc/4)),int(np.round(minr/4))])   
        ax4.add_patch(Rectangle([int(np.round(minc/4)),int(np.round(minr/4))], 12, 12, color = colours[idx], alpha = 0.4))
        idx = idx + 1
        
    # ax3.imshow(np.sum(imgs,0), cmap = 'gray')
    # ax3.axis('off')
    
    fig.tight_layout() 
    fig.suptitle(filename[:-4])   
    plt.show()
    plt.savefig(data_path + '\\' + filename[:-4] +'.png', bbox_inches='tight',pad_inches=0)
    
    return coords, filename
    


