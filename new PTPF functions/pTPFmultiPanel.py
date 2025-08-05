# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 22:33:24 2021

@author: Ben
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from binary_mask import binary_mask
from skimage.morphology import disk
from skimage.filters import median
import os

def pTPFmultiPanel(allSum, phi2, I2, I4, I4a, I4s, threshLow, threshHigh, data_path):
    
    mask = binary_mask( threshLow, threshHigh, allSum)   
    im = median(mask*phi2, disk(3))
    im = np.ma.masked_where(im == 0.0, im)
    
    mycm = plt.cm.get_cmap("hsv").copy()
    mycm.set_bad(color='black')    
    
    fig, (ax1, ax3, ax4, ax5, ax6) = plt.subplots(1, 5, figsize=(24, 8))
    
    p1 = ax1.imshow(allSum*mask, cmap = 'gray', vmin = 0, vmax = np.max(allSum))
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax1. axis('off')
    ax1.set_title('TPF intensity')
    plt.colorbar(p1, cax=cax, format = '%.0e')     
    
    p3 = ax3.imshow(im, cmap = mycm,  vmin = 0, vmax = 180)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax3. axis('off')
    ax3.set_title('Phi2')
    plt.colorbar(p3, cax=cax)
        
    im = median(mask*I2, disk(3))
    im = np.ma.masked_where(im == 0, im)
    
    mycm = plt.cm.get_cmap("hot").copy()
    mycm.set_bad(color='black') 
    
    p4 = ax4.imshow(im, cmap = mycm, vmin = 0, vmax = np.max(im))
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax4. axis('off')
    plt.colorbar(p4, cax=cax)
    ax4. set_title('I2')
    
    im = median(mask*I4a, disk(1))
    im = np.ma.masked_where(mask == 0, im)
    
    mycm = plt.cm.get_cmap("jet").copy()
    mycm.set_bad(color='black') 
    
    p5 = ax5.imshow(im, cmap = mycm, vmin = -0.3, vmax = 0.3)
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax5. axis('off')
    plt.colorbar(p5, cax=cax)
    ax5. set_title('I4a')
    
    im = median(mask*I4s, disk(1))
    im = np.ma.masked_where(mask == 0, im)  
    
    p6 = ax6.imshow(im, cmap = mycm, vmin = -0.3, vmax = 0.3)
    divider = make_axes_locatable(ax6)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax6. axis('off')
    plt.colorbar(p6, cax=cax)
    ax6. set_title('I4s')
    
    plt.savefig(data_path + '\\pTPF multipanel.png', dpi=200, bbox_inches='tight',pad_inches=0)
    
    if not os.path.exists(data_path + '/pTPF stack parameters'):
        os.mkdir(data_path + '/pTPF stack parameters')
    
    # np.save(data_path + '/pSHG stack parameters' + '/I2 mean.npy' , np.mean(I2*mask))
    # np.save(data_path + '/pSHG stack parameters' + '/I2 stdev.npy' , np.std(I2*mask))
    # np.save(data_path + '/pSHG stack parameters' + '/I4 mean.npy' , np.mean(I4*mask))
    # np.save(data_path + '/pSHG stack parameters' + '/I4 stdev.npy' , np.std(I4*mask))
    # np.save(data_path + '/pSHG stack parameters' + '/I4a mean.npy' , np.mean(I4a*mask))
    # np.save(data_path + '/pSHG stack parameters' + '/I4a stdev.npy' , np.std(I4a*mask))
    # np.save(data_path + '/pSHG stack parameters' + '/I4s mean.npy' , np.mean(I4s*mask))
    # np.save(data_path + '/pSHG stack parameters' + '/I4s stdev.npy' , np.std(I4s*mask))
    
    import xlsxwriter
    
    wb = xlsxwriter.Workbook(data_path + '\\pTPF stack parameters\\'  + os.path.split(data_path)[1] + 'pTPF stack parameters.xlsx')
    I2WS = wb.add_worksheet('pTPF stack parameters')
    I2WS.write(1,2,'Total TPF pixel counts')
    I2WS.write(1,3,'I2 mean')
    I2WS.write(1,4,'I2 stdev')
    I2WS.write(1,5,'I4 mean')
    I2WS.write(1,6,'I4 stdev')
    I2WS.write(1,7,'I4a mean')
    I2WS.write(1,8,'I4a stdev')
    I2WS.write(1,9,'I4s mean')
    I2WS.write(1,10,'I4s stdev')
    
    I2WS.write(2,2,np.sum(allSum))
    I2WS.write(2,3,np.nanmean(np.ma.masked_equal(I2*mask, 0)))
    I2WS.write(2,4,np.nanstd(np.ma.masked_equal(I2*mask, 0)))
    I2WS.write(2,5,np.nanmean(np.ma.masked_equal(I4*mask, 0)))
    I2WS.write(2,6,np.nanstd(np.ma.masked_equal(I4*mask, 0)))
    I2WS.write(2,7,np.nanmean(np.ma.masked_equal(I4a*mask, 0)))
    I2WS.write(2,8,np.nanstd(np.ma.masked_equal(I4a*mask, 0)))
    I2WS.write(2,9,np.nanmean(np.ma.masked_equal(I4s*mask, 0)))
    I2WS.write(2,10,np.nanstd(np.ma.masked_equal(I4s*mask, 0)))

    
    
    wb.close()
    
    return mask