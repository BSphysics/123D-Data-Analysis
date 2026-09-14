# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 21:49:31 2021
@author: Ben

 Code to analyse polarisation resolved second harmonic generation data set acquired in Physics 123D

    Added I4a and I4s to pshgProjector, but not sure the values are right. These need to be checked.
"""
import os
scriptDir = os.getcwd()
import sys
sys.path.append(os.path.join(scriptDir,"new pSHG functions" ))
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from tiff_loader import tiff_loader
import cv2
import matplotlib.patches as patches
plt.close('all')

#%% Select folder where pSHG data is saved and read all images into memory
from pSHG_GUI import pSHGGUI


defaultDir = scriptDir + '\\Test data'
initialDir = r'C:\Users\bs426\OneDrive - University of Exeter\!Work\Work.2024\Lab 2024\123D'
defaultThreshold = 1e2

[data_path, plotHistograms, pSHGFitViewer, plotPolarHistogram, arrowPlot,
 usepolygonalROI, useSlider, thresh, polygonal_roi_image,
 maskMode, sigmaPhi2Max, poolSigma, minFeature, closeRadius] = pSHGGUI(
    initialDir, defaultDir, defaultThreshold, 'transmission')    
#%%
[filenames, imgs] = [[],[]]
filenames, imgs = tiff_loader(data_path)

from imageMetaData import imageMetaData
[zoomFactor , frames, zStackStep, allMetaData] = imageMetaData(data_path, filenames)

# Create Transmission pSHG results folder inside data_path
results_folder = os.path.join(data_path, "Transmission pSHG results")

if not os.path.exists(results_folder):
    os.makedirs(results_folder)

# Save directory
folderName = results_folder

#%% 
[pshg, ptpf, sumSHG, sumTPF, I2, Phi2, pshg_raw] = [[],[],[],[],[],[],[]]

for idx in range(len(imgs)-2): #MAIN LOOP

    im = imgs[idx]
    imShape = im.shape
    if min(imShape)>2:
       im = np.transpose(im, axes=(2,0,1)) 
       
    if im.ndim>2:
        shg = im[2, :, :].astype(float)
        shg_raw = shg.copy()          
        shg[shg < 0] = 0              
   

        tpf = im[1,:,:]
        shg[shg<0] = 0
        tpf[tpf<0] = 0   
    else:
        shg = im
        shg[shg < 0] = 0
        tpf=0
    pshg.append(shg)
    pshg_raw.append(shg_raw)  
    ptpf.append(tpf)
    sumSHG.append(np.sum(shg))
    sumTPF.append(np.sum(tpf))
pshg = np.asarray(pshg)   
ptpf = np.asarray(ptpf)    
allSum = np.sum(pshg,0)
allSumTPF = np.sum(ptpf,0)

from pshgProjector import pshgProjector
delta = -30 #NOTE: This only affects Phi2, not I2 or any other parameters.
angles = (np.arange(0,len(pshg))*15 - delta) * np.pi/180
[I2, Phi2, I4, Phi4, I4a, I4s] = pshgProjector(data_path, pshg, np.flip(angles))

#%%
from pSHG_gof import pshg_gof
from pSHG_SNR_mask import snr_trust_mask
import copy
if maskMode == 'trust':

    res  = pshg_gof(np.asarray(pshg_raw), angles)
    SNR_mask = snr_trust_mask(res, sigmaPhi2Max, poolSigma,
                          min_feature=minFeature, close_radius=closeRadius)['mask']
    # --- bail out cleanly if too few pixels survive the SNR filter ---
    n_pass    = int(np.count_nonzero(SNR_mask))
    frac_pass = n_pass / SNR_mask.size
    if frac_pass < 1e-3:                       # < 0.1 % of pixels (covers the zero case)
        print(f"\nSNR filter left {n_pass} pixels "
              f"({frac_pass:.3%} of {SNR_mask.size}) — too sparse for a meaningful mask.")
        print("Loosen the filter and re-run: raise sigmaPhi2Max, drop minFeature/"
              "closeRadius, or lower the intensity threshold.")
        sys.exit()
    # -----------------------------------------------------------------

    gray_red = copy.copy(plt.cm.gray)      # grayscale, but...
    gray_red.set_bad('red')                # ...masked-out pixels drawn bright red
    vmin, vmax = np.percentile(allSum, [1, 99.5])   # shared scaling for both panels

    fig, ax = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    ax[0].imshow(allSum, cmap='gray', vmin=vmin, vmax=vmax)
    ax[0].set_title('Summed SHG intensity')
    ax[1].imshow(np.ma.masked_where(~SNR_mask, allSum), cmap=gray_red, vmin=vmin, vmax=vmax)
    ax[1].set_title('SNR-masked (rejected pixels = red)')
    for a in ax: a.axis('off')
    plt.show()
    plt.pause(0.1)
    while True:
        if plt.waitforbuttonpress(): break  
    
    plt.close(fig = 1)

else:
    SNR_mask = None
    from slider_thresh import sliderThresh
    plt.ion()
    slide = sliderThresh(allSum)
    plt.show()
    plt.pause(0.1)
    while True:
        if plt.waitforbuttonpress(): break  # Exit loop if user presses a key.
    print('\n Threshold = ' + str(np.round(slide.val)))
    plt.close(fig = 1)
    thresh = slide.val


#%%
def micronsPerPixel(zoomFactor):
    mpp = 1.9634*zoomFactor**(-0.987)
    return mpp
scaleBarinMicrons = float(50)
mpp = micronsPerPixel(zoomFactor)
scaleBarLength = scaleBarinMicrons / mpp
scaleBarWidth = 15


#%%
threshHigh=1e9
from pSHGmultiPanel import pSHGmultiPanel
mask = pSHGmultiPanel(allSum, Phi2, I2, I4, I4a, I4s, thresh, threshHigh, folderName, SNR_mask, maskMode)

#%%
scaleBar = patches.Rectangle((400, 475), scaleBarLength, scaleBarWidth, linewidth = 1, edgecolor = 'm', facecolor = 'w')
fig, ax = plt.subplots(figsize = (12,10))
im = (allSum**mask)/np.max(allSum) 
plt.imshow(np.power(im, .7), cmap = 'gray')#, vmin = 0, vmax = np.max(allSum)/1.0)
plt.axis('off') 
ax.add_patch(scaleBar) 
plt.title(os.path.basename(data_path) + ',  scale bar = ' + str(scaleBarinMicrons) + r' $\mu$m', fontsize = 6, loc="right")
plt.savefig(folderName + '\\' + 'all sum SHG.png', dpi = 400, bbox_inches='tight',pad_inches=0)
plt.close('all')

#%% Merge image with scale bar

allSum = np.sum(pshg,0)
allSumTPF = np.sum(ptpf,0)

color_image = np.zeros((allSum.shape[0], allSum.shape[1], 3), dtype=float)
color_image[:, :, 1] = allSum/(np.mean(allSum)*2.5)  # Green channel
color_image[:, :, 0] = allSumTPF/(np.mean(allSumTPF)*4)  # Red channel 
color_image = np.clip(color_image, 0, 1)   
fig, ax1 = plt.subplots()
plt.imshow(color_image[:,:,:])
plt.axis('off')
scaleBar1 = patches.Rectangle((400, 475), scaleBarLength, scaleBarWidth, linewidth = 1, edgecolor = 'm', facecolor = 'w')
ax1.add_patch(scaleBar1)
plt.title(os.path.basename(data_path) + ',  scale bar = ' + str(scaleBarinMicrons) + r' $\mu$m', fontsize = 6, loc="right")
plt.savefig(folderName + '\\' + 'merge_Image.png', dpi = 400, bbox_inches='tight',pad_inches=0)
plt.close('all')


#%%
if arrowPlot == True:
    from pSHG_arrows import pSHGArrows
    arrowColourMax = 0.5
    pSHGArrows(I2, Phi2, allSum, mask, folderName, arrowColourMax)
#%% View individual pixel fits
if pSHGFitViewer == True:
    from pSHGfitviewer import pSHGFitViewer
    points = pSHGFitViewer(pshg, I2, Phi2, mask, angles, folderName)
    
#%% Plot pSHG histograms

from pSHG_histograms_NEW import pSHGhistogramsNEW

if plotHistograms == True:
     
  
    [bins , binCounts, data] = pSHGhistogramsNEW(Phi2, mask,'Phi2', folderName)
    [bins , binCounts, data] = pSHGhistogramsNEW(I2, mask,'I2', folderName)
    [bins , binCounts, data] = pSHGhistogramsNEW(I4a, mask,'I4a', folderName)
    [bins , binCounts, data] = pSHGhistogramsNEW(I4s, mask,'I4s', folderName)

#%%
if plotPolarHistogram == True: 
    from pSHGpolar import pSHGpolar
    pSHGpolar(binCounts, folderName)         # Plot plot of the Phi2 histogram

#%%
plt.close('all')

from pSHG_histograms_NEW import pSHGhistogramsNEW
if usepolygonalROI:
   
    from polygonal_roi import polygonalROI
    
    if polygonal_roi_image == 'Phi2':
        roi_image = Phi2
    elif polygonal_roi_image == 'I2':
        roi_image = I2
    elif polygonal_roi_image == 'All SHG':
        roi_image = allSum
        
    [ROI , polygonalSavePath] = polygonalROI(roi_image, Phi2, I2, mask, folderName , plot_mode = polygonal_roi_image) #ROI is a mask than is 0 outside, and 1 inside the hand draw region
    [bins , binCounts , roi_data] = pSHGhistogramsNEW(Phi2*ROI, mask, 'Phi2', polygonalSavePath)
    [bins , binCounts , roi_data] = pSHGhistogramsNEW(I2*ROI, mask, 'I2', polygonalSavePath)
    [bins , binCounts , roi_data] = pSHGhistogramsNEW(I4a*ROI, mask, 'I4a', polygonalSavePath)
    [bins , binCounts , roi_data] = pSHGhistogramsNEW(I4s*ROI, mask, 'I4s', polygonalSavePath)
   
plt.close('all')
path = os.path.realpath(folderName)
os.startfile(path)   

#%%
filename = 'I2_masked'
fullFileName = os.path.join(folderName, filename)
np.save(fullFileName, I2*mask)

filename = 'Phi2_masked'
fullFileName = os.path.join(folderName, filename)
np.save(fullFileName, Phi2*mask)

filename = 'I4a_masked'
fullFileName = os.path.join(folderName, filename)
np.save(fullFileName, I4a*mask)

filename = 'I4s_masked'
fullFileName = os.path.join(folderName, filename)
np.save(fullFileName, I4s*mask)

#%%

from densityPlot import density_plot
filename = 'I4s_vs_I4a.png'
fullFileName = os.path.join(folderName, filename)
density_plot(I4a, I4s, I2=I2, draw_arc=True, out_path=fullFileName);
plt.close('all')




