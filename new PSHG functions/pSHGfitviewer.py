# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 22:40:09 2021

@author: Ben
"""
import sys
sys.path.append(r"C:\Users\Ben\Desktop\Ben\Work\Software\fns")

import numpy as np
import matplotlib.pyplot as plt
from pshg_I_Alpha import I_alpha
from skimage.morphology import disk
from skimage.filters import median

def pSHGFitViewer(pshg, I2, Phi2, mask, angles, data_path):
    
   mycm = plt.cm.hsv.copy()
   mycm.set_bad(color='black') 
   
   im = median(mask*Phi2, disk(3))
   im = np.ma.masked_where(im == 0, im)

   roi=2
   plt.figure(figsize = (14,12)) 
   plt.imshow(im, cmap = mycm)
   plt.clim(vmin=0, vmax=180)
   plt.colorbar()
   plt.axis('off')
   plt.title('Phi2 map - Click on pixel to see fit')
   props = dict(boxstyle='round', facecolor='white', alpha=1.0)
   ax = plt.gca()
   plt.savefig(data_path + '\\' + 'fit viewer points.png', bbox_inches='tight',pad_inches=0)
   
   anglesFit = np.arange(0,180,2)*np.pi/180
   angles = (np.arange(0,len(pshg))*15) *np.pi/180
   # angle

   pts = np.round(plt.ginput(6, 0, show_clicks = True))

   idx = 0
   for idx in range(0,len(pts)):
       pt0 = int(pts[idx][0])
       pt1 = int(pts[idx][1])
       ax.text(pt0, pt1, str(idx+1), fontsize=10, bbox=props)
   plt.savefig(data_path + '\\' + 'fit viewer points.png', bbox_inches='tight',pad_inches=0)

   fig = plt.figure(figsize = (16,8))   
   for idx in range(0,len(pts)):
       pt0 = int(pts[idx][0])
       pt1 = int(pts[idx][1])    
       pshgTrace = np.sum(pshg[:,pt1-roi:pt1+roi,pt0-roi:pt0+roi], axis=(1,2))
       
       pk2pk = np.round(np.max(pshgTrace) - np.min(pshgTrace),1)
       
       A0trace = np.mean(pshgTrace,0)
       A2trace = 2*np.mean(pshgTrace*np.cos(2*angles),0)
       A4trace = 2*np.mean(pshgTrace*np.cos(4*angles),0)
       B2trace = 2*np.mean(pshgTrace*np.sin(2*angles),0)
       B4trace = 2*np.mean(pshgTrace*np.sin(4*angles),0)

       fit = A0trace + A2trace*np.cos(2*anglesFit) + B2trace*np.sin(2*anglesFit) + A4trace*np.cos(4*anglesFit) + B4trace*np.sin(4*anglesFit)
       
       ax = fig.add_subplot(2, 3, idx+1)
       ax.plot(angles *180/np.pi , pshgTrace, 'ro')
       ax.plot(anglesFit*180/np.pi, fit, 'b--')
       ax.set_title('Pt ' + str(idx + 1) + ', pk-to-pk = ' + str(pk2pk))
       
   plt.savefig(data_path + '\\' + 'fit viewer.png', bbox_inches='tight',pad_inches=0)
   
   return pts