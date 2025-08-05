# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 09:09:16 2022

@author: Ben
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def pSHGhistogramsANA(data,name,data_path,ROI_num):
    
     if name == 'Phi2':
         
        data = data[data != np.inf]
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
        
               
        n, bins, patches = ax1.hist(np.ravel(data), bins = range(-90,90,2), edgecolor = 'gray', facecolor='yellow', alpha=0.5)     
        bin_centers = bins[:-1] + np.diff(bins) / 2
        ax1.set_title(name + ' histogram')
        ax1.set_xticks(range(-90,76,15))        
        
        circShift = int(len(n)/2 - (np.argmax(n)))
        nr = np.roll(n, circShift)
        w = (bins[1]-bins[0])*0.9
        
        xTickLabels = np.linspace(-90,80,18, dtype=str)
        for idx in range(0,len(xTickLabels)):
            q = xTickLabels[idx]
            xTickLabels[idx] = q[:-2]
        
        
        ax2.bar(bin_centers, nr, width = w, align = 'center', edgecolor = 'gray', facecolor = 'green', alpha=0.5)    
        ax2.set_xticks(range(-90,81,10))
        ax2.set_xticklabels(np.roll(xTickLabels, int(np.round(circShift/4.5))+1))
        ax2.set_title('Shifted ' + ' histogram')
        
        def gaussian_1D(x, a, b, c):
            return a * np.exp(-(x-b)**2/(2*c**2))
         
        popt, pcov = curve_fit(gaussian_1D, bin_centers, nr, bounds=([0, -np.inf,-np.inf], [np.inf,np.inf,np.inf]))
        
        xTickLabels = np.linspace(-90,80,18, dtype=str)
        for idx in range(0,len(xTickLabels)):
            q = xTickLabels[idx]
            xTickLabels[idx] = q[:-2]
            
        
        xDataforFit = np.linspace(bins[0], bins[-1], 10000)
        ax3.plot(xDataforFit, gaussian_1D(xDataforFit, *popt),'--r') 
        ax3.bar(bin_centers, nr, width = w, align = 'center', edgecolor = 'gray',facecolor = 'green', alpha= 0.1)    
        ax3.set_xticks(range(-90,81,10))
        ax3.set_xticklabels(np.roll(xTickLabels, int(np.round(circShift/4.5))+1))
        ax3.set_title('FWHM = ' + str(np.round(abs(popt[2]*2.355),1)) + ' degrees')
        print('Phi2 mean = ' + str(np.round(np.mean(np.ravel(data)),1)) + ' +/- ' + str(np.round(np.std(np.ravel(data)),1)) + ' degrees')
        
        plt.savefig(data_path + '\\' + name + '  ROI ' + ROI_num + ' histogram.png', bbox_inches='tight',pad_inches=0)
        return n 
        
     elif name == 'I2':
         
         data = data[data != np.inf]
         
         fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 4))
        
         n2, bins, patches = ax1.hist(np.ravel(data), bins = np.arange(0,1.5,0.02), edgecolor = 'black', facecolor='blue', alpha=0.5)
         bin_centers = bins[:-1] + np.diff(bins) / 2
         ax1.set_title(name + ' histogram')
         w = (bins[1]-bins[0])*0.9
        
        
             
         ax2.bar(bin_centers, n2, width = w, align = 'center', edgecolor = 'black', facecolor = 'blue', alpha=0.1)    
        
         def gaussian_1D(x, a, b, c):
             return a * np.exp(-(x-b)**2/(2*c**2))
         
         popt, pcov = curve_fit(gaussian_1D, bin_centers, n2)
        
         xDataforFit = np.linspace(bins[0], bins[-1], 10000)
         ax2.plot(xDataforFit, gaussian_1D(xDataforFit, *popt),'--r') 
         ax2.set_title('FWHM = ' + str(np.round(abs(popt[2]*2.355),3)) +' Centre = '+ str(np.round((popt[1]),3)) )
         print('I2 mean = ' + str(np.round(np.mean(np.ravel(data)),3)) + ' +/- ' + str(np.round(np.std(np.ravel(data)),3)))
         plt.savefig(data_path + '\\' + name + '  ROI ' + ROI_num + ' histogram.png', bbox_inches='tight',pad_inches=0)
         
     else:
         print('Name needs to be either Phi2 or I2')
        
               
         
