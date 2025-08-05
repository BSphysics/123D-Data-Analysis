# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 09:09:16 2022

@author: Ben
"""

import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import xlsxwriter

def pSHGhistograms(data,mask,name,data_path):
    
     mask = mask.astype('float')
     mask[mask==0]=np.inf
     data = data*mask
     
     data = data[data != np.inf]
     data = data[data != -np.inf]
     data[-1] = 0
     
         
     if name == 'Phi2':  
         
        print('\n Number of pixels in histogram = ' + str(len(data)) + '\n')
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        
        n, bins, patches = ax1.hist(np.ravel(data), bins = range(0,180,2), edgecolor = 'black', facecolor='green', alpha=0.5)     
        bin_centers = bins[:-1] + np.diff(bins) / 2
        ax1.set_xticks(range(0,165,15))  
        ax1.set_title(name + ' histogram')
   
        n, bins, patches = ax2.hist(np.ravel(data), bins = range(0,180,2), edgecolor = 'black', facecolor='green', alpha=0.1)     
        bin_centers = bins[:-1] + np.diff(bins) / 2
        ax2.set_xticks(range(0,165,15))  
        
                
        def gaussian_1D(x, a, b, c):
            return a * np.exp(-(x-b-bins[np.argmax(n)])**2/(2*c**2))
        try: 
            pop1, pcov = curve_fit(gaussian_1D, bin_centers, n)
           
        except RuntimeError:
            print('Error - curve_fit failed :(')
            pop1 = 'Empty'

           
        if pop1.all() == 'Empty':
           print('Phi2 histogram fit failed')
        else:
            
            xDataforFit = np.linspace(bins[0], bins[-1], 1000)
            ax2.plot(xDataforFit, gaussian_1D(xDataforFit, *pop1),'--r') 
            ax2.set_title('FWHM = ' + str(np.round(abs(pop1[2]*2.355),1)) + ' degs , Centre = '+ str(bins[np.argmax(n)] + np.round((pop1[1]),1)) )
                    

        print('\n Phi2 mean = ' + str(np.round(np.mean(np.ravel(data)),1)) + ' +/- ' + str(np.round(np.std(np.ravel(data)),1)) + ' degrees')
        
        plt.savefig(data_path + '\\' + name + ' histogram.png', bbox_inches='tight',pad_inches=0)
                        
        wb = xlsxwriter.Workbook(data_path + '\\' + name + ' histogramData.xlsx')
        phi2WS = wb.add_worksheet('Phi2')
        phi2WS.write(1,4,'Histogram bins')
        phi2WS.write(1,5,'Histogram counts')
        phi2WS.write(1,7,'Phi2 pixel values')
        phi2WS.write_column(2,4,bins)
        phi2WS.write_column(2,5,n)
        phi2WS.write_column(2,7,data)
        wb.close()
        
     elif name == 'I2':
               
         fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        
         n, bins, patches = ax1.hist(np.ravel(data), bins = np.arange(0,1.5,0.02), edgecolor = 'black', facecolor='blue', alpha=0.5)
         bin_centers = bins[:-1] + np.diff(bins) / 2
         ax1.set_title(name + ' histogram')
         w = (bins[1]-bins[0])*0.9
        
        
             
         ax2.bar(bin_centers, n, width = w, align = 'center', edgecolor = 'black', facecolor = 'blue', alpha=0.1)    
        
         def gaussian_1D(x, a, b, c, d):
             return a * np.exp(-(x-b)**2/(2*c**2)) + d
          
         try: 
            pop2, pcov = curve_fit(gaussian_1D, bin_centers, n)
        
         except RuntimeError:
            print("Error - curve_fit failed :(")
            pop2 = 'Empty'
            
         if pop2.all() == 'Empty':
            print('I2 histogram fit failed')
         else:
        
             xDataforFit = np.linspace(bins[0], bins[-1], 10000)
             ax2.plot(xDataforFit, gaussian_1D(xDataforFit, *pop2),'--r') 
             ax2.set_title('FWHM = ' + str(np.round(abs(pop2[2]*2.355),3)) +' Centre = '+ str(np.round((pop2[1]),3)) )
             print('\n I2 mean = ' + str(np.round(np.mean(np.ravel(data)),3)) + ' +/- ' + str(np.round(np.std(np.ravel(data)),3)))
         plt.savefig(data_path + '\\' + name + ' histogram.png', bbox_inches='tight',pad_inches=0)
         
         
         wb = xlsxwriter.Workbook(data_path + '\\' + name + ' histogramData.xlsx')
         I2WS = wb.add_worksheet('I2')
         I2WS.write(1,4,'Histogram bins')
         I2WS.write(1,5,'Histogram counts')
         I2WS.write(1,7,'I2 pixel values')
         I2WS.write_column(2,4,bins)
         I2WS.write_column(2,5,n)
         I2WS.write_column(2,7,data)
         wb.close()
         
     elif name == 'I4a':
              
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
       
        n, bins, patches = ax1.hist(np.ravel(data), bins = np.arange(-0.5,0.5,0.01), edgecolor = 'black', facecolor='red', alpha=0.5)
        bin_centers = bins[:-1] + np.diff(bins) / 2
        ax1.set_title(name + ' histogram')
        w = (bins[1]-bins[0])*0.9
       
       
            
        ax2.bar(bin_centers, n, width = w, align = 'center', edgecolor = 'black', facecolor = 'red', alpha=0.1)    
       
        def gaussian_1D(x, a, b, c, d):
            return a * np.exp(-(x-b)**2/(2*c**2)) + d
         
        try: 
           pop2, pcov = curve_fit(gaussian_1D, bin_centers, n)
       
        except RuntimeError:
           print("Error - curve_fit failed :(")
           pop2 = 'Empty'
           
        if pop2.all() == 'Empty':
           print('I4a histogram fit failed')
        else:
       
            xDataforFit = np.linspace(bins[0], bins[-1], 10000)
            ax2.plot(xDataforFit, gaussian_1D(xDataforFit, *pop2),'--r') 
            ax2.set_title('FWHM = ' + str(np.round(abs(pop2[2]*2.355),3)) +' Centre = '+ str(np.round((pop2[1]),3)) )
            print('\n I4a mean = ' + str(np.round(np.mean(np.ravel(data)),3)) + ' +/- ' + str(np.round(np.std(np.ravel(data)),3)))
        plt.savefig(data_path + '\\' + name + ' histogram.png', bbox_inches='tight',pad_inches=0)
        
        
        wb = xlsxwriter.Workbook(data_path + '\\' + name + ' histogramData.xlsx')
        I2WS = wb.add_worksheet('I4a')
        I2WS.write(1,4,'Histogram bins')
        I2WS.write(1,5,'Histogram counts')
        I2WS.write(1,7,'I4a pixel values')
        I2WS.write_column(2,4,bins)
        I2WS.write_column(2,5,n)
        I2WS.write_column(2,7,data)
        wb.close()
         
     elif name == 'I4s':
                 
           fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
          
           n, bins, patches = ax1.hist(np.ravel(data), bins = np.arange(-0.5,0.5,0.01), edgecolor = 'black', facecolor='green', alpha=0.5)
           bin_centers = bins[:-1] + np.diff(bins) / 2
           ax1.set_title(name + ' histogram')
           w = (bins[1]-bins[0])*0.9
          
          
               
           ax2.bar(bin_centers, n, width = w, align = 'center', edgecolor = 'black', facecolor = 'green', alpha=0.1)    
          
           def gaussian_1D(x, a, b, c, d):
               return a * np.exp(-(x-b)**2/(2*c**2)) + d
            
           try: 
              pop2, pcov = curve_fit(gaussian_1D, bin_centers, n)
          
           except RuntimeError:
              print("Error - curve_fit failed :(")
              pop2 = 'Empty'
              
           if pop2.all() == 'Empty':
              print('I4s Histogram fit failed')
           else:
          
               xDataforFit = np.linspace(bins[0], bins[-1], 10000)
               ax2.plot(xDataforFit, gaussian_1D(xDataforFit, *pop2),'--r') 
               ax2.set_title('FWHM = ' + str(np.round(abs(pop2[2]*2.355),3)) +' Centre = '+ str(np.round((pop2[1]),3)) )
               print('\n I4s mean = ' + str(np.round(np.mean(np.ravel(data)),3)) + ' +/- ' + str(np.round(np.std(np.ravel(data)),3)))
           plt.savefig(data_path + '\\' + name + ' histogram.png', bbox_inches='tight',pad_inches=0)
           
           
           wb = xlsxwriter.Workbook(data_path + '\\' + name + ' histogramData.xlsx')
           I2WS = wb.add_worksheet('I4s')
           I2WS.write(1,4,'Histogram bins')
           I2WS.write(1,5,'Histogram counts')
           I2WS.write(1,7,'I4s pixel values')
           I2WS.write_column(2,4,bins)
           I2WS.write_column(2,5,n)
           I2WS.write_column(2,7,data)
           wb.close()
     else:
         print('\n Name needs to be either Phi2, I2, I4a or I4s')
         bin_centers, n = [[],[]]
         
     # data = []    
     return bin_centers, n, data
        
               
         
