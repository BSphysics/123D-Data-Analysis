# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 22:55:05 2021

@author: Ben
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl

def pTPFArrows(I2, Phi2, allSum, mask, data_path, arrowColourMax):

    cmap = plt.cm.jet
    # cNorm  = colors.Normalize(vmin = 0, vmax = np.max(I2*mask))
    cNorm  = colors.Normalize(vmin = 0, vmax = arrowColourMax)
    scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=cmap)
    
    fig = plt.figure(figsize = (12,10)) 
    ax  = fig.add_axes([0.1, 0.1, 0.80, 0.80]) # [left, bottom, width, height]
    axc = fig.add_axes([0.85, 0.10, 0.05, 0.8])
    
    ax.imshow(allSum*mask, cmap = 'gray', vmin = 0, vmax =np.max(allSum))
    
    Phi2 = np.where(mask == 0 , np.inf, Phi2)
    n = 6 # Add every nth arrow to intensity image
    arrowLength = 4
    angularOffset = -20
            
    for rows in range(0, int(Phi2.shape[1] / n) - 0):
        # print('row = '+ str(rows) + '\n')
        for cols in range(0, int(Phi2.shape[0] / n) - 0):
            pt0 = rows*n + int(n/2)
            pt1 = cols*n + int(n/2)
            if Phi2[pt1,pt0] == np.inf:
    
                arrowX = 0
                arrowY = 0
                arrowAlpha = 0
                colorVal= [0,0,0]
            else:
    
                arrowX = arrowLength*np.cos((Phi2[pt1,pt0]+ angularOffset) / 180*np.pi)
                arrowY = arrowLength*np.sin((Phi2[pt1,pt0]+ angularOffset) / 180*np.pi)
    
                arrowAlpha = 0.5       
                colorVal = scalarMap.to_rgba(I2[pt1, pt0])
            
            ax.arrow(pt0,pt1,arrowX, arrowY, color=colorVal, linewidth = 1, alpha = arrowAlpha, head_width=1, head_length = 1)
    
            
    ax.axis('off')
    mpl.colorbar.ColorbarBase(axc, cmap=cmap,norm=cNorm,orientation='vertical')
    
    plt.show()
    fig.suptitle('pTPF - Arrow color showing local I2 value')
    
    plt.savefig(data_path + '\\' + 'pTPF Arrows.png', dpi = 400, bbox_inches='tight',pad_inches=0)