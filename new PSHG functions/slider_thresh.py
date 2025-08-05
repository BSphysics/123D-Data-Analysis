# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 22:44:38 2020

@author: Ben
"""

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def sliderThresh(im):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.25, bottom=0.25)
    ax.set_title('Press any key to exit')
    min0 = 0
    max0 = np.max(im)
  

    im1 = ax.imshow(im, cmap = 'gray', vmin = min0, vmax= max0*0.9)
    fig.colorbar(im1)   
    axcolor = 'lightgoldenrodyellow'
    
    thresh  = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor = axcolor)    
    sthresh = Slider(thresh, 'Threshold', min0, max0*0.9, valinit = 1e2)
    
    def update(val):
        a = np.copy(im)
        t = sthresh.val
        a[a<t] = 0
        im1.set_data(a)
        fig.canvas.draw()
    sthresh.on_changed(update)
    
    plt.show()
    
    return sthresh