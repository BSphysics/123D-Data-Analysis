# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 13:13:28 2022

@author: Ben
"""
import numpy as np
import matplotlib.pyplot as plt


def pSHGpolar(data, data_path):

    # plt.figure()
    plt.figure(figsize = (14,12)) 
    N = len(data)
    bottom = 0
    
    theta = np.linspace(-np.pi*0.5, 0.5 * np.pi, N, endpoint=True)
    radii = data
    width = (np.pi) / N
    
    ax = plt.subplot(111, polar=True, theta_offset= np.pi/2)
    bars = ax.bar(theta, radii, width=width, bottom=bottom)
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    ax.get_yaxis().set_visible(False)
    
    # Use custom colors and opacity
    for r, bar in zip(radii, bars):
        bar.set_facecolor(plt.cm.viridis(r/np.max(radii)))
        bar.set_alpha(0.8)
    plt.title('Phi2 polar plot')
    plt.show()
    
    plt.savefig(data_path + '\\' + 'Phi2 polar plot.png', bbox_inches='tight',pad_inches=0)
