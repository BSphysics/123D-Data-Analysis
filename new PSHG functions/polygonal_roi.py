# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 12:13:20 2026

@author: bs426
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path as Path1
from matplotlib.path import Path as Path2
from skimage.filters import median
from skimage.morphology import disk
import os
import matplotlib
matplotlib.use('Qt5Agg') 

def polygonalROI(im, phi2, I2, mask, data_path , plot_mode = None):
    
    sPath = data_path + '\\polygon ROI results'
      
    save_path = Path1(sPath)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    numberOfExistingFolders = len(next(os.walk(sPath))[1])
    polySavePath = sPath + '\\ROI' + str(numberOfExistingFolders + 1)
    os.mkdir(Path1(polySavePath))
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    if plot_mode == 'Phi2':
        mycm = plt.cm.hsv.copy()
        mycm.set_bad(color='black') 
        
        im = median(mask*im, disk(3))
        im = np.ma.masked_where(im == 0, im)

        ax.imshow(im, cmap = mycm , vmin=0, vmax=180) 
    
    if plot_mode == 'I2':
        im = median(mask*I2, disk(3))
        im = np.ma.masked_where(im == 0, im)
        mycm = plt.cm.get_cmap("hot").copy()
        mycm.set_bad(color='black') 
        ax.imshow(im,cmap = mycm, vmin = 0, vmax = np.mean(I2)+np.std(I2))
    
    if plot_mode == 'All SHG':
        im = (im**mask)/np.max(im) 
        ax.imshow(np.power(im, .7), cmap = 'gray')
        
        # ax.imshow(im * mask, cmap='gray', vmin = 0, vmax = np.max(im) * 0.75)
    plt.title('left-click to add points, right-click to close polygon')
    plt.axis('off')
    
    clicked_points = []
    lines = []
    closing_line = [None]
    poly_path = [None]

    def onClick(event):
        if event.inaxes != ax:
            return
        
        # Right-click closes the polygon
        if event.button == 3:
            if len(clicked_points) >= 3:
                # Draw closing line back to first point
                xs = [clicked_points[-1][0], clicked_points[0][0]]
                ys = [clicked_points[-1][1], clicked_points[0][1]]
                closing_line[0], = ax.plot(xs, ys, 'y--', linewidth=2, alpha=0.8)
                
                verts = clicked_points + [clicked_points[0]]  # close the path
                poly_path[0] = Path2(verts)
                fig.canvas.draw()
                print(f'Polygon closed with {len(clicked_points)} vertices.')
            else:
                print('Need at least 3 points to close a polygon.')
            return
        
        # Left-click adds a point
        if event.button == 1:
            x, y = event.xdata, event.ydata
            clicked_points.append((x, y))
            # Draw highly visible marker
            ax.plot(x, y, 'o', markersize=10, color='black', zorder=3)
            ax.plot(x, y, 'o', markersize=6, color='yellow', zorder=4)
            
            # Draw line from previous point to this one
            if len(clicked_points) > 1:
                xs = [clicked_points[-2][0], clicked_points[-1][0]]
                ys = [clicked_points[-2][1], clicked_points[-1][1]]
                line, = ax.plot(xs, ys, 'r--', linewidth=2, alpha=0.8)
                lines.append(line)
            
            fig.canvas.draw()

    fig.canvas.mpl_connect('button_press_event', onClick)
    plt.show(block=False)
    while poly_path[0] is None:
        plt.pause(0.1)
    plt.close()
    
    # Fall back to a full-image mask if user didn't close the polygon
    if poly_path[0] is None:
        print('Warning: polygon was not closed — returning empty mask.')
        return np.zeros_like(im, dtype=float), str(polySavePath)
    
    p = poly_path[0]
    width, height = len(im), len(im)
    x, y = np.mgrid[:height, :width]
    coors = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    m = p.contains_points(coors)
    n = np.reshape(m, (width, height))
    n = n.transpose(1, 0) * 1.0

    # ── Results figure (identical layout to lassoROI) ──────────────────────
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
    ax1.imshow(n * im, cmap='gray')
    ax1.axis('off')
    ax1.set_title('Intensity ROI')

    w = np.ma.masked_where(n * phi2 * mask == 0, n * phi2)
    mycm = plt.cm.get_cmap("hsv").copy()
    mycm.set_bad(color='black')
    ax2.imshow(w, cmap=mycm, vmin=0, vmax=180)
    ax2.axis('off')
    ax2.set_title('Phi2 ROI')

    p3 = ax3.imshow(n * I2 * mask, vmin=0, vmax=0.3, cmap='hot')
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax3.axis('off')
    ax3.set_title('I2 ROI')
    plt.colorbar(p3, cax=cax)

    filename = 'ROI mask.png'
    plt.savefig(polySavePath + '\\' + filename)

    filename = 'ROI mask'
    np.save(polySavePath + '\\' + filename, n)

    return n, str(polySavePath)
