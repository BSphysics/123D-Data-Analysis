# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 15:14:09 2018

@author: Ben
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector

def region_selector(im, mn, mx, pltitle, sPath, filename):

    fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize=(10, 8))
    plt.imshow(im, cmap = 'gray',vmin = mn, vmax = mx)
    plt.title(pltitle + '  (in < 5 seconds!)')
    plt.pause(1.0)
    def line_select_callback(eclick, erelease):
        'eclick and erelease are the press and release events'
#        global x1, y1, x2, y2
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
#        print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
#        print(" The button you used were: %s %s" % (eclick.button, erelease.button))
    
    def toggle_selector(event):
        print(' Key pressed.')
        if event.key in ['Q', 'q'] and toggle_selector.RS.active:
            print(' RectangleSelector deactivated.')
            toggle_selector.RS.set_active(False)
        if event.key in ['A', 'a'] and not toggle_selector.RS.active:
            print(' RectangleSelector activated.')
            toggle_selector.RS.set_active(True)
            
    # drawtype is 'box' or 'line' or 'none'
    toggle_selector.RS = RectangleSelector(axes, line_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1, 3],  # don't use middle button
                                           minspanx=5, minspany=5,
                                           spancoords='pixels',
                                           interactive=True)
    plt.connect('key_press_event', toggle_selector)
    plt.show(block=False)
    plt.pause(8.0) 
 
    rect_selection_coords = toggle_selector.RS.extents
    plt.savefig(sPath + '\\' + filename +'.png', bbox_inches='tight',pad_inches=0)
    plt.close(fig)
#    input('Press Enter')
    
    return (rect_selection_coords)