# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 11:58:03 2022

@author: Ben
"""

from ScanImageTiffReader import ScanImageTiffReader

def imageMetaData(data_path, files_tiff):    
    
    filename = data_path + r'/' + files_tiff[0]
    # vol = ScanImageTiffReader(filename).data();
    # des = ScanImageTiffReader(filename).description(0);
    md = ScanImageTiffReader(filename).metadata();
    
    frames = md.find('framesPerSlice')
    slices = md.find('numSlices')
    zoom = md.find('scanZoomFactor')
    acqState = md.find('acqState')
    trig = md.find('trigAcqEdge')

    zStackStep = md.find('stackZStepSize')

    lines = [md[slices : slices+14] , md[zoom:zoom+20] , md[acqState:acqState+18], md[trig:trig+24], md[frames:frames+18], md[zStackStep:zStackStep+21]]

    filename = data_path + '/Image meta data.txt'
    with open(filename, 'w') as f:
        f.writelines(lines)
    f.close()
    
    zoom = md[zoom:zoom+20]
    frames = md[frames:frames+18]
    
    import re
    numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
    rx = re.compile(numeric_const_pattern, re.VERBOSE)
    rx.findall(zoom)
    zoomFactor = float(rx.findall(zoom)[0])
    print(' Zoom factor = ' + str(zoomFactor))

    rx.findall(frames)
    frames = float(rx.findall(frames)[0])
    print(' Number of frames = ' + str(frames))
    
    return zoomFactor, frames,  md[zStackStep:zStackStep+21] , md