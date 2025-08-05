# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 11:34:58 2022

@author: Ben
"""
import numpy as np
import matplotlib.pyplot as plt
import os

plt.close('all')

data_path = r'C:\Users\Ben\OneDrive - University of Exeter\!Work\Work.2022\Data\123D\2022_03_09 pSHG fresh hydrated tendon birefringence depth\New FOV\Depth'

files = os.listdir(data_path)
files_np = [i for i in files if i.endswith('.npy')]

meanIntensity = []
angles = np.arange(0,14)*16

plt.figure(figsize = (12,6))
for filename in files_np:

    data = np.load(data_path+'\\'+filename)
    plt.plot(angles,data, label = filename[-7:-4] + ' um')
    meanIntensity.append(data)
# plt.legend(loc='lower left')    
plt.legend(bbox_to_anchor=(1.0,1), loc="upper left")
plt.xlabel('Polarisation angle (deg)')
plt.ylabel('Norm. intensity')
plt.title('SHG intensity vs polarisation angle for different depths in rat tail tendon')
meanIntensity = np.asarray(meanIntensity)

plt.figure()
plt.imshow(np.transpose(meanIntensity))
# plt.axis('off')
plt.colorbar()