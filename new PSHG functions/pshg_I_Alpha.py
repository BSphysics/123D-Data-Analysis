# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 22:44:36 2021

@author: Ben
"""
import numpy as np

def I_alpha(alpha, a0, a2, b2, a4, b4):    
    return a0 + a2*np.cos(2*alpha/180*np.pi) + b2*np.sin(2*alpha/180*np.pi) + a4*np.cos(4*alpha/180*np.pi) + b4*np.sin(4*alpha/180*np.pi)