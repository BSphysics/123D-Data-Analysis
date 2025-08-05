# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 15:54:57 2022

@author: Ben
"""

import xlsxwriter
import numpy as np

def BBOIntensityDataWriter(data_path, angle, angles, pshg):
    # workbook = xlsxwriter.Workbook(data_path + '/' + (data_path[-17:-14]) + '_BBO.xlsx')
    workbook = xlsxwriter.Workbook(data_path + '/' + angle + '_BBO.xlsx')
    worksheet = workbook.add_worksheet()
    bold = workbook.add_format({'bold': 1})
    
    worksheet.write('A1', 'BBO angle', bold)
    worksheet.write('B1', angle, bold)
    worksheet.write('B2', 'Angles (deg)', bold)
    worksheet.write('C2', 'Intensity', bold)
    
    row = 2
    col = 1
    for item in angles:
        worksheet.write_number(row, col, item)
        row += 1
    
    sumData = np.sum(np.sum(pshg,2),1)
    row = 2
    col = 2    
    for item in sumData:
        worksheet.write_number(row, col, item)
        row += 1
        
    workbook.close()