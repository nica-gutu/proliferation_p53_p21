#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 17:48:55 2023

@author: nicagutu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyboat import WAnalyzer
import seaborn as sns

plt.rcParams.update({'font.size': 22})    
plt.rcParams['svg.fonttype'] = 'none'

def wavelet(signal):
    detrended_signal=wAn.sinc_detrend(signal, T_c = 50)
    wAn.compute_spectrum(detrended_signal, do_plot=False) #computes the detrended signal wavelet spectrum
    wAn.get_maxRidge(power_thresh=0, smoothing_wsize=4) #gets the maximum ridge
    rd=wAn.ridge_data
    return rd

path = '/Users/nicagutu/Nextcloud/GranadaLab/Users/Nica/p53_Long_recordings_Data_analysis/'
path2 = '/Users/nicagutu/Nextcloud/Manuscripts/p53Dynamics/Figures/Fig4/'

file_list = ['0G','2G','4G','10G_MD']
labels = ['0Gy','2Gy','4Gy','10Gy']

colors_1=['tab:blue','tab:orange','tab:green','tab:red']

#Setting the parameters necessary for the Wavelet spectrum analysis
dt = 0.5 # the sampling interval, 0.5hours
lowT = 2
highT = 10
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

properties={'amplitude': [], 'period': []}
for i in range(1):
    data = pd.read_csv(path+'p53_'+file_list[i]+r'.csv', header=None) 
    # data =pd.read_csv(path+'p21_'+file_list[i]+r'.csv', header = 0, index_col = 0) 
    time = data.index.values*0.5
            
    for column in data:
        signal = data[column]
        rd=wavelet(signal)
        
        properties['amplitude'].append(np.mean(rd['amplitude']))
        properties['period'].append(np.mean(rd['periods']))
    
fig = plt.figure(figsize = (12,10))
sns.histplot(properties['amplitude'], kde=True, stat='density')
plt.xlabel('Amplitude')
# plt.savefig(path2+'Amplitude_distr_p21.svg')
plt.show()
        
fig = plt.figure(figsize=(12,10))
sns.histplot(properties['period'], kde=True, stat='density')
plt.xlabel('Period [hours]')
# plt.savefig(path2+'Period_distr.svg')
plt.show()
        



