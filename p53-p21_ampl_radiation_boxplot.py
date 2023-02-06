#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:26:32 2023

@author: nicagutu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyboat import WAnalyzer

plt.rcParams.update({'font.size': 26})    
plt.rcParams['svg.fonttype'] = 'none'

path = '/Users/nicagutu/Nextcloud/GranadaLab/Users/Nica/p53_Long_recordings_Data_analysis/'
path2 = '/Users/nicagutu/Nextcloud/Manuscripts/p53Dynamics/Figures/Fig2/'

file_list = ['0G','2G','4G','10G_MD','10G_NMD']
labels = ['0Gy','2Gy','4Gy','10Gy']

#Setting the parameters necessary for the Wavelet spectrum analysis
dt = 0.5 # the sampling interval, 0.5hours
lowT = 2
highT = 10
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

dict_ampl_p53 = {}
dict_ampl_p21 = {}
for i in range(4):
    file = file_list[i]
    ampl_p53 = []
    ampl_p21 = []
    
    datap53 = pd.read_csv(path+'p53_'+file+r'.csv', header = None) 
    datap21 =pd.read_csv(path+'p21_'+file+r'.csv', header = 0, index_col = 0) 
    datap53.columns = datap21.columns
            
    for column in datap53:
        signalp53 = datap53[column]
        detrended_signalp53 = wAn.sinc_detrend(signalp53, T_c = 50)
        wAn.compute_spectrum(detrended_signalp53, do_plot = False) 
        wAn.get_maxRidge(power_thresh = 0, smoothing_wsize = 4) 
        rdp53 = wAn.ridge_data

        signalp21 = datap21[column]
        wAn.compute_spectrum(signalp21, do_plot = False) 
        wAn.get_maxRidge(power_thresh = 0, smoothing_wsize = 4) 
        rdp21 = wAn.ridge_data
        
        ampl_p53.append(np.median(rdp53['amplitude']))
        ampl_p21.append(np.median(rdp21['amplitude']))
    
    dict_ampl_p53[i] = ampl_p53
    dict_ampl_p21[i] = ampl_p21

df_ampl_p53 = pd.DataFrame.from_dict(dict_ampl_p53, orient='index')
df_ampl_p53 = df_ampl_p53.transpose()
df_ampl_p53.columns = labels
# print(df_ampl_p53)

df_ampl_p21 = pd.DataFrame.from_dict(dict_ampl_p21, orient = 'index')
df_ampl_p21 = df_ampl_p21.transpose()
df_ampl_p21.columns = labels
# print(df_ampl_p21)

fig = plt.figure(figsize = (10,10))
df_ampl_p53.boxplot(showfliers = False)
plt.ylabel('p53 amplitude (a.u.)')
# plt.savefig(path2+'p53_boxplot_ampl_doses.svg')
plt.show()

fig = plt.figure(figsize = (10,10))
df_ampl_p53.boxplot(showfliers = False)
plt.ylabel('p21 amplitude (a.u.)')
# plt.savefig(path2+'p21_boxplot_ampl_doses.svg')
plt.show()

