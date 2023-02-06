#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:23:37 2023

@author: nicagutu
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyboat import WAnalyzer

plt.rcParams.update({'font.size': 22})
plt.rcParams['svg.fonttype'] = 'none'

path = '/Users/nicagutu/Nextcloud/GranadaLab/Users/Nica/p53_Long_recordings_Data_analysis/'
path2 = '/Users/nicagutu/Nextcloud/Manuscripts/p53Dynamics/Figures/Fig2/'

files = ['0G', '2G', '4G', '10G_MD']
labels = ['0Gy','2Gy','4Gy','10Gy']
colors = ['tab:blue','tab:orange','tab:green','tab:purple']

#Setting the parameters necessary for the Wavelet spectrum analysis
dt = 0.5 # the sampling interval, 0.5hours
lowT = 2
highT=10
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

df=pd.DataFrame()
for i in range(1):#len(files)):
    data = pd.read_csv(path+'p21_'+files[i]+r'.csv', header = 0, index_col = 0) 
    data2 = pd.read_csv(path+'p53_'+files[i]+r'.csv', header = None, index_col = None) 
    division_matrix=pd.read_csv(path+'division_matrix_'+files[i]+r'.csv', header = None)
    division_matrix=division_matrix.T
    data2.columns = data.columns
    division_matrix.columns = data.columns
    
    amp_dict = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    std_dict = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    
    for column in data:
        signalp53 = data2[column]
        detrended_signalp53 = wAn.sinc_detrend(signalp53, T_c = 50)
        wAn.compute_spectrum(detrended_signalp53, do_plot = False) 
        wAn.get_maxRidge(power_thresh = 0, smoothing_wsize = 4) 
        rdp53 = wAn.ridge_data

        signalp21 = data[column]
        detrended_signalp21 = wAn.sinc_detrend(signalp21, T_c = 50)
        
        division_times = (np.where(division_matrix[column] == 1)[0])
        num_div = len(division_times)
        amp_dict[num_div].append(np.mean(rdp53['amplitude']))
        std_dict[num_div].append(np.std(detrended_signalp21)) #or amplitude of p21
                       
    df_list=pd.DataFrame(columns=['0','1','2','3','4','5']).from_dict(amp_dict, orient='index')
    df_list=df_list.T
    
    fig=plt.figure(figsize = (12,10))
    df_list.boxplot(showfliers = False)
    plt.ylabel('p53 amplitude (a.u.)')
    plt.xlabel('Number divisions')
    plt.legend(loc = 'best')
    # plt.savefig(path2+'Boxplot_amplitude_p21.svg')
    plt.show() 
        
    df_list=pd.DataFrame(columns=['0','1','2','3','4','5']).from_dict(std_dict, orient='index')
    df_list=df_list.T
    
    fig = plt.figure(figsize = (12,10) )
    df_list.boxplot(showfliers = False)#color=colors_1[i])
    plt.ylabel('p21 STD of detrended signal (a.u.)')
    plt.xlabel('Number divisions')
    fig.set_figheight(10)
    fig.set_figwidth(12)
    plt.legend(loc = 'best')
    # plt.savefig(path2+'Boxplot_std_p21.svg')
    plt.show() 

    