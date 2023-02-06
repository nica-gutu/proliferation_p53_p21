#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 15:49:03 2023

@author: nicagutu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({'font.size': 22})
plt.rcParams['svg.fonttype'] = 'none'
plt.style.use('seaborn-dark-palette')

path = '/Users/nicagutu/Nextcloud/GranadaLab/Users/Nica/p53_Long_recordings_Data_analysis/'
path2 = '/Users/nicagutu/Nextcloud/Manuscripts/p53Dynamics/Figures/Fig2/'

file_list = ['0G', '2G', '4G', '10G_MD']
labels = ['0Gy','2Gy','4Gy','10Gy']

imt_dict={}
for i in range(4):  
    file='division_matrix_'+str(file_list[i])    
    data=pd.read_csv(path+file+r'.csv')
    data=data.T
    
    imt = []
    for column in data:
        division_time=((np.where(data[column] == 1)[0]))
        if len(division_time)>=0:
            for ii in range(len(division_time)-1):
                diff = (division_time[ii+1]-division_time[ii])*0.5
                if 10 <= diff<= 50:
                    imt.append(diff)
    imt_dict[i] = imt
    
df = pd.DataFrame.from_dict(imt_dict,orient = 'index')
df = df.T
df.columns = labels

fig = plt.figure(figsize = (12,10))
sns.violinplot(data = df)
plt.ylabel('Cell cycle length (h)')
# plt.savefig(path2+'Cell_cycle_length.svg')
plt.show()  
