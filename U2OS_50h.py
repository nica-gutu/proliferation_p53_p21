#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 12:22:42 2023

@author: nicagutu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import pearsonr

plt.rcParams.update({'font.size': 22})
plt.rcParams['svg.fonttype'] = 'none'

path = '/Users/nicagutu/Nextcloud/Manuscripts/p53Dynamics/Data/IFLiveData/untreated_data_GFP_proliferation/'
path2 = '/Users/nicagutu/Nextcloud/Manuscripts/p53Dynamics/Figures/Fig1/'

dishes = ['2','3','4','5']
color_list=['tab:blue','tab:green','tab:orange','tab:purple']

gfp_values = []
IDs = []
metrics1 = {'gfp':[], 'imt':[]}
metrics2 = {'time_last_division': {0:[],1:[],2:[],3:[]}, 'gfp':{0:[],1:[],2:[],3:[]}}
df_list = []
     
for i in dishes:
    #Input data
    gfp_file = 'gfp_dish'+str(i)+'_untreated'
    divisions_file = 'divisions_dish'+str(i)+'_untreated'
    
    data_gfp = pd.read_excel(path+gfp_file+r'.xlsx',header=None)
    data_gfp = data_gfp.T
    divisions = pd.read_excel(path+divisions_file+r'.xlsx',header=None) 
    divisions = divisions.T
    
    #Dataframe cleaning and filtering
    for col in data_gfp:
        if data_gfp[col][0] == -1:
            data_gfp=data_gfp.drop(col, axis = 1)
            divisions=divisions.drop(col, axis = 1)
            
    for col in divisions:
        division_times=((np.where(divisions[col] == 1)[0]))
        if len(division_times) != 0:
            if len(division_times) > 1:
                for ii in range(len(division_times)-1):
                    if (division_times[ii+1]-division_times[ii]) < 18:
                        if col in divisions:
                            divisions=divisions.drop(col, axis = 1) 
                                                        
    #Looking for the IMT, time of last division and GFP values
    for col in divisions:
        gfp_values.append(data_gfp[col][0])
        IDs.append(str(i)+'_'+str(col))
        division_times=((np.where(divisions[col] == 1)[0]))
        if len(division_times)>1:
            imt=[]
            for ii in range(len(division_times)-1):
                imt.append(division_times[ii+1]-division_times[ii])
            metrics1['imt'].append(np.mean(imt))
            metrics1['gfp'].append(data_gfp[col][0])
            
        division_count = len(np.where(divisions[col] == 1)[0])
        if division_count>=1:
            metrics2['time_last_division'][division_count].append(54 - division_times[-1])
            metrics2['gfp'][division_count].append(data_gfp[col][0])
        else:            
            metrics2['time_last_division'][division_count].append(54)
            metrics2['gfp'][division_count].append(data_gfp[col][0])
            
    divisions.columns = [str(i) + '_' + str(col) for col in divisions.columns]
    df_list.append(divisions)

divisions = pd.concat(df_list, axis=1)

df_gfp = pd.DataFrame().from_dict(metrics2['gfp'], orient = 'index')
df_gfp = df_gfp.T

# #%%Plotting

fig = plt.figure(figsize = (12,10))
sns.boxplot(data = df_gfp,showfliers=False)
plt.ylabel('GFP intensity (a.u.)')
plt.xlabel('Number of divisions')
# plt.savefig(path2+'GFP_proliferation_boxplot.svg')
plt.show()

print('Mean of the GFP values ',np.mean(gfp_values))
print('Standard deviation of the GFP values ',np.std(gfp_values))
fig = plt.figure(figsize=(12,10))
sns.histplot(data = gfp_values,kde=True,stat='density')
plt.xlabel('GFP values (a.u.)')
# plt.savefig(path2+'GFP_histogram_untreated.svg')
plt.show()

fig = plt.figure(figsize = (12,10))
for i in range(4):
    plt.plot(metrics2['time_last_division'][i], metrics2['gfp'][i],'o',label = str(i)+' divisions')
plt.ylabel('GFP intensity (a.u.)')
plt.xlabel('Time since last division (h)')
plt.ylim([0,600])
plt.legend(loc = 'best')
# plt.savefig(path2+'GFP_vs_TimeLastDivision.svg')
plt.show()

corr_coef, _ = pearsonr(metrics1['imt'], metrics1['gfp'])
fig = plt.figure(figsize = (12,10))
plt.plot(metrics1['imt'], metrics1['gfp'],'o')
plt.text(25, 500, 'Correlation = %0.5f'%(corr_coef))
plt.ylabel('GFP intensity (a.u.)')
plt.xlabel('IMT (h)')
plt.ylim([0,600])
# plt.savefig(path2+'GFP_vs_IMT.svg')
plt.show()

fig = plt.figure(figsize = (12,10))
for i in range(1,4):
    bins_num = int(max(metrics2['time_last_division'][i])-min(metrics2['time_last_division'][i]))
    print(bins_num)
    sns.histplot(metrics2['time_last_division'][i],stat='density',kde = True, bins = bins_num,label = str(i)+' division',color=color_list[i])
plt.xlabel('Cell age (h)')
plt.legend(loc = 'best')
# plt.savefig(path2+'Cell_age_number_divisions_hist.svg')
plt.show()


#%%Division profile
    
num_Ids = len(divisions.columns)
time_points = len(divisions.index.values) 
properties = {'num_divisions':[],'time_1st':[]}

for col in IDs:
    properties['num_divisions'].append(sum(divisions[col]))
    count=0
    for jj in range(len(divisions[col])):
        if sum(divisions[col]) == 0 and count <1:
            count+=1
            properties['time_1st'].append(0)
        elif sum(divisions[col]) != 0:
            if divisions[col][jj] == 1 and count <1:
                properties['time_1st'].append(divisions.index.values[jj])
                count+=1
            
properties['num_divisions'],properties['time_1st'],IDs=zip(*sorted(zip(properties['num_divisions'],properties['time_1st'],IDs)))

divisions = divisions.reindex(columns=IDs)
division_profile = (divisions.to_numpy())
division_profile = division_profile.T

for xx in range(num_Ids):
    divisions = []
    for yy in range(time_points):
        if division_profile[xx,yy] == 1:
            divisions.append(yy)
    if len(divisions) >0:
        for n in range(len(divisions)-1):
            division_profile[xx,divisions[n]+1:divisions[n+1]] = division_profile[xx,divisions[n]+1:divisions[n+1]]+(n+1)   
        division_profile[xx,divisions[-1]+1:time_points] = division_profile[xx,divisions[-1]+1:time_points]+len(divisions)
    if len(divisions) >1:    
        division_profile[xx,divisions[-1]] = division_profile[xx,divisions[-1]]+1
 
plt.figure(figsize=(8, 14))
color_map=plt.imshow(division_profile)
color_map.set_cmap("plasma")
plt.xlabel('Time(h)')
plt.ylabel('Cell number')
ax=plt.gca()
divider=make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cb=plt.colorbar(cax=cax)
labels_list=np.arange(0,sum(divisions),1)
loc=labels_list+0
cb.set_ticks(loc)
cb.set_ticklabels(labels_list)
cb.set_label('No. divisions',rotation=270, labelpad=25)
# plt.savefig(path2+'Division_profile_U2OS.svg')
plt.show()

 