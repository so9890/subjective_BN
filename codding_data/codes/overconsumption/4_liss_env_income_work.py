#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  2 17:28:26 2022

@author: sonja

This script merges the waves of the working data set to the environmental and
income data set merged in 3_liss.

The script prepares figures of
    1) the evolution of the percentage of the working population below 36hours (parttime or below) 
        for reasons of stronger preferences for leisure
    2) the evolution of the percentage which would want to stop working for leisure reasons
    3) What determines whether a household wants to work less voluntarily? => Probit
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from itertools import compress


""" read in data from pickles (env +income) and working data set;
    merge.
    """
    
# Get list of all files in work folder
path="../../liss_data/work"
dir_list = os.listdir(path)

# drop zip elements in list
res =  [True for i in range(len(dir_list))]
for i in range(0, len(dir_list)):
    res[i]=str.__contains__(dir_list[i], '.dta')

list_waves_work = list(compress(dir_list, res))

# initiate dictionary to save data
dic_data=dict.fromkeys(list_waves_work)
# also save only panel data
dic_panel_con=dict.fromkeys(list_waves_work)

# loop over waves to merge year specific datasets
# save to pickles
for i in list_waves_work:
    
    # read in prepared dataset (env+income)
    data= pd.read_pickle('../../liss_data/merged/env_income_wave_20'+i[2:4])
    #read in working time data
    df_help = pd.read_stata('../../liss_data/work/'+i)
    
    # rename variables to be same across dataframes
    variable_names=df_help.columns.to_list()
    for j in range(0,len(variable_names)):
        #strr=list_occur[j]
        "variable of interest same across waves"
        variable_names[j]=variable_names[j].replace(i[2:5], "") 
    df_help.columns=variable_names
   
    # merge panel wave to environmental dataset
    # keep all observations also if not in wave \ar missing date!
    helper=data.merge(df_help,  left_on= 'nomem_encr', right_on= 'nomem_encr', how= 'left', validate='1:1', indicator='source_work')
    
    # save to pickles as env_income_bywave
    helper.to_pickle('../../liss_data/merged/env_income_work_wave_20'+i[2:4])
    
    dic_data[i]=helper
    
    # save wave
    dic_panel_con[i]=df_help
    
# concatenate from dictionary to long dataset 
data_long=pd.concat(list(dic_data.values()), join="outer")


""" Plot variables over time: percentage working voluntarily below 36 hours 
    or willing to reduce

    1) total sample 
    2) by environmental variable?

    """

# construct percentage variable
cats=['actual', 'willing']

# dictionary with meaning of each variable
list_variables=[['cw399', 'cw400'], \
         ['cw292', 'cw294', 'cw295']]
dicc=dict(zip(cats, list_variables))    

# drop unemployed, i.e. not receiving payment or nan working status
helper=data_long[(data_long['cw001']=='yes') | (data_long['cw001']=='Yes')] # no info on whether paid or not

# generate indicator for joint share (all 2(3) variables in dicc jointly)
for i in cats:
    
    varrs=dicc[i] 

    # generate indicator
    helper[i+'_ind']=0

    # loop over variables which determine voluntary reduction => indicator=1 if any of the 2(3) variables is affirmed
    for b in varrs:
        if i==cats[0]:
            helper.actual_ind[(helper[b]=='yes') | (helper[b]=='Yes')]=1
        elif i== cats[1]:
            helper.willing_ind[(helper[b]=='yes') | (helper[b]=='Yes')]=1
    

""" 1) plot total sample"""
 
# total share evolution

data_total=helper[[ 'year', 'willing_ind', 'actual_ind']].groupby(['year'], as_index=False).sum()
total_all=helper[[ 'year', 'willing_ind', 'actual_ind']].groupby(['year'], as_index=False).size()
data_total=data_total.merge(total_all[[ 'year','size']],  left_on= ['year'], right_on= [ 'year'], how= 'left', validate='1:1', indicator='source')

# generate share for both variables
for i in cats:
    data_total[i+'share']=data_total[i+'_ind'].values/data_total['size']

    fig = plt.figure()
    plt.plot(data_total.year, data_total[i+'share'])
    plt.title(i)
    #plt.legend([])
    plt.xticks(rotation=30, ha='right', fontsize = 12)
    plt.yticks( fontsize = 12)
    plt.savefig('../../results/liss/total_share_voluntary_work_reduction_'+i+'.png', format='png', bbox_inches='tight')
    plt.show()
    plt.clf()


""" 2) plot share by environmental indicator """

for v in ['qk20a175', 'qk20a135', 'qk20a181', 'qk20a183', 'qk20a144', 'qk20a148']:
    
    # only two groups agree and disagree
    helper['group_broad']='Others'
    helper.group_broad[(helper[v]=='Completely disagree') | (helper[v]=='Disagree')]='Disagree'
    helper.group_broad[(helper[v]=='Completely agree') | (helper[v]=='Agree')]='Agree'
    
    # collapse dataset by v
    data_col=helper[['group_broad', 'year', 'willing_ind', 'actual_ind']].groupby(['group_broad', 'year'], as_index=False).sum()
    total=helper[['group_broad', 'year', 'willing_ind', 'actual_ind']].groupby(['group_broad', 'year'], as_index=False).size()
    data_col=data_col.merge(total[['group_broad', 'year','size']],  left_on= ['group_broad', 'year'], right_on= ['group_broad', 'year'], how= 'left', validate='1:1', indicator='source')
    
    #generate share for both variables
    for i in cats:
        data_col[i+'share']=data_col[i+'_ind'].values/data_col['size']
     
    # list of categories in grouping variable
    groups=list(data_col.group_broad.drop_duplicates())
    
    # list of dataframes split
    frames=[data_col[data_col.group_broad==groups[0]]]
    for s in range(1,len(groups)):
        frames.append(data_col[data_col.group_broad==groups[s]])
    
    # plot time series
    for i in cats:
        fig = plt.figure()
        for frame in frames:
            plt.plot(frame['year'], frame[i+'share'])
        plt.legend(groups)
        #plt.title(dicc[i]+': no, not necessary (in percent)')
        plt.xticks(rotation=30, ha='right', fontsize = 12)
        plt.yticks( fontsize = 12)
    
        plt.savefig('../../results/liss/broad_groups_work_redcuction'+v+'_'+i+'.png', format='png', bbox_inches='tight')
        plt.show()
        plt.clf()
    
