#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 12:27:44 2021

@author: sonja

script to connect environmental data set and income data sets (all waves)

Goal: 1) How does the variable of thinking furniture replacememnt and new clothes vary over time by different income groups?

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from itertools import compress

"""read in data and merge to environmental data set """

# environmental data 
dtafile= '../../liss_data/environment/qk20a_EN_1.0p.dta'
df = pd.read_stata(dtafile)

""" read in and merge data on income. 

    - only keep variables relevant to answer the question in income dataset
    - ensure appending correct variables to one another (changes in variable names)

"""

# merge individually to environmental data file and later append to get 2-indices structure
# save to dictionary
# import OS module

# Get list of all files and directories
path="../../liss_data/income"
dir_list = os.listdir(path)

# drop zip elements in list
res =  [True for i in range(len(dir_list))]
for i in range(0, len(dir_list)):
    res[i]=str.__contains__(dir_list[i], '.dta')

list_waves = list(compress(dir_list, res))

# initiate dictionary to save data
dic_data=dict.fromkeys(list_waves)
# also save only panel data
dic_panel_con=dict.fromkeys(list_waves)

for i in list_waves:
    df_help = pd.read_stata('../../liss_data/income/'+i)
    
    # rename variables to be same across dataframes
    variable_names=df_help.columns.to_list()
    for j in range(0,len(variable_names)):
        #strr=list_occur[j]
        "variable of interest same across waves"
        variable_names[j]=variable_names[j].replace(i[2:5], "") 
    df_help.columns=variable_names
    
    # merge panel wave to environmental dataset
    # keep all observations also if not in wave \ar missing date!
    helper=df.merge(df_help[['nomem_encr','ci307', 'ci306']],  left_on= 'nomem_encr', right_on= 'nomem_encr', how= 'left', validate='1:1', indicator='source')
    # create year variable
    helper['year']=df_help.ci_m.astype(str).str[0:4]
    #helper['month']=df_help.ci_m.astype(str).str[4:6]
    dic_data[i]=helper
    
    # save wave
    dic_panel_con[i]=df_help
    
# concatenate from dictionary to long dataset 
data_long=pd.concat(list(dic_data.values()), join="outer")

# move year of panel interview to second column
helpper = data_long.pop('year') # removes column b from frame and saves to new dataframe
data_long.insert(loc = 1, column = 'year', value = list(helpper))


""" Plot variables over time by environmental indicator. 

    for categorical variables: percentage not necessary
    """

# generate share of those which think new items are not necessary by group
# ensure not to include nans: drop nans in ci306 and ci307
# save to dictionary

cons=['ci306', 'ci307']

# dictionary with meaning of each variable
meaning=['Do you buy new clothes regularly?', \
         'Do you replace worn furniture with new furniture?']
dicc=dict(zip(cons, meaning))    

for i in cons:
   data_long[i+'_cat'] = data_long[i].astype('category').cat.codes
   # no i dont think it is necessary is encoded as 1

# drop nans => same for both variables
helper=data_long[~pd.isna(data_long['ci307'])]

helper.ci306_cat[helper.ci306_cat!=1]=0
helper.ci307_cat[helper.ci307_cat!=1]=0

# drop 2008-2010 as no observation from environmental dataset in panel
helper=helper[helper.year.astype(float)>2010]

""" collapse data set by distinct env. variable indicators. """

for v in ['qk20a175', 'qk20a135', 'qk20a181', 'qk20a183', 'qk20a144', 'qk20a148']:
    
    # collapse dataset by willingness to change
    data_col=helper[[v, 'year', 'ci306_cat', 'ci307_cat']].groupby([v, 'year'], as_index=False).sum()
    total=helper[[v, 'year', 'ci306_cat', 'ci307_cat']].groupby([v, 'year'], as_index=False).size()
    data_col=data_col.merge(total[[v, 'year','size']],  left_on= [v, 'year'], right_on= [v, 'year'], how= 'left', validate='1:1', indicator='source')
    
    #generate share for both variables
    for i in cons:
        data_col[i+'share']=data_col[i+'_cat'].values/data_col['size']
     
    
    # list of categories in grouping variable
    groups=list(data_col[v].drop_duplicates())
    
    # list of dataframes split
    frames=[data_col[data_col[v]==groups[0]]]
    for s in range(1,len(groups)):
        frames.append(data_col[data_col[v]==groups[s]])
        
        
    # plot time series
    for i in cons:
        fig = plt.figure()
        for frame in frames:
            plt.plot(frame['year'], frame[i+'share'])
        plt.legend(groups)
        plt.title(dicc[i])
        #plt.show()
    
        plt.savefig('../../results/liss/share_notnecessary'+v+'_'+i+'.png', format='png', bbox_inches='tight')
        plt.clf()
        
    """ plot total"""
    
    # total share evolution
    data_total=helper[[ 'year', 'ci306_cat', 'ci307_cat']].groupby(['year'], as_index=False).sum()
    total_all=helper[[ 'year', 'ci306_cat', 'ci307_cat']].groupby(['year'], as_index=False).size()
    data_total=data_total.merge(total_all[[ 'year','size']],  left_on= ['year'], right_on= [ 'year'], how= 'left', validate='1:1', indicator='source')
    
    # generate share for both variables
    for i in cons:
        data_total[i+'share']=data_total[i+'_cat'].values/data_total['size']
    
        fig = plt.figure()
        plt.plot(data_total.year, data_total[i+'share'])
        plt.title(dicc[i])
        plt.legend(['no, not necessary (in percent)'])
        #plt.show()
        plt.savefig('../../results/liss/total_share_notnecessary'+v+'_'+i+'.png', format='png', bbox_inches='tight')
        plt.clf()
        
    """ plot broader categories"""
    
    # only two groups agree and disagree
    helper['group_broad']='Others'
    helper.group_broad[(helper[v]=='Completely disagree') | (helper[v]=='Disagree')]='Disagree'
    helper.group_broad[(helper[v]=='Completely agree') | (helper[v]=='Agree')]='Agree'
    
    # collapse dataset by willingness to change
    data_col=helper[['group_broad', 'year', 'ci306_cat', 'ci307_cat']].groupby(['group_broad', 'year'], as_index=False).sum()
    total=helper[['group_broad', 'year', 'ci306_cat', 'ci307_cat']].groupby(['group_broad', 'year'], as_index=False).size()
    data_col=data_col.merge(total[['group_broad', 'year','size']],  left_on= ['group_broad', 'year'], right_on= ['group_broad', 'year'], how= 'left', validate='1:1', indicator='source')
    
    #generate share for both variables
    for i in cons:
        data_col[i+'share']=data_col[i+'_cat'].values/data_col['size']
     
    # list of categories in grouping variable
    groups=list(data_col.group_broad.drop_duplicates())
    
    # list of dataframes split
    frames=[data_col[data_col.group_broad==groups[0]]]
    for s in range(1,len(groups)):
        frames.append(data_col[data_col.group_broad==groups[s]])
        
        
    # plot time series
    for i in cons:
        fig = plt.figure()
        for frame in frames:
            plt.plot(frame['year'], frame[i+'share'])
        plt.legend(groups)
        plt.title(dicc[i]+': no, not necessary (in percent)')
        #plt.show()
    
        plt.savefig('../../results/liss/broad_groups_notnecessary'+v+'_'+i+'.png', format='png', bbox_inches='tight')
        plt.clf()
    
