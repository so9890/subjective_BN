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

    - only keep variables relevant to answer the question
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
for i in cons:
   data_long[i+'_cat'] = data_long[i].astype('category').cat.codes
   # no i dont think it is necessary is encoded as 1

# drop nans => same for both variables
helper=data_long[~pd.isna(data_long['ci307'])]

helper.ci306_cat[helper.ci306_cat!=1]=0
helper.ci307_cat[helper.ci307_cat!=1]=0

# collapse dataset by willingness to change
data_col=helper[['qk20a175', 'year', 'ci306_cat', 'ci307_cat']].groupby(['qk20a175', 'year'], as_index=False).sum()
total=helper[['qk20a175', 'year', 'ci306_cat', 'ci307_cat']].groupby(['qk20a175', 'year'], as_index=False).size()
data_col=data_col.merge(total[['qk20a175', 'year','size']],  left_on= ['qk20a175', 'year'], right_on= ['qk20a175', 'year'], how= 'left', validate='1:1', indicator='source')

#generate share for both variables
for i in cons:
    data_col[i+'share']=data_col[i+'_cat'].values/data_col['size']
 

# dictionary with meaning of each variable
meaning=['Do you buy new clothes regularly?', \
         'Do you replace worn furniture with new furniture?']
dicc=dict(zip(cons, meaning))    

# list of categories in grouping variable
groups=list(data_col.qk20a175.drop_duplicates())

# list of dataframes split
frames=[data_col[data_col.qk20a175==groups[0]]]
for s in range(1,len(groups)):
    frames.append(data_col[data_col.qk20a175==groups[s]])
    
    
    # plot time series
for i in cons:
    fig = plt.figure()
    for frame in frames:
        plt.plot(frame['year'], frame[i+'share'])
    plt.legend(groups)
    plt.title(dicc[i])
    plt.show()

   # plt.savefig('../../results/liss/conditional_heatmap'+i+'_'+j+'labels'+str(p)+'.png', format='png', bbox_inches='tight')
    #plt.clf()