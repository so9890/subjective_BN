#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 11:53:22 2021

@author: sonja

This file merges the dataset on environmental preferences and actions
to the liss waves on consumption, income, and demographics.

The resulting data set is used to
    a) (file: 3_liss) run a regression of changes in consumption on observables
        => the residuals are used to document unexplained changes in consumption
        => are these unexplained changes explained by willingness to change and take action?
    b) (file: 4_liss) descriptives on who are these people which are willing to change their
       livestyle or lower their consumption due to environmental reasons
"""

import pandas as pd
import numpy as np

"""read in data and merge to environmental data set """

# environmental data 
dtafile= '../../liss_data/environment/qk20a_EN_1.0p.dta'
df = pd.read_stata(dtafile)

""" merge consumption data waves and concatenate """

# merge individually to environmental data file and later append to get 2-indices structure
# save to dictionary
list_consumption_data=['09a', '10b', '12c', '15d', '17e']
dic_data=dict.fromkeys(list_consumption_data)
# also save only panel data
dic_panel_con=dict.fromkeys(list_consumption_data)

for i in list_consumption_data:
    df_help = pd.read_stata('../../liss_data/timeUse_con/bf'+i+'_EN_1.0p.dta')
    
    # rename variables to be same across dataframes
    variable_names=df_help.columns.to_list()
    for j in range(0,len(variable_names)):
        #strr=list_occur[j]
        variable_names[j]=variable_names[j].replace(str(i), "")
    df_help.columns=variable_names
    
    # merge panel wave to environmental dataset
    # keep all observations also if not in wave \ar missing date!
    helper=df.merge(df_help,  left_on= 'nomem_encr', right_on= 'nomem_encr', how= 'left', validate='1:1', indicator='source')
    # create year variable
    helper['year']='20'+i[:2]
    dic_data[i]=helper
    
    # save panel
    dic_panel_con[i]=df_help
    
# concatenate from dictionary to long dataset 
data_long=pd.concat(list(dic_data.values()), join="outer")

# move year of panel interview to second column
helpper = data_long.pop('year') # removes column b from frame and saves to new dataframe
data_long.insert(loc = 1, column = 'year', value = list(helpper))


""" Plot variables over time by environmental indicator. """

# other household expenditures: bf17e089 => durables?






























