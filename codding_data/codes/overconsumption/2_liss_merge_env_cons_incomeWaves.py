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

# consumption data => merge individually to environmental data file and later append to get 2-indices structure
# save to dictionary
list_consumption_data=['09a', '10b', '12c', '15d', '17e']
dic_data=dict.fromkeys(list_consumption_data)
for i in list_consumption_data:
    df_help = pd.read_stata('../../liss_data/timeUse_con/bf'+i+'_EN_1.0p.dta')
    
    # generate year variable
    df_help.year=
    dic_data[i]=df.merge(df_help,  left_on= 'nomem_encr', right_on= 'nomem_encr', how= 'left', validate='1:1', indicator='source'+i)



