#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 16:46:14 2022

@author: sonja

Script to run 
    
    (a) average characteristics of voluntary reducers
    (b) POOLED logit/probit regression of choices:
    - more consumption not necessary
    - voluntary working time reduction
     => informative on typical features of those which reduce consumtpion or work
     => also check if this is the same group of people
"""
import pandas as pd
import matplotlib.pyplot as plt
import os
from itertools import compress
import seaborn as sns
from matplotlib.colors import ListedColormap

from functions_liss import dont
import math

import warnings

from mpl_toolkits.mplot3d import Axes3D

from scipy import stats
from statsmodels.formula.api import logit
from statsmodels.formula.api import probit
from statsmodels.formula.api import glm


import statsmodels.api as sm

""" read in data """

data= pd.read_pickle('../../liss_data/merged/env_income_work_long')
indics= pd.read_pickle('../../liss_data/merged/indic_consumption_reduc') # note: does not contain unemployed!!

data=data.merge(indics,  left_on= ['nomem_encr', 'year'], right_on= ['nomem_encr', 'year'], how= 'left', validate='1:1', indicator='source_data')

""" Average economic characteristics of those which at some point 
    voluntarily consume or work less
    
    - static: pooled
    - dynamic: those who reduce/ transition to consume less
        -at time of reduction
        - on average before reduction! (after reduction can have changed)
    """
"""Static: 
    1) who are part time workers?/ not necessary
      a) all (at some time in panel)
      b) those which become part time workers and remain till end of panel before they reduce => is there a role for environmental concerns?
   2) variables: 
        - wage (annual net wage: ci08a012, ci08a021, ci08a03; wage period in last job: cw08a322; cw08a323 through cw08a325 wage in last job )
        - eduction: cw08a005 (ordinary); if other then next variable needed: cw08a006
        - sector: cw08a402; 
        - profession: cw08a404 => informative on skills!
        - family size
        - sex 
        - married
        - head of household
        - age cw08a003, year of birth cw08a002
    """
 
""" Sort skills into green, neutral or brown in 2020 (info on env. attitudes)

    Based on Consoli et al. 2016: they find that an occupation which is green has a higher level of 
    non-routine analytic, non-routine interactive, non-routine manual work
    
    What matters is not if the actual work is green but rather if they could work in a green sector
    """
    
# skill variable (cw404) coded as
#          1/2 higher academic, independent profession, supervisory;
#          3+4 intermediate academic independent; 
#           5 other mental work
#           6 skilled manual work, 7 semiskilled manual work, 8 unskilled manual work
#           9 agrarian

helpp = data[data.year=='2020']
# index of categories in skill variable
cats=helpp['cw404'].cat.categories
# list of green skills
list_green=[cats[0], cats[1] ,cats[3], cats[5]]
list_neutral=[cats[2], cats[8]]
list_lessgreen=[cats[4], cats[6], cats[7]]
# less green job =0 : unskilled mental work, 
helpp['green_skill']=0
# set indicator equal to 1 if in neutral skill, and equal to 2 if green skill
helpp.green_skill[helpp['cw404'].isin(list_green)]=2
helpp.green_skill[helpp['cw404'].isin(list_neutral)]=1
# test
helpp.green_skill[helpp['cw404'].isin(list_lessgreen)].unique()==0

# code as categorical  
helpp.green_skill=helpp.green_skill.astype('category')  

"""Correlations skill and 
    1) environmental attitudes in Crosssection,
    2) part time work, reduction consumption
    """

# Heatmaps for skills (cw404): conditional and joint distributions
# only keep 2020 -> current state


list_env= ['qk20a175', 'qk20a135', 'qk20a181', 'qk20a183', 'qk20a144', 'qk20a148']

# frequency headmap: marginal distribution: given skill what do they choose?

dic_kind={'index': 'conditional', 'all': 'joint'}
for varrs in ['cw404', 'green_skill']: # over separate skills and over categoroies green, neutral, less green

    for i in list(dic_kind.keys()):
        
        for v in list_env:
            helpps= pd.crosstab( helpp[varrs], helpp[v],  normalize=i, margins = i=='all') # margins false if i==index

           # convert columns and row names to srings
            helpps.columns = helpps.columns.astype(str)
            helpps.index = helpps.index.astype(str)
            if varrs == 'green_skill': 
                helpps.index =['less green', 'neutral', 'green']
                
     
            helpps=dont(helpps)
           
            data1=helpps.copy()
            if i=='all':
                data1['All'] = float('nan')   # columns
            data1.loc['All'] = float('nan')    # rows
            ax = sns.heatmap(data1, annot=True, cmap="BuPu")
            
             # Greens
            data2 = helpps.copy()
            #data2.columns = data2.columns.add_categories('All')
            if i=='index':
                data2['All']=1
                data2.iloc[:,:-1] = float('nan')
            elif i=='all':
                data2.loc['All'] = float('nan')    # rows
                data2.iloc[:-1,:-1] = float('nan')
                
            sns.heatmap(data2, annot=True, cbar=False, cmap=ListedColormap(['white']))
            plt.savefig('../../results/liss/heatmap_skills_'+dic_kind[i]+v+varrs+'.png', format='png', bbox_inches='tight')
            plt.show()
            plt.clf()

"""evolution over time skill and voluntary reduction """

# collapse data set by skills
data_col=data[['cw404', 'year', 'actual_ind', 'ci306_cat', 'ci307_cat']].groupby(['cw404', 'year'], as_index=False).sum()
total=data[['cw404', 'year', 'actual_ind', 'ci306_cat', 'ci307_cat']].groupby(['cw404', 'year'], as_index=False).size()
data_col=data_col.merge(total[['cw404', 'year','size']],  left_on= ['cw404', 'year'], right_on= ['cw404', 'year'], how= 'left', validate='1:1', indicator='source')

#generate share for indicator variables
for i in ['actual_ind', 'ci306_cat', 'ci307_cat']:
    data_col[i+'share']=data_col[i].values/data_col['size']
 
# list of categories in grouping variable
groups=list(data_col.cw404.drop_duplicates())

# list of dataframes split
frames=[data_col[data_col.cw404==groups[0]]]
for s in range(1,len(groups)):
    frames.append(data_col[data_col.cw404==groups[s]])

# plot time series
for i in 'actual_ind', 'ci306_cat', 'ci307_cat':
    fig = plt.figure()
    for frame in frames:
        plt.plot(frame['year'], frame[i+'share'])
    plt.legend(groups, bbox_to_anchor=(1,1), loc="upper left")
    #plt.title(dicc[i]+': no, not necessary (in percent)')
    plt.xticks(rotation=30, ha='right', fontsize = 12)
    plt.yticks( fontsize = 12)

    plt.savefig('../../results/liss/voluntary_behaviour_byskill_'+i+'.png', format='png', bbox_inches='tight')
    plt.show()
    plt.clf()

"""Heatmaps skill and behaviour"""
# generate indicator if member reduces at some point in time
# within panelist take sum of indicators if not necessary furniture, clothes, or part time work for leisure
indics_reduc=data[['nomem_encr','ci306_cat', 'ci307_cat', 'actual_ind']].groupby(['nomem_encr'], as_index=False).max()
indics_reduc.columns=['nomem_encr', 'nec_furniture_ind', 'nec_clothes_ind', 'parttime_ind']

for v in indics_reduc.columns[1:]:
    helpp=indics_reduc[v]
    helpp[helpp==1]='voluntarily consuming less'
    helpp[helpp==0]='nope'
    indics_reduc[v]=indics_reduc[v].astype('category')

#merge indicator to dataset
data=data.merge(indics_reduc,  left_on= ['nomem_encr'], right_on= ['nomem_encr'], how= 'left', validate='m:1')

# correlations with skills=> heatmaps of whether household voluntarily consumes less and parttime work
# for 2020 and 2021

dic_kind={'index': 'conditional', 'all': 'joint'}
for y in ['2020', '2021']:
    helpp=data[data.year==y]
    for i in list(dic_kind.keys()):
        
        # loop over indicators of behaviour
        for v in list(indics_reduc.columns[1:]):
            helpps= pd.crosstab( helpp['cw404'], helpp[v],  normalize=i, margins = i=='all') # margins false if i==index
            #helpps=dont(helpps)
            
            data1=helpps.copy()
            #if i=='joint':
            #    data1['All'] = float('nan')   # columns
            #data1.loc['All'] = float('nan')     # rows
            sns.heatmap(data1, annot=True, cmap="BuPu")
            
             # Greens
            data2 = helpps.copy()
            #data2.columns = data2.columns.add_categories('All')
            #if i=='index':
             #   data2['All']=1
              #  data2.iloc[:,:-1] = float('nan')
            #elif i=='all':
            #    data2.iloc[:-1,:-1] = float('nan')
                
            #sns.heatmap(data2, annot=True, cbar=False, cmap=ListedColormap(['white']))
            plt.savefig('../../results/liss/heatmap_skills_'+dic_kind[i]+v+y+'.png', format='png', bbox_inches='tight')
            plt.show()
            plt.clf()
            
"""Pooled Panel: Probit/Logit: What explains the probability to xxx
    heteroskedasticity of error terms due to panel structure
    Outcome variables: ci306_cat, ci307_cat, actual_ind, hours worked
    """


"""Dynamic: Local projection Panel
    1) create indicator to before and after reduction on individual level
        How to treat those which increase hours again?
        Need to treat those which work more hours again differently
        Include hours worked dynamically as outcome variable and voluntary as a regressor
        
        - actually hours worked, used to work in last job: cw08a127 in most important job
        - hours per week per contract in most important job
        """
        




    

""" What explains the probability to think new clothes/furniture are not necessary?

    1. step: economic/ demographics observables:
    regressors: income, education, family size, change in family size
    
    2. step: errors on environmental concerns (not sure this works with categorical variable)
    """
