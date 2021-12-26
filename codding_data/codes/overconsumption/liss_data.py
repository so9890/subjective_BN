
"""
Created on Sat Dec 25 14:12:05 2021

@author: sonja

File to read in Liss data
and generate descriptive statistics
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

"read in data"

dtafile= '../../liss_data/qk20a_EN_1.0p.dta'
varfile= '../../liss_data/variable_explanation.xls'
df = pd.read_stata(dtafile)
df.tail()
var_names=pd.read_excel(varfile, usecols=("A:B"), header=None)
var_names=var_names[~pd.isna(var_names[1])]

"get numbers of variable names which are to be plotted"
list_numbers=var_names[0].str[-3:].tolist()

"visualising data"

for i in list_numbers:
    "collapse dataset for variable of interest"
    occur=df.groupby(['qk20a'+str(i)]).size()
    total=occur.sum()
     
    # this is to check whether the sum in occur does not contain nans
    #df['s']=pd.isna(df['qk20a'+str(i)])
    #total=df.groupby(['s']).size()
    
    #percentages
    occur=occur/total*100
    
    # replace apostrophe with \' in index
    list_occur=occur.index.tolist()
    for j in range(0,len(occur.index)):
       # strr=list_occur[j]
        list_occur[j]=list_occur[j].replace("Don\x92t know", "Don\'t know")
    occur.index=list_occur
    
    "plot bar chart"    
    plt.bar(occur.index,occur)    
    # extract title
    helper=var_names[var_names[0]=='qk20a'+str(i)]
    plt.title(str(helper.iloc[0,1]), fontsize=15)
    #plt.xlabel('categories')
    plt.ylabel('percent', fontsize = 18)
    plt.xticks(rotation=30, ha='right', fontsize = 12)
    plt.yticks(fontsize = 18)
    plt.savefig('../../results/liss/'+str(helper.iloc[0,1]),format='png', bbox_inches='tight')
    plt.clf()
