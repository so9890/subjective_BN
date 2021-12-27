
"""
Created on Sat Dec 25 14:12:05 2021

@author: sonja

File to read in Liss data
and generate descriptive statistics
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
"read in data"

dtafile= '../../liss_data/qk20a_EN_1.0p.dta'

varfile= '../../liss_data/variable_explanation.xls'
df = pd.read_stata(dtafile)

# some stats
df.tail()
df.info()
df.describe()

# read in list of variable names
var_names=pd.read_excel(varfile, usecols=("A:B"), header=None)
var_names=var_names[~pd.isna(var_names[1])]

#get numbers of variable names which are to be plotted
list_numbers=var_names[0].str[-3:].tolist()

""" visualising data """

#- Bar shart univariate
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
    #plt.xlabel('categories')
    plt.ylabel('percent', fontsize = 18)
    plt.xticks(rotation=30, ha='right', fontsize = 12)
    plt.yticks(fontsize = 18)
    for s in range(0,2):
        if s==1:
            plt.title(str(helper.iloc[0,1]), fontsize=15)
            
        plt.savefig('../../results/liss/'+str(helper.iloc[0,0])+'title'+str(s)+'.png',format='png', bbox_inches='tight')
    plt.clf()

# contingency tables and multivariate barcharts
#   split for each opinion how they behave (I buy second hand,
#   I want to own things, I prefer new products,  I would be open to leasing )

#pre-specify dictionary to store crosstables
dic=dict.fromkeys(['175', '181', '183'])
dic['175']=dict.fromkeys(['135', '141', '144', '147', '148'])
dic['181']=dict.fromkeys(['135', '141', '144', '147', '148'])
dic['183']=dict.fromkeys(['135', '141', '144', '147', '148'])

for i in ['175', '181', '183']:
    dicc=dic[i]
    for j in ['135', '141', '144', '147', '148']:
        dats=df[['qk20a'+i,'qk20a'+j]].copy()
        helper1=var_names[var_names[0]=='qk20a'+str(i)]
        helper2=var_names[var_names[0]=='qk20a'+str(j)]
        
        dats=dats.set_axis([str(helper1.iloc[0,1]), str(helper2.iloc[0,1])], axis=1)
        dicc[j] = pd.crosstab(dats[str(helper1.iloc[0,1])], dats[helper2.iloc[0,1]], normalize='all', margins = False)
        
        # plot bar chart
            
        ax = dicc[j].plot(kind='bar', stacked=True, rot=30)
        ax.legend(title=str(helper2.iloc[0,1]), bbox_to_anchor=(1, 1.02), loc='upper left')
        plt.xticks(rotation=30, ha='right', fontsize = 12)

        plt.savefig('../../results/liss/conditional_'+str(helper1.iloc[0,0])+'_'+str(helper2.iloc[0,0])+'.png', format='png', bbox_inches='tight')
    dic[i]=dicc


""" Intention for behaviour: Second-hand shopping, leasing, buying recycled products

    1) what is driving reducing behaviour? What is preventing it? 
        a) quality, price, environment, social perception, availability
        b) add opinions on environment => is it significant in explaining behaviour?; that is, do people take action?
    2) Who are those people which buy second hand/ lease, buy recycles products? 
        a) Income => motive might be either for too low income or for environmental concerns; correlated through education
        b) Willingness to reduce and skills=> relevant for effect on macroeconomy
    """
    
#-------------------------
"1a) Drivers of reducing behaviour; "
#




# encoding categorical variables
for i in list_numbers:
    df['qk20a'+str(i)+'_cat'] = df['qk20a'+str(i)].cat.codes

# make binary CONTINUE
df.binary_175=0
# Logit
log_reg = smf.logit("qk20a135_cat ~ qk20a175_cat + qk20a181_cat + qk20a183_cat", data=df).fit()
