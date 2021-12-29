
"""
Created on Sat Dec 25 14:12:05 2021

@author: sonja

File to read in Liss data
and generate descriptive statistics
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

from functions_liss import dont

"read in data"

dtafile= '../../liss_data/qk20a_EN_1.0p.dta'
df = pd.read_stata(dtafile)

# replace dont know
df.replace("Don\x92t know", "Don\'t know")
    
# some stats
df.tail()
df.info()
df.describe()

# read in list of variable names
varfile= '../../liss_data/variable_explanation.xls'

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
    
     #replace apostrophe with \' in index
    list_occur=occur.index.tolist()
    for j in range(0,len(occur.index)):
        #strr=list_occur[j]
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
        helper1=var_names[var_names[0]=='qk20a'+str(i)] # opinion => rows
        helper2=var_names[var_names[0]=='qk20a'+str(j)] # action  => columns
        
        dats=dats.set_axis([str(helper1.iloc[0,1]), str(helper2.iloc[0,1])], axis=1)
        # cross table for heatmaps normalize by opinion = index => distribution of action for each opinion
        dicc[j] = pd.crosstab(dats[str(helper1.iloc[0,1])], dats[helper2.iloc[0,1]], normalize='index', margins = False)
        
        #replace apostrophe with \' in index
        dicc[j]=dont(dicc[j])
        
        # plot bar chart => normalise by total number of observations, no margin
        helpps= pd.crosstab(dats[str(helper1.iloc[0,1])], dats[helper2.iloc[0,1]], normalize='all', margins = False)
        helpps=dont(helpps)         # update column and row names using function 
            
        ax = helpps.plot(kind='bar', stacked=True, rot=30)
        ax.legend(title=str(helper2.iloc[0,1]), bbox_to_anchor=(1, 1.02), loc='upper left')
        plt.xticks(rotation=30, ha='right', fontsize = 12)

        plt.savefig('../../results/liss/conditional_bar_'+str(helper1.iloc[0,0])+'_'+str(helper2.iloc[0,0])+'.png', format='png', bbox_inches='tight')
        plt.clf()

        # frequency headmap: joint distribution
        helpps= pd.crosstab(dats[str(helper1.iloc[0,1])], dats[helper2.iloc[0,1]], normalize='all', margins = True)
        helpps=dont(helpps)
        
        data1=helpps.copy()
        data1['All'] = float('nan')   # columns
        data1.loc['All']=float('nan') # rows
        ax = sns.heatmap(data1, annot=True, cmap="BuPu")

        data2=helpps.copy()
        data2.loc[:6,:6]=float('nan')
        sns.heatmap(data2, annot=True, cbar=False, cmap=ListedColormap(['white']))
        
        for p in range(0,2):
            if p==1:
    
                plt.xlabel(str(helper2.iloc[0,1]), fontsize = 15) # x-axis label with fontsize 15
                plt.ylabel(str(helper1.iloc[0,1]), fontsize = 15) # y-axis label with fontsize 15
        
            plt.savefig('../../results/liss/joint_heatmap'+i+'_'+j+'labels'+str(p)+'.png', format='png', bbox_inches='tight')
        plt.clf()
        
    dic[i]=dicc

# heatmaps from cross table with conditional probabilities: distribution of actions withing opinion categories
for i in ['175', '181', '183']:
    s=dic[i]
    for j in ['135', '141', '144', '147', '148']:
        
        # joints
        data1 = s[j].copy()
        #data1[7] = float('nan')
        ax = sns.heatmap(data1.iloc[:6], annot=True, cmap="BuPu")
        
        # Greens
        data2 = s[j].copy()
        #data2.columns = data2.columns.add_categories('All')
        data2['All']=1
        data2.iloc[:,:6] = float('nan')
        sns.heatmap(data2, annot=True, cbar=False, cmap=ListedColormap(['white']))

        for p in range(0,2):
            if p==1:
    
                helper1=var_names[var_names[0]=='qk20a'+str(i)] # opinion => rows
                helper2=var_names[var_names[0]=='qk20a'+str(j)] # action  => columns
                
                plt.xlabel(str(helper2.iloc[0,1]), fontsize = 15) # x-axis label with fontsize 15
                plt.ylabel(str(helper1.iloc[0,1]), fontsize = 15) # y-axis label with fontsize 15

            plt.savefig('../../results/liss/conditional_heatmap'+i+'_'+j+'labels'+str(p)+'.png', format='png', bbox_inches='tight')
        plt.clf()
        
        
# bubble plot using frequency as size indicator

#he=df.groupby(['qk20a'+str(i), 'qk20a'+str(j)]).size().reset_index()
#plt.scatter(he['qk20a'+str(i)], he['qk20a'+str(j)], s=he[0])

""" Intention for behaviour: Second-hand shopping, leasing, buying recycled products

    1) what is driving reducing behaviour? What is preventing it? 
        a) quality, price, environment, social perception, availability
        b) add opinions on environment => is it significant in explaining behaviour?; that is, do people take action?
In next steps    2) Who are those people which buy second hand/ lease, buy recycles products? 
        a) Income => motive might be either for too low income or for environmental concerns; correlated through education
        b) Willingness to reduce and skills=> relevant for effect on macroeconomy
    """
    
#-------------------------
"1a) Drivers of reducing behaviour"
#

# encoding categorical variables
#for i in list_numbers:
 #   df['qk20a'+str(i)+'_cat'] = df['qk20a'+str(i)].cat.codes

#plt.scatter(df.qk20a175_cat, df.qk20a135_cat)
#s=df[['qk20a175_cat', 'qk20a135_cat']].corr().style.background_gradient(cmap="Blues")

# make binary CONTINUE
#df.binary_175=0
# Logit
#log_reg = smf.logit("qk20a135_cat ~ qk20a175_cat + qk20a181_cat + qk20a183_cat", data=df).fit()
