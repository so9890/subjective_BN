""" This file is to 
    1) construct household percentiles of income distribution  
 
    i)   Read in ITBI/ITII files. ITBI until 2004, startind from 2004 ITII, that's 
         when they are first available.
    ii)  Merge sampling weights from FMLI file. 
    -----Tests for how many Consumer Units (CUs) have not reported income ------
    iii) Derive income distribution percentiles. 

"""
##############################################################################

import pandas as pd
import numpy as np

#from functions import  weights_percentiles
##############################################################################

""" i) and ii) Read in data. """

data = pd.read_csv('../original_data/CEX/diary19/expd191.csv')
 
weights=pd.read_csv('../../original_data/CEX_Data/intrvw96/fmli961x.csv')[['NEWID', 'FINLWT21']]
# Note to myself
# the data set weights only contains each CU once, thus, weights are the same for each quarter! Make sense as for each month within a quarter
# the sample is the same , test: 
# unique_weights =pd.unique(weights['NEWID'])
    
data=data.merge(weights, left_on= 'NEWID', right_on= 'NEWID', how= 'left')


""" Tests.

1) check whether all NEWID in data has been matched

 """

if len( data[data['NEWID'].isin(weights['NEWID'])].index)==len(data.index):
    print("Length of data matches. All CUs in data got a sampling weight.")
else:
    print("Error: There are CUs without weight.")
  
    
""" 2) check how many households will be missing due to no income reported."""

if len(pd.unique(data['NEWID']))==len(pd.unique(data[data['UCC']==980000]['NEWID'])):
    print('All households reported income')
else:
    s = len(pd.unique(data['NEWID']))-len(pd.unique(data[data['UCC']==980000]['NEWID']))
    print(s, ' households did not report income and will be missing.' )
   
#################################################################################
    
""" iii) Derive income distribution percentiles for pre-tax income. 
    
Derive percentiles for each month separately. This 
    
"""
 
# Note to myself: the UCC-item 'Income before taxes' and 'Income after taxes' don't 
# need to be divided by 4! a specified for other income variables in the ITBI/ITII files.
      
income_data_before_tax=data[data['UCC']==980000]
    
income_12_1995=income_data_before_tax[income_data_before_tax['REFMO']==12 ]

""" Test. 

Ensure there is only one year, i.e. all observations stem from the same month-year combination.

"""

class MyError(LookupError):
    '''to be looked up.'''


if len(pd.unique(income_12_1995['REFYR']))==1:
    print('test passed: only one year')
else:
   raise MyError('test failed: several years although there should only be one!')
  

""" Derive cummulative distribution function and percentiles. """

d=income_12_1995
d_percentiles = weights_percentiles(d)
d_percentiles_12_1995 = d_percentiles[['NEWID','FINLWT21', 'Percentile']]

""" Save files. """ 
d_percentiles_12_1995.to_pickle('../../data/Percentiles/12_1995')


    