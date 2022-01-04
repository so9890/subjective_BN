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
    
"""Dynamic:
    1) create indicator to before and after reduction on individual level
        How to treat those which increase hours again?
        Need to treat those which work more hours again differently
        Include hours worked dynamically as outcome variable and voluntary as a regressor"""
    



    

""" What explains the probability to think new clothes/furniture are not necessary?

    1. step: economic/ demographics observables:
    regressors: income, education, family size, change in family size
    
    2. step: errors on environmental concerns (not sure this works with categorical variable)
    """
