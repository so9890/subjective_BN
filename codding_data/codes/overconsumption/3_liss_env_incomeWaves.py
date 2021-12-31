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

"""read in data and merge to environmental data set """

# environmental data 
dtafile= '../../liss_data/environment/qk20a_EN_1.0p.dta'
df = pd.read_stata(dtafile)
