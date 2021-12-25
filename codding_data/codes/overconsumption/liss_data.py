
"""
Created on Sat Dec 25 14:12:05 2021

@author: sonja

File to read in Liss data
"""

import pandas as pd
import numpy as np

"read in data"

"dtafile = '../../liss_data/bf17e_EN_1.0p.dta'"

dtafile= '../../liss_data/ad08a_1p_EN.dta'

df = pd.read_stata(dtafile)
df.tail()
