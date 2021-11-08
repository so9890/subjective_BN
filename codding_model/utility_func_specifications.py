# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

# functions
def utility(c, B):
    """ Calculate the percentile of each household.
    
    1) Sum up weights grouped by income value.
    2) Derive cummulative distribution and probibility distribution functions.
    3) Assign percentiles to each household.
            
    d is the data set

    """
    
    utils = -(c-B)**2
    #utils = np.log(c)-B*c**2
    #utils = c**(1/3)-B*c
    
    return utils


# plot utility
B=4 # if B big enough there is a satiation point! 
c= np.linspace(0,10,100)

u=utility(c, B)

# plot figures
fig =plt.figure()

plt.plot(c,u, 'r')
umax=np.amax(u)
i=np.where(u== np.amax(u))
