# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

# functions
def utility(c, B):
  
    utils = -(c-B)**2
    #utils = np.log(c)-B*c**2
    #utils = c**(1/3)-B*c
    
    return utils


def habits(c, etta, thetta):
  
    utils = (c-thetta)**(1-1/etta)
    
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
