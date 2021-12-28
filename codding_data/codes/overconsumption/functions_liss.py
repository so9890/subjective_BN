#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 15:03:54 2021

@author: sonja
"""
def dont(helpps):
    """ replace column and row names if contains dont

    """

    list_index=helpps.index.tolist()
    
    for t in range(0,len(helpps.index)):
        list_index[t]=list_index[t].replace("Don\x92t know", "Don\'t know")
    helpps.index=list_index
    
    #replace apostrophe with \' in column names
    list_columns=helpps.columns.tolist()

    for t in range(0,len(helpps.columns)):
        list_columns[t]=list_columns[t].replace("Don\x92t know", "Don\'t know")
    helpps.columns=list_columns
    
        
    return helpps
