#!/usr/bin/env python
# coding: utf-8

# In[37]:


import numpy as np
import astropy.units as u


# In[44]:


#takes file name as an input
def Read(filename):
    #open the file
    file = open(filename,'r')
    #Read the first line
    line1 = file.readline()
    label, value = line1.split()
    #store time in units Myr
    time = float(value)*u.Myr
    
    #read the second line
    line2 = file.readline()
    label, value = line2.split()
    total = float(value)
    
    #close the file
    file.close()
    
    #storing the remainder of the file
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    #filename is input from user
    #dtype means line is split using space
    #skip header skips first 3 lines
    #names=True creates arrays with given labels in file
    
    return time,total,data

