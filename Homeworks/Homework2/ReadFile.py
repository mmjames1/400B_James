#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Homework 2
#Due: Jan 30 2020
#Mackenzie James

#import modules
import numpy as np
import astropy.units as u


# In[2]:


#(added 2/3)
#Define a function that will read in a file
#: Use: time,total,data= Read('filename')

#takes file name as an input

def Read(filename):
#(added 2/3)
#Input: 
    #filename- for this example MW_000.txt
#Return:
    # time (Myr), total number of particles, and array that holds all the data
    
    
    #open the file
    file = open(filename,'r')
    
    #(added 2/3)
    #read header info line by line (with line as string)
    #read two lines first, store as variable
    
    #Read the first line 
    #(added 2/3) read and store the time
    line1 = file.readline()
    label, value = line1.split()
    #store time in units Myr
    time = float(value)*u.Myr
    
    #read the second line
    #(added 2/3): read and store the number of particles
    line2 = file.readline()
    label, value = line2.split()
    total = float(value)
    
    #close the file
    file.close()
    
    #storing the remainder of the fi;e
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    #filename is input from user
    #(added 2/3): dtype = None specifies data type. None is default Float
    #(added 2/3): default delimiter is line is split using white space
    #skip header skips first 3 lines
    #names=True creates arrays with given labels in file
    
    return time,total,data


# In[3]:


#Added 2/3, did not include in my original submission to check if my code worked
time, total, data = Read("MW_000.txt")


# In[4]:


time
#output 0Myr


# In[5]:


total
#output: 135000.0


# In[6]:


#mass of first particle
data['m'][0]*u.Msun*1e10

#output: 39498500Msun


# In[8]:


#type of first particle
data['type'][0]

#output: 1.0


# In[ ]:




