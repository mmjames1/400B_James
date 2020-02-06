#!/usr/bin/env python
# coding: utf-8

# In[49]:


#Homework 3
#Due: 2/6/2020
#Mackenzie James

#calling modules needed
import numpy as np
from ReadFile import Read
import astropy.units as u


# In[89]:


#creating the program ComponentMass
#this program will return the total mass of any desired galaxy component

def ComponentMass(filename,ParticleType):
#inputs:
    #filename: for this Hw it is either MW_000.txt,M31_000.txt, or M33_000.txt
    #particle type: 1.0 for Halo, 2.0 for Disk, or 3.0 for Bulge
#return:
    #Mass: units of 10e12Msun, rounded to 3 decimal places

#reading in the return values from ReadFile code (part 2 of HW)
    time,total,data= Read(filename)
    
    #this creates an index based on the particle type input
    index = np.where(data['type'] == ParticleType)
    
    #creates an array with all of the mass values from a specific particle type input
    mass= data['m'][index]
    
    mass2 = mass * 1e-2
    #this is to correct for the units because the given mass value is in 1e10
    
    #To find the total mass of a specific type of particle, using the sum
    #function to add up every mass value created from the array
    TotalMass = np.around((sum(mass2)*1e12*u.Msun),3)
    #multiplying constant at the end to put in units of 10e12Msun
    
    return TotalMass


# Below is all of the code that I used to generate the table of values for masses. I calculated the total mass of each component of each galaxy, made a table of those using the pandas dataframe package for python, and then added on to that main table as I went. The final table was created in LaTeX and is attached seperately. 

# In[90]:


ComponentMass("MW_000.txt",1.0)


# In[79]:


ComponentMass("MW_000.txt",2.0)


# In[80]:


ComponentMass("MW_000.txt",3.0)


# In[81]:


ComponentMass("M31_000.txt",1.0)


# In[82]:


ComponentMass("M31_000.txt",2.0)


# In[83]:


ComponentMass("M31_000.txt",3.0)


# In[84]:


ComponentMass("M33_000.txt",1.0)


# In[85]:


ComponentMass("M33_000.txt",2.0)


# In[86]:


ComponentMass("M33_000.txt",3.0)


# Part 3 of the hw, storing the results in a table

# In[210]:


import pandas as pd


# In[211]:


"""Units note, want this all in 1e12Msun"""
#creating the arrays for Galaxy Name
Names = ['Milky Way','M31','M33']
#Halo Mass Column
HaloMass = [1.975,1.921,0.188]
#Disk Mass Column
DiskMass = [0.07500,0.1200,0.009300]
#Bulge Mass Column
BulgeMass = [0.01001,0.01905,0.000]

#this will create the dataframe for my inital table to work with
d = {'Galaxy Name': Names,'Halo Mass (1e12Msun) ': HaloMass, 'Disk Mass (1e12Msun)': DiskMass, 'Bulge Mass (1e12Msun)': BulgeMass}
df = pd.DataFrame(data=d)
df.set_index('Galaxy Name', inplace=True)


# In[212]:


df


# In[213]:


#by summing up all of the rows, I will get a total mass for each galaxy
df.sum(axis=1)


# In[214]:


#creating the array for the total mass of each galaxy 
Total = [2.060,2.060,0.197]
#adding this to my original dataframe so now the table will show the total values as well
dTotal = {'Galaxy Name': Names,'Halo Mass (1e12Msun) ': HaloMass, 'Disk Mass (1e12Msun)': DiskMass, 'Bulge Mass (1e12Msun)': BulgeMass,'Total (1e12Msun)':Total}
dfTotal = pd.DataFrame(data=dTotal)
dfTotal.set_index('Galaxy Name', inplace=True)


# In[215]:


dTotal = {'Galaxy Name': Names,'Halo Mass (1e12Msun) ': HaloMass, 'Disk Mass (1e12Msun)': DiskMass, 'Bulge Mass (1e12Msun)': BulgeMass,'Total (1e12Msun)':Total}
dfTotal = pd.DataFrame(data=dTotal)
dfTotal.set_index('Galaxy Name', inplace=True)


# In[216]:


dfTotal


# In[217]:


#By summing the columns, I will get a value for each mass component for the Local Group
dfTotal.sum(axis=0)


# In[218]:


#Because I am lengthing the columns, I am redefining my columns for the dataframe by adding the values I calculated 
#for the local group at the end of each array

#creating the arrays for Galaxy Name
Names = ['Milky Way','M31','M33','Local Group']
#Halo Mass Column
HaloMass = [1.975,1.921,0.188,4.084]
#Disk Mass Column
DiskMass = [0.07500,0.1200,0.009300,0.204]
#Bulge Mass Column
BulgeMass = [0.01001,0.01905,0.000,0.0291]
#Total Mass Column
Total = [2.060,2.060,0.197,4.317]

#redefining my dataframe to include the calculated values for the local group
dTotal = {'Galaxy Name': Names,'Halo Mass (1e12Msun) ': HaloMass, 'Disk Mass (1e12Msun)': DiskMass, 'Bulge Mass (1e12Msun)': BulgeMass,'Total (1e12Msun)':Total}
dfTotal = pd.DataFrame(data=dTotal)
dfTotal.set_index('Galaxy Name', inplace=True)


# In[219]:


dfTotal


# In[221]:


#Doing the last addition to my table by calculating the values for the fbar column
#fbar = disk+bulge / total

#calculating fbar for each row by using the previously created arrays for mass componets (instead of hard coding with numbers)
MWbary = (DiskMass[0] + BulgeMass[0])/Total[0]
M31bary = (DiskMass[1] + BulgeMass[1])/Total[1]
M33bary = (DiskMass[2] + BulgeMass[2])/Total[2]
Totalbary = (DiskMass[3] + BulgeMass[3])/Total[3]


# In[220]:


#creating the array to add to my table
FBar=np.round([MWbary,M31bary,M33bary,Totalbary],3)


# In[224]:


#adding the fBar column, editing my dataframe for the last time, this will create my final table with all the 
#needed calculated values for the galaxies in the local group
#this final table is also created in LaTeX and is attached in a seperate pdf

dFinal = {'Galaxy Name': Names,'Halo Mass (1e12Msun) ': HaloMass, 'Disk Mass (1e12Msun)': DiskMass, 'Bulge Mass (1e12Msun)': BulgeMass,'Total (1e12Msun)':Total, 'fbar':FBar}
dfFinal = pd.DataFrame(data=dFinal)
dfFinal.set_index('Galaxy Name', inplace=True)


# In[225]:


dfFinal

