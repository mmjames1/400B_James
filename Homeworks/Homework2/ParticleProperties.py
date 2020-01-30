#!/usr/bin/env python
# coding: utf-8

# In[23]:


import numpy as np
import astropy.units as u
from ReadFile import Read


# In[93]:


def ParticleInfo(filename,ParticleType,ParticleNumber):
    
    #reading in the return values from ReadFile code (part 2 of HW)
    time,total,data= Read(filename)
    
    #set aside data for a specific particle type, use index for this
    index = np.where(data['type'] == ParticleType)
    
    #getting all of my variables from this new index
    xnew = data['x'][index]
    ynew = data['y'][index]
    znew = data['z'][index]
    Vxnew = data['vx'][index]
    Vynew = data['vy'][index]
    Vznew = data['vz'][index]
    mass_new = data['m'][index]
    
    #read in the x,y, and z position of specific particle
    x = xnew[ParticleNumber]*u.kpc
    y = ynew[ParticleNumber]*u.kpc
    z = znew[ParticleNumber]*u.kpc
    
    #1. get the magnitude of the distance by finding the vector based on all of these componets
    mag_dist = ((x**2)+(y**2)+(z**2))**(1/2)
    #rounding value to 3 decimal places
    mag_dist_3 = np.around(mag_dist,3)
    
    #read in the of the velocity values of the specific particle
    Vx = Vxnew[ParticleNumber]*(u.km/u.s)
    Vy = Vynew[ParticleNumber]*(u.km/u.s)
    Vz = Vznew[ParticleNumber]*(u.km/u.s)
    
    #2. magnitude of velocity in km/s
    mag_vel = ((Vx**2)+(Vy**2)+(Vz**2))**(1/2)
    #rounding value to 3 decimal places
    mag_vel_3 = np.around(mag_vel,3)
    
    #3. Finding the Specific Mass Value
    mass_particle = mass_new[ParticleNumber]*u.Msun*10e10
    #round value to 3 decimal places
    mass_particle_3 = np.around(mass_particle,3)
    
    #converting kpc distance to ly
    mag_dist_ly = mag_dist_3.to(u.lyr)
    #round to 3 decimal places
    mag_dist_ly_3 = np.around(mag_dist_ly,3)
    
    #printing out values needed to answer question 4
    print("The 3D Distance:" , mag_dist_3)
    print("The 3D Velocity:", mag_vel_3)
    print("The Mass:", mass_particle_3)
    print("The Distance in LightYears:",mag_dist_ly_3)
    
    
    return mag_dist_3,mag_vel_3,mass_particle_3


# For part 4 of this HW I'll be calculating the 3D distance, 3D velocity, mass. Additionally the 3D distance will be given in kpc then converted to light years

# In[95]:


#prints out values for 100 particle, use index 99 because particles 
#start at 0
#index 2.0 for disk particle
mag_dist,mag_vel,mass= ParticleInfo("MW_000.txt",2.0,99)


# In[ ]:




