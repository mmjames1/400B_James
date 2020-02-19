#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Homework 4 - COM Calculations
#Mackenzie James
#due: 2/13

#importing modules needed
from ReadFile import Read
import numpy as np
import astropy.units as u
import astropy.table as tbl


# In[2]:


class CenterOfMass:
#creates a class to define the center of mass position and velocity properties

    #initializing a class so each objected will store data from the simulation-based on particle type
    def __init__(self,filename,particletype):
    #inputs:
        #self
        #filename- for this we are using MW_000.txt, M31_000.txt, and M33_000.txt
        #particle type: 1.0=halo, 2.0=disk, 3.0=bulge
    #returns:
        #sets indices
    
        #read data from user given file (detailed in filename above)
        self.time,self.total,self.data=Read(filename)
    
        #creates an array with the user specified particle type
        self.index = np.where(self.data['type']==particletype)
    
        #storing the mass, position, and velocities of a particle 
        #Same idea as in HW2 ParticleProperties- this time just defining with self 
    
        #mass
        self.m = self.data['m'][self.index]
        #x,y,z position
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        #x,y,z velocity
        self.Vx = self.data['vx'][self.index]
        self.Vy = self.data['vy'][self.index]
        self.Vz = self.data['vz'][self.index]
    

    #function that generically returns the 3D cooridinates of the center of mass
    #(position or velocity ) of any galaxy
    def COMdefine(self,a,b,c,m):
    #inputs: 
        #a,b,c: correspond to the x,y,z position or velocity coordinates of a particle
        #m: mass of given particle
        #self
    #Returns: the x,y,z position or velocity coordinates of the center of mass


        #Calculating x component of center of mass vector
        Xcom = (np.sum(a*m))/(np.sum(m))
        #Calculating y component of center of mass vector
        Ycom = (np.sum(b*m))/(np.sum(m))
        #Calculating z component of center of mass vector
        Zcom = (np.sum(c*m))/(np.sum(m))

        return Xcom,Ycom,Zcom

    #function that will determine the center of mass of position vectors
    #of a user given galaxy and particle type
    def COM_P(self,delta):
    #inputs: 
        #self
        #delta: the dolerance that decides whether the center of mass position has convered
    #returns:
        # position vector for the center of mass (units: kpc)

        ##step 1: first guess COM position##
        
        Xcomp,Ycomp,Zcomp = self.COMdefine(self.x, self.y, self.z, self.m)
        
        #computing mag of the vector
        Rcom = np.sqrt((Xcomp**2)+(Ycomp**2)+(Zcomp**2))
        
        #changing particle reference frame to COM frame
        #computing difference between vecors
        Xnew = self.x-Xcomp
        Ynew = self.y-Ycomp
        Znew = self.z-Zcomp
        #computing the new magnitude of the vector from new components
        Rnew = np.sqrt((Xnew**2)+(Ynew**2)+(Znew**2))
        
        #finding the max 3D seperation
        Rmax = max(Rnew)/2.0 #restarting at half the radius
        
        ##step 2: refining guess##
        
        #pick initial value for change in COM position - want to be larger, using suggested 1000kpc value
        change=1000.0 #units: kpc
        
        #starting a while loop to determine the center of mass position
        #this is where we introduce tolerance term
        #delta is the tolerance for the difference in the old and new center of mass
        #the while loop continues while the change in Rcom is larger than delta
        
        while (change > delta):
            #selects all the particles in the reduced radius, starting from original 
            index2 = np.where(Rnew <= Rmax)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]  
            
            #refine COM posotion
            #compute the center of mass using particles in the new reduced radius
            Xcom2,Ycom2,Zcom2 = self.COMdefine(x2,y2,z2,m2)
            #computing new 3d COM position
            Rcom2 = np.sqrt((Xcom2**2)+(Ycom2**2)+(Zcom2**2))
            
            #determine the difference between the center of mass positions
            change = abs(Rcom-Rcom2)
            
            #reset the Rmax value before looping again (divide volume by factor of 2 again)
            Rmax = Rmax/2.0
            
            #change the frame of reference to the new COM calcuated
            
            ###Correction 2/29: changing Xnew=x2-Xcom2 to Xnew=self.x-Xcom2
            ###from Prof. Besla's corrections/solutions. This will now select the 
            ###correct particles
            
            Xnew = self.x-Xcom2
            Ynew = self.y-Ycom2
            Znew = self.z-Zcom2
            #computing the new magnitude of the vector from new components
            Rnew = np.sqrt((Xnew**2)+(Ynew**2)+(Znew**2))
            
            #setting up center of mass to the new refined values
            Xcomp = Xcom2
            Ycomp = Ycom2
            Zcomp = Zcom2
            Rcom = Rcom2
            
            #creating vector to store the center of mass position
            COM_P = [Xcomp,Ycomp,Zcomp]
            #putting vector into correct units, kpc   
            COM_P_units = COM_P*u.kpc
            #roudning values to 2 decimal places
            COM_P_round = np.round(COM_P_units,2)
            
        return COM_P_round
            
       
    #function that will determine the center of mass of velocity vectors
    def COM_V(self,comX,comY,comZ):
    #inputs:
        #comX,comY,comZ are the x,y,z positions of the center of mass (units:kpc)
        #self
    #returns:
        #Vx,Vy,Vz for center of mass (units: km/s)
        
        #max distance from the center that we store all the particles for
        RVmax = 15.0*u.kpc
        
        #position of all particles relative to COM position
        xV = (self.x*u.kpc)-comX
        yV = (self.y*u.kpc)-comY
        zV = (self.z*u.kpc)-comZ
        
        #calculating vector
        RV = np.sqrt((xV**2)+(yV**2)+(zV**2))
        
        #determine new index 
        indexV = np.where(RV <= RVmax)
        
        #velocity and mass of particles in the max radius
        Vxnew = self.Vx[indexV]
        Vynew = self.Vy[indexV]
        Vznew = self.Vz[indexV]
        mnew = self.m[indexV]
        
        #compute the COM velocity using the components above
        Vxcom, Vycom, Vzcom = self.COMdefine(Vxnew,Vynew,Vznew,mnew)
        
        #create the vector to store COM velocity
        COM_V = [Vxcom,Vycom,Vzcom]
        #putting velocity in correct units: km/s
        COM_V_units = COM_V*(u.km/u.s)
        #round component values to 2 decimal places   
        COM_V_round = np.round(COM_V_units,2) 
    
        return COM_V_round
    


# ## Question 6, Testing the Code.

# In[3]:


#finding the COM Position and Velocity Vector for MW, M31, and M33

##Milky Way
MWCOM = CenterOfMass("MW_000.txt", 2)
MW_COMP = MWCOM.COM_P(0.1)
MW_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2])
print("Milky Way COM Position:",MW_COMP)
print("Milky Way COM Velocity:",MW_COMV)

##M31
M31COM = CenterOfMass("M31_000.txt", 2)
M31_COMP = M31COM.COM_P(0.1)
M31_COMV = M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2])
print("M31 COM Position:",M31_COMP)
print("M31 COM Velocity:",M31_COMV)

##M33
M33COM = CenterOfMass("M33_000.txt", 2)
M33_COMP = M33COM.COM_P(0.1)
M33_COMV = M33COM.COM_V(M33_COMP[0], M33_COMP[1], M33_COMP[2])
print("M33 COM Position:",M33_COMP)
print("M33 COM Velocity:",M33_COMV)


# In[4]:


#What is the magnitude of current separation and velocity of MW and M31?

#in order to find the seperation vector, take the difference between the x,y,and z components
Xdif = MW_COMP[0]-M31_COMP[0]
Ydif = MW_COMP[1]-M31_COMP[1]
Zdif = MW_COMP[2]-M31_COMP[2]

#magnitude of position separation
Pos_sep = np.sqrt((Xdif**2)+(Ydif**2)+(Zdif**2))
print("current separation between the Milky way and M31 is:",Pos_sep)

Vxdif = MW_COMV[0]-M31_COMV[0]
Vydif = MW_COMV[1]-M31_COMV[1]
Vzdif = MW_COMV[2]-M31_COMV[2]
velocity = np.sqrt((Vxdif**2)+(Vydif**2)+(Vzdif**2))
print("current velocity difference between Milky way and M31 is:",velocity)


# Code output: 
# current separation between the Milky way and M31 is: 769.0 kpc
# current velocity difference between Milky way and M31 is: 118.0 km/ s

# In[87]:


#what is the magnitude of current separation and velocity of M31 and M33?

#this will be the same process as before, just referencing different galaxies
Xdif2 = M31_COMP[0]-M33_COMP[0]
Ydif2 = M31_COMP[1]-M33_COMP[1]
Zdif2 = M31_COMP[2]-M33_COMP[2]

#magnitude of position separation
Pos_sep2 = np.sqrt((Xdif2**2)+(Ydif2**2)+(Zdif2**2))
print("current separation between the M31 and M33 is:",Pos_sep2)

Vxdif2 = M31_COMV[0]-M33_COMV[0]
Vydif2 = M31_COMV[1]-M33_COMV[1]
Vzdif2 = M31_COMV[2]-M33_COMV[2]
velocity2 = np.sqrt((Vxdif2**2)+(Vydif2**2)+(Vzdif2**2))
print("current velocity difference between M31 and M33 is:",velocity2)


# Code output: 
# current separation between the M31 and M33 is: 201.0 kpc
# 
# current velocity difference between M31 and M33 is: 199.0 km / s

# Given that M31 and MW are about to merge, why is the iterative process to determine the COM important?
# 
# This is important because we should have a very accurate measurement of where the center of mass is in the galaxy. When the galaxies merge, the movement of the collision is going to be centered around the interaction between those two COM. If we know how each center of mass is being affected, we can more accurately model what the collision for MW and M31 will look like

# In[ ]:




