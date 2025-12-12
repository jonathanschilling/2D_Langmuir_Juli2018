# -*- coding: utf-8 -*-
"""
read 2d Langmuir measurement and plot profile of floating potential
"""

import numpy as np
import matplotlib.pyplot as plt

datapath="/home/jonathan/Uni/Forschung/01_Messdaten/Messdaten_2D_Langmuir_Juli2018/"
year="2018/"
monthday="0712/"
series="002/"

positions_subdir="2D_Control_Piezo_Linearencoder/"
probechar_subdir="Langmuir_Kennlinie_SM2400/"

datadir=datapath+year+monthday+series
print("looking for data in "+datadir)

# read positions
posfile=datadir+positions_subdir+"path.dat"
positions=np.round(np.loadtxt(posfile))
num_pos = np.shape(positions)[0]

num_r = len(np.unique(positions[:,0]))
num_z = len(np.unique(positions[:,1]))

r_meshgrid = np.reshape(positions[:,0], [num_z,num_r])
z_meshgrid = np.reshape(positions[:,1], [num_z,num_r])

#plt.plot(positions[:,0], positions[:,1])


# read probe characteristics
phif = np.zeros([num_z, num_r])
i=0
chars=[]
for z in range(num_z):
    for r in range(num_r):    
        
        char_filename=datadir+probechar_subdir+("%04d-00.dat"%i)
        i+=1
        #print(char_filename)
        char_data=np.loadtxt(char_filename)
        chars.append(char_data)
        phif_here=np.interp(0,char_data[:,1],char_data[:,0])
        phif[z,r]=phif_here
        
#%%
        
plt.imshow(np.flipud(phif), interpolation="nearest", extent=[r_meshgrid[0,0], r_meshgrid[0,-1], z_meshgrid[0,0], z_meshgrid[-1,0]])
cb=plt.colorbar()
cb.set_label("floating potential / V")
plt.xlabel("r / mm")
plt.ylabel("z / mm")
plt.tight_layout()
