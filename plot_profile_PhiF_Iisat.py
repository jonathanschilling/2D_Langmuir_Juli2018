# -*- coding: utf-8 -*-
"""
read 2d Langmuir measurement and plot profile of floating potential
"""

import numpy as np
import matplotlib.pyplot as plt

datapath="/home/jonathan/Uni/Forschung/01_Messdaten/Messdaten_2D_Langmuir_Juli2018/"
#datapath="/run/media/jonathan/Stick16G/Messdaten_2D_Langmuir/"
year="2018/"
monthday="0716/"
series="002/"

positions_subdir="2D_Control_Piezo_Linearencoder/"
probechar_subdir="Langmuir_Kennlinie_SM2400/"

datadir=datapath+year+monthday+series
vacuum_meas_filename = datapath+year+monthday+"500mT_no_plasma.dat"


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

vacuum_char = np.loadtxt(vacuum_meas_filename)

# read probe characteristics
phif = np.zeros([num_z, num_r])
iisat = np.zeros([num_z, num_r])
i=0
chars=[]
for z in range(num_z):
    for r in range(num_r):    
        
        char_filename=datadir+probechar_subdir+("%04d-00.dat"%i)
        i+=1
        #print(char_filename)
        char_data=np.loadtxt(char_filename)
        
        char_data[:,1] -= vacuum_char[:,1]
        
        # add 99.6V to voltage (HP 6516A power supply in series to SourceMeter)
#        char_data[:,0]+=99.6
        
        chars.append(char_data)
        phif_here=np.interp(0,char_data[:,1],char_data[:,0])
        phif[z,r]=phif_here
        
        iisat_here=np.interp(phif_here-30.0, char_data[:,0], char_data[:,1])
        iisat[z,r]=iisat_here
        
#%%
   
z0=0


def onclick(event):
#    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
#          ('double' if event.dblclick else 'single', event.button,
#           event.x, event.y, event.xdata, event.ydata))
    r_idx=np.int(np.round(event.xdata)/2.0)
    z_idx=np.int((z_meshgrid[0,0]-np.round(event.ydata-1.0))/2.0)
    #print(r_idx)
    #print(z_idx)
    print("r=%g z=%g"%(r_meshgrid[z_idx,r_idx],z_meshgrid[z_idx,r_idx]))
    char_idx=z_idx*num_r+r_idx
    
    fig42 = plt.figure(42)
    plt.plot(chars[char_idx][:,0], chars[char_idx][:,1]*1e6, label="r=%g z=%g"%(r_meshgrid[z_idx,r_idx],z_meshgrid[z_idx,r_idx]))
#    plt.plot(chars[char_idx][:,0], chars[char_idx][:,2], label="dI/dU")
    plt.xlabel("U / V")
    plt.ylabel("I / uA")
    plt.legend(loc="best")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    

fig=plt.figure()     
plt.imshow(np.flipud(phif[z0:,:]), interpolation="nearest", extent=[r_meshgrid[z0,0]-1, r_meshgrid[z0,-1]+1, z_meshgrid[z0,0]+1, z_meshgrid[-1,0]-1])
cb=plt.colorbar()
cb.set_label("floating potential / V")
plt.xlabel("r / mm")
plt.ylabel("z / mm")
plt.tight_layout()
cid = fig.canvas.mpl_connect('button_press_event', onclick)

#%%
   
fig2=plt.figure()     
plt.imshow(np.flipud(iisat[z0:,:]*-1e6), interpolation="nearest", extent=[r_meshgrid[z0,0], r_meshgrid[z0,-1], z_meshgrid[z0,0], z_meshgrid[-1,0]])
cb=plt.colorbar()
cb.set_label("ion saturation current (at 30V below phif) / uA")
plt.xlabel("r / mm")
plt.ylabel("z / mm")
plt.tight_layout()
cid = fig2.canvas.mpl_connect('button_press_event', onclick)
