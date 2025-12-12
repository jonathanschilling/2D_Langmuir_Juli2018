# -*- coding: utf-8 -*-
"""
read 2d Langmuir measurement and plot profile of floating potential
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as s


datapath="/home/jonathan/Uni/Forschung/01_Messdaten/Messdaten_2D_Langmuir_Juli2018/"
#datapath="/run/media/jonathan/Stick16G/Messdaten_2D_Langmuir/"
year="2018/"
monthday="0716/"
series="001/"

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


def zero_crossing(xp, fp):
    last_before=np.where(fp>0)[0][-1]
    first_after=np.where(fp<0)[0][0]
    return np.interp(0, fp[last_before:first_after], xp[last_before:first_after])

#Ioff=-2e-6
Ioff=0.0

# read probe characteristics
phif = np.zeros([num_z, num_r])
iisat = np.zeros([num_z, num_r])
phip = np.zeros([num_z, num_r])
i=0
chars=[]
for z in range(num_z):
    for r in range(num_r):    
        
        char_filename=datadir+probechar_subdir+("%04d-00.dat"%i)
        i+=1
#        print(char_filename)
        char_data=np.loadtxt(char_filename)
        
        char_data[:,1] -= vacuum_char[:,1]+Ioff
        
        # add 99.6V to voltage (HP 6516A power supply in series to SourceMeter)
#        char_data[:,0]+=99.6
        
        chars.append(char_data)
        phif_here=np.interp(0,char_data[:,1],char_data[:,0])
        phif[z,r]=phif_here
        
        iisat_here=np.interp(phif_here-30.0, char_data[:,0], char_data[:,1])
        iisat[z,r]=iisat_here
        
        U=char_data[:,0]
        I=char_data[:,1]
        
        # calc derivate
        window_length = 11
        poly_order=2
        smoothed_I = s.savgol_filter(I, window_length, poly_order)
        dI_dU = s.savgol_filter(I, window_length, poly_order, deriv=1)
        d2I_dU2 = s.savgol_filter(I, window_length, poly_order, deriv=2)
        
        # find plasma potential: look for first zero of 2nd derivative after floating potential
#        phif=np.interp(0,smoothed_I,U)
        first_idx_after = np.where(U>phif_here)[0][0]
        
        # find relative maxima
        relmax = s.argrelextrema(dI_dU, np.greater)
        first_relmax=np.where(relmax[0]>first_idx_after)[0][0]
        phip_here=U[relmax[0][first_relmax]]
        
        phip[z,r]=phip_here
        
#%%
   
z0=0


def plot_char(r_idx, z_idx):
    char_idx=z_idx*num_r+r_idx
    print("r=%g z=%g ==> %04d-00.dat"%(r_meshgrid[z_idx,r_idx],z_meshgrid[z_idx,r_idx],char_idx))
   
        
    
    
    U=chars[char_idx][:,0]
    I=chars[char_idx][:,1]
    
    # calc derivate
    window_length = 11
    poly_order=2
    smoothed_I = s.savgol_filter(I, window_length, poly_order)
    dI_dU = s.savgol_filter(I, window_length, poly_order, deriv=1)
    d2I_dU2 = s.savgol_filter(I, window_length, poly_order, deriv=2)
    
    # find plasma potential: look for first zero of 2nd derivative after floating potential
    phif_here=np.interp(0,smoothed_I,U)
    first_idx_after = np.where(U>phif_here)[0][0]
    
    # find relative maxima
    relmax = s.argrelextrema(dI_dU, np.greater)
    first_relmax=np.where(relmax[0]>first_idx_after)[0][0]
    phip_here=U[relmax[0][first_relmax]]
        
        
    fig42=plt.figure(42)
    plt.plot(U, I*1e6, label="raw char")
    plt.plot(U, smoothed_I*1e6, label="smoothed char")
    plt.xlabel("U / V")
    plt.ylabel("I / uA")
    plt.legend(loc="best")
    plt.grid(True)
    plt.axvline(phif_here, color="r")
    plt.axvline(phip_here, color="b")
    plt.tight_layout()
    
    fig43=plt.figure(43)
    plt.plot(U, dI_dU*1e6, label="derivative")
    plt.xlabel("U / V")
    plt.ylabel("dI/dU / uA/V")
    plt.legend(loc="best")
    plt.grid(True)
    plt.axvline(phif_here, color="r")
    plt.axvline(phip_here, color="b")
    plt.tight_layout()
    
    fig44=plt.figure(44)
    plt.plot(U, d2I_dU2*1e6, label="2nd derivative")
    plt.xlabel("U / V")
    plt.ylabel("d^2I/dU^2 / uA/V^2")
    plt.legend(loc="best")
    plt.grid(True)
    plt.axvline(phif_here, color="r")
    plt.axvline(phip_here, color="b")
    plt.tight_layout()

    
    
#    fig42 = plt.figure(42)
#    plt.plot(chars[char_idx][:,0], chars[char_idx][:,1]*1e6, label="r=%g z=%g"%(r_meshgrid[z_idx,r_idx],z_meshgrid[z_idx,r_idx]))
##    plt.plot(chars[char_idx][:,0], chars[char_idx][:,2], label="dI/dU")
#    plt.xlabel("U / V")
#    plt.ylabel("I / uA")
#    plt.legend(loc="best")
#    plt.grid(True)
#    plt.tight_layout()
#    plt.show()

    



def onclick_profile(event):
#    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
#          ('double' if event.dblclick else 'single', event.button,
#           event.x, event.y, event.xdata, event.ydata))
    r_idx=np.int(np.round(event.xdata)/2.0)
    z_idx=np.int((z_meshgrid[0,0]-np.round(event.ydata-1.0))/2.0)
    #print(r_idx)
    #print(z_idx)
    
    plot_char(r_idx, z_idx)
    

fig=plt.figure()     
plt.imshow(np.flipud(phif[z0:,:]), interpolation="nearest", extent=[r_meshgrid[z0,0]-1, r_meshgrid[z0,-1]+1, z_meshgrid[z0,0]+1, z_meshgrid[-1,0]-1], vmin=0)
cb=plt.colorbar()
cb.set_label("floating potential / V")
plt.xlabel("r / mm")
plt.ylabel("z / mm")
plt.tight_layout()
cid = fig.canvas.mpl_connect('button_press_event', onclick_profile)

#%%
   
fig2=plt.figure()     
plt.imshow(np.flipud(iisat[z0:,:]*-1e6), interpolation="nearest", extent=[r_meshgrid[z0,0], r_meshgrid[z0,-1], z_meshgrid[z0,0], z_meshgrid[-1,0]])
cb=plt.colorbar()
cb.set_label("ion saturation current (at 30V below phif) / uA")
plt.xlabel("r / mm")
plt.ylabel("z / mm")
plt.tight_layout()
cid = fig2.canvas.mpl_connect('button_press_event', onclick_profile)


#%%
   
fig3=plt.figure()     
plt.imshow(np.flipud(phip[z0:,:]), interpolation="nearest", extent=[r_meshgrid[z0,0], r_meshgrid[z0,-1], z_meshgrid[z0,0], z_meshgrid[-1,0]], vmin=0)
cb=plt.colorbar()
cb.set_label("plasma potential / V")
plt.xlabel("r / mm")
plt.ylabel("z / mm")
plt.tight_layout()
cid = fig3.canvas.mpl_connect('button_press_event', onclick_profile)


#%%
#
char_idx=116

U=chars[char_idx][:,0]
I=chars[char_idx][:,1]

# calc derivate
window_length = 11
poly_order=2
smoothed_I = s.savgol_filter(I, window_length, poly_order)
dI_dU = s.savgol_filter(I, window_length, poly_order, deriv=1)
d2I_dU2 = s.savgol_filter(I, window_length, poly_order, deriv=2)

# find plasma potential: look for first zero of 2nd derivative after floating potential
phif_here=np.interp(0,smoothed_I,U)
first_idx_after = np.where(U>phif_here)[0][0]

# find relative maxima
relmax = s.argrelextrema(dI_dU, np.greater)
first_relmax=np.where(relmax[0]>first_idx_after)[0][0]
phip_here=U[relmax[0][first_relmax]]
#
##plasma_pot_idx = np.where(d2I_dU2[first_idx_after:]<0)[0][0]+first_idx_after
#plasma_pot_idx = np.argmax(dI_dU)
##phip_here=np.interp(0, d2I_dU2[first_idx_after:plasma_pot_idx+1], U[first_idx_after:plasma_pot_idx+1])
#phip_here=zero_crossing(U[plasma_pot_idx-10:plasma_pot_idx+10], d2I_dU2[plasma_pot_idx-10:plasma_pot_idx+10])

fig42=plt.figure(42)
plt.plot(U, I*1e6, label="raw char")
plt.plot(U, smoothed_I*1e6, label="smoothed char")
plt.xlabel("U / V")
plt.ylabel("I / uA")
plt.legend(loc="best")
plt.grid(True)
plt.axvline(phif_here, color="r")
plt.axvline(phip_here, color="b")
plt.tight_layout()

fig43=plt.figure(43)
plt.plot(U, dI_dU*1e6, label="derivative")
plt.xlabel("U / V")
plt.ylabel("dI/dU / uA/V")
plt.legend(loc="best")
plt.grid(True)
plt.axvline(phif_here, color="r")
plt.axvline(phip_here, color="b")
plt.tight_layout()

fig44=plt.figure(44)
plt.plot(U, d2I_dU2*1e6, label="2nd derivative")
plt.xlabel("U / V")
plt.ylabel("d^2I/dU^2 / uA/V^2")
plt.legend(loc="best")
plt.grid(True)
plt.axvline(phif_here, color="r")
plt.axvline(phip_here, color="b")
plt.tight_layout()