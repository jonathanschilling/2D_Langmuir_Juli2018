#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:59:47 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt

# positions from img_to_coord
posfile = "/home/jonathan/Uni/Forschung/01_Messdaten/Messdaten_2D_Langmuir_Juli2018/2018/0716/005/ximea_QMH_main/positions.dat"

# read positions from manual identification
allImg = []
allX   = []
allY   = []
with open(posfile, "r") as f:
    for l in f.readlines():
        
        lineparts = l.strip().split(" ")
        imgname = lineparts[0]
        xPos = np.double(lineparts[1])
        yPos = np.double(lineparts[2])

        allImg.append(imgname)
        allX.append(xPos)
        allY.append(yPos)

# read positions from piezo positioning unit
posfilePiezo = "/home/jonathan/Uni/Forschung/01_Messdaten/Messdaten_2D_Langmuir_Juli2018/2018/0716/005/2D_Control_Piezo_Linearencoder/path.dat"
posPiezo = np.loadtxt(posfilePiezo)

def scale(xLinEnc, off, scale):
    return np.add(np.multiply(xLinEnc, scale), off)

#%%
plt.figure()
plt.plot(posPiezo[:,0], allX, '+')
plt.plot(posPiezo[:,0], scale(posPiezo[:,0], 28.0, 1.0), '+')
plt.title("X scaling")
plt.xlabel("pos from linear encoders")
plt.ylabel("pos from image")

#%%
plt.figure()
plt.title("Y scaling")
plt.plot(posPiezo[:,1], allY, '+')
plt.plot(posPiezo[:,1], scale(posPiezo[:,1], 54, -1.0), '+')
plt.xlabel("pos from linear encoders")
plt.ylabel("pos from image")

#%%

xOff = 28.5-30.0 # r
xScale = 0.9

yOff = 49 # z
yScale = -0.9

realX = scale(posPiezo[:,0], xOff, xScale)
realY = scale(posPiezo[:,1], yOff, yScale)

plt.figure()
plt.plot(realX, realY, 'rx', label="centers")


gridX = realX.reshape([12,23])
gridY = realY.reshape([12,23])

dx = np.mean(np.diff(gridX, axis=1))
dy = np.mean(np.diff(gridY, axis=0))

# compute corners for pcolormesh
extX = np.hstack([gridX-dx/2, np.transpose([gridX[:,-1]+dx/2])])
cornersX = np.vstack([extX, [extX[-1,:]]])

extY=np.vstack([gridY-dy/2, [gridY[-1,:]+dy/2]])
cornersY=np.hstack([extY, np.transpose([extY[:,-1]])])

plt.plot(cornersX, cornersY, 'b+', label="corners")

plt.xlabel("r / mm")
plt.ylabel("z / mm")
plt.title("measured locations")




