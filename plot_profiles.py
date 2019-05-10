# -*- coding: utf-8 -*-
#
# read 2d Langmuir measurement and plot profile of floating potential
# 
# Jonathan Schilling (jonathan.schilling@ipp.mpg.de)


import numpy as np
import matplotlib.pyplot as plt

listOfAnalysisParameters = []

# common parameters
analysisParameters = {}
analysisParameters["datapath"]="/home/jonathan/Uni/Forschung/01_Messdaten/Messdaten_2D_Langmuir_Juli2018"
analysisParameters["outDir"] = "/home/jonathan/Uni/Forschung/03_Auswertung/2D_Langmuir_Juli2018/out"
analysisParameters["cmap"] = "jet"
analysisParameters["contourLevels"] = 20
analysisParameters["minPhiF"] = 45 # V
analysisParameters["maxPhiF"] = 80 # V
analysisParameters["minIisat"] = 0 # uA
analysisParameters["maxIisat"] = 20 # uA
analysisParameters["minR"] = -5 # mm
analysisParameters["maxR"] = 40 # mm
analysisParameters["minZ"] =  0 # mm
analysisParameters["maxZ"] = 30 # mm
analysisParameters["year"]="2018"

# build individual runs
analysisParameters["monthday"]="0716"
analysisParameters["series"]="002"
analysisParameters["magFieldLabel"] = "500mT"
listOfAnalysisParameters.append(analysisParameters)





########## execution below ##########

def getOutputFilePrefix(analysisParameters):
    year          = analysisParameters["year"]
    monthday      = analysisParameters["monthday"]
    series        = analysisParameters["series"]
    magFieldLabel = analysisParameters["magFieldLabel"]
    my_cmap       = analysisParameters["cmap"]
    contourLevels = analysisParameters["contourLevels"]
    return magFieldLabel+"_"+year+"_"+monthday+"_"+series+"_"+my_cmap+"_"+str(contourLevels)

def plotElectrodes(plt):
    
    
    

           
def runAnalysis(listOfAnalysisParameters):
    positions_subdir="2D_Control_Piezo_Linearencoder"
    probechar_subdir="Langmuir_Kennlinie_SM2400"
    
    for analysisParameters in listOfAnalysisParameters:
        
        datapath      = analysisParameters["datapath"]
        year          = analysisParameters["year"]
        monthday      = analysisParameters["monthday"]
        series        = analysisParameters["series"]
        magFieldLabel = analysisParameters["magFieldLabel"]
        outDir        = analysisParameters["outDir"]
        my_cmap       = analysisParameters["cmap"]
        contourLevels = analysisParameters["contourLevels"]
        minPhiF       = analysisParameters["minPhiF"]
        maxPhiF       = analysisParameters["maxPhiF"]
        minIisat      = analysisParameters["minIisat"]
        maxIisat      = analysisParameters["maxIisat"]
        minR          = analysisParameters["minR"]
        maxR          = analysisParameters["maxR"]
        minZ          = analysisParameters["minZ"]
        maxZ          = analysisParameters["maxZ"]
        
        outputSuffix  = getOutputFilePrefix(analysisParameters)
        
        datadir=datapath+"/"+year+"/"+monthday+"/"+series+"/"
        
        print("looking for data in "+datadir)
        print("       => saving to "+outDir+"/"+outputSuffix+"*")
        
        # read positions
        posfilePiezo=datadir+"/"+positions_subdir+"/path.dat"
        posPiezo = np.loadtxt(posfilePiezo)
        
        xOff = 28.5-30.0 # r
        xScale = 0.9
        
        yOff = 49 # z
        yScale = -0.9
        
        def scale(xLinEnc, off, scale):
            return np.add(np.multiply(xLinEnc, scale), off)
        
        realX = scale(posPiezo[:,0], xOff, xScale)
        realY = scale(posPiezo[:,1], yOff, yScale)
        
        num_r = 23
        num_z = 12
        
        gridX = realX.reshape([num_z,num_r])
        gridY = realY.reshape([num_z,num_r])
        
#        dx = np.mean(np.diff(gridX, axis=1))
#        dy = np.mean(np.diff(gridY, axis=0))
#        
#        # compute corners for pcolormesh
#        extX = np.hstack([gridX-dx/2, np.transpose([gridX[:,-1]+dx/2])])
#        cornersX = np.vstack([extX, [extX[-1,:]]])
#        
#        extY=np.vstack([gridY-dy/2, [gridY[-1,:]+dy/2]])
#        cornersY=np.hstack([extY, np.transpose([extY[:,-1]])])
        
        
        
        #plt.plot(positions[:,0], positions[:,1])
        
        # read vacuum probe characteristics (no plasma)
        vacuum_meas_filename = datapath+"/"+year+"/0716/500mT_no_plasma.dat"
        vacuum_char = np.loadtxt(vacuum_meas_filename)

        # read probe characteristics
        phif = np.zeros([num_z, num_r])
        iisat = np.zeros([num_z, num_r])
        i=0
        chars=[]
        for z in range(num_z):
            for r in range(num_r):    
                
                char_filename=datadir+"/"+probechar_subdir+("/%04d-00.dat"%i)
                #print(char_filename)
                char_data=np.loadtxt(char_filename)
                
                char_data[:,1] -= vacuum_char[:,1]
                
                # add 99.6V to voltage (HP 6516A power supply in series to SourceMeter)
        #        char_data[:,0]+=99.6
                
                chars.append(char_data)
                phif_here=np.interp(0,char_data[:,1],char_data[:,0])
                phif[z,r]=phif_here
                
                #iisat_here=np.interp(phif_here-30.0, char_data[:,0], char_data[:,1])
                iisat_here=np.interp(0.0, char_data[:,0], char_data[:,1])
                iisat[z,r]=iisat_here
                
                i+=1
                

        
        
#        def onclick(event):
#        #    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
#        #          ('double' if event.dblclick else 'single', event.button,
#        #           event.x, event.y, event.xdata, event.ydata))
#            r_idx=np.int(np.round(event.xdata)/2.0)
#            z_idx=np.int((gridY[0,0]-np.round(event.ydata-1.0))/2.0)
#            #print(r_idx)
#            #print(z_idx)
#            print("r=%g z=%g"%(gridX[z_idx,r_idx],gridY[z_idx,r_idx]))
#            char_idx=z_idx*num_r+r_idx
#            
#            fig42 = plt.figure(42)
#            plt.plot(chars[char_idx][:,0], chars[char_idx][:,1]*1e6, label="r=%g z=%g"%(gridX[z_idx,r_idx],gridY[z_idx,r_idx]))
#        #    plt.plot(chars[char_idx][:,0], chars[char_idx][:,2], label="dI/dU")
#            plt.xlabel("U / V")
#            plt.ylabel("I / uA")
#            plt.legend(loc="best")
#            plt.grid(True)
#            plt.tight_layout()
#            plt.show()
        
        
        
        
        #fig=plt.figure()     
        #plt.imshow(np.flipud(phif[z0:,:]), interpolation="nearest", extent=[r_meshgrid[z0,0]-1, r_meshgrid[z0,-1]+1, z_meshgrid[z0,0]+1, z_meshgrid[-1,0]-1])
        #plt.contourf(np.flipud(phif[z0:,:]), interpolation="nearest", extent=[gridX[z0,0]-1, gridX[z0,-1]+1, gridY[z0,0]+1, gridY[-1,0]-1], levels=20, cmap=my_cmap)
        #plt.pcolormesh(cornersX, cornersY, phif, cmap=my_cmap)
        plt.figure()
        plt.contourf(gridX, gridY, phif, cmap=plt.get_cmap(my_cmap), levels=contourLevels, vmin=minPhiF, vmax=maxPhiF)
        cb=plt.colorbar()
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-2,2))
        cb.update_ticks() 
        cb.set_label("floating potential / V")
        plt.xlabel("r / mm")
        plt.ylabel("z / mm")
        plt.xlim([minR, maxR])
        plt.ylim([minZ, maxZ])
        plt.title(magFieldLabel)
        plt.axis("equal")
        plt.tight_layout()
        plt.savefig(outDir+"/PhiF_"+outputSuffix+".png")
#        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        
        #fig2=plt.figure()     
        #plt.imshow(np.flipud(iisat[z0:,:]*-1e6), interpolation="nearest", extent=[r_meshgrid[z0,0], r_meshgrid[z0,-1], z_meshgrid[z0,0], z_meshgrid[-1,0]])
        #plt.contourf(np.flipud(iisat[z0:,:]*-1e6), interpolation="nearest", extent=[gridX[z0,0]-1, gridX[z0,-1]+1, gridY[z0,0]+1, gridY[-1,0]-1], levels=20, cmp=my_cmap)
        plt.figure()
        plt.contourf(gridX, gridY, iisat*-1e6, cmap=plt.get_cmap(my_cmap), levels=contourLevels, vmin=minIisat, vmax=maxIisat)
        cb=plt.colorbar()
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-2,2))
        cb.update_ticks() 
        cb.set_clim(vmin=minIisat, vmax=maxIisat)
        cb.set_label("ion saturation current (at 30V below phif) / uA")
        plt.xlabel("r / mm")
        plt.ylabel("z / mm")
        plt.xlim([minR, maxR])
        plt.ylim([minZ, maxZ])
        plt.title(magFieldLabel)
        plt.axis("equal")
        plt.tight_layout()
        plt.savefig(outDir+"/Iisat_"+outputSuffix+".png")
#        cid = fig2.canvas.mpl_connect('button_press_event', onclick)


runAnalysis(listOfAnalysisParameters)
