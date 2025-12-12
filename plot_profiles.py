# -*- coding: utf-8 -*-
#
# read 2d Langmuir measurement and plot profile of floating potential
#
# Jonathan Schilling (jonathan.schilling@ipp.mpg.de)


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from lmfit import Model

from copy import deepcopy

from pathlib import Path

# folder = Path("/home/jonathan/Uni/Forschung")
folder = Path("/data/jonathan/Uni/Forschung")

listOfAnalysisParameters = []

# common parameters
analysisParameters = {}
analysisParameters["datapath"]= folder / "01_Messdaten" / "Messdaten_2D_Langmuir_Juli2018"
analysisParameters["outDir"] = folder / "03_Auswertung" / "2D_Langmuir_Juli2018" / "out"
analysisParameters["cmap"] = "jet"
analysisParameters["contourLevels_PhiF"] = 40 + 1
analysisParameters["contourLevels_Iisat"] = 20 + 1
analysisParameters["contourLevels_Te"] = 20 + 1
analysisParameters["minPhiF"] = 60 # V
analysisParameters["maxPhiF"] = 80 # V
analysisParameters["minIisat"] =  0 # uA
analysisParameters["maxIisat"] = 20 # uA
analysisParameters["minTe"] = 2.0 # eV
analysisParameters["maxTe"] = 4.0 # eV
analysisParameters["minR"] = -5 # mm
analysisParameters["maxR"] = 40 # mm
analysisParameters["minZ"] =  0 # mm
analysisParameters["maxZ"] = 30 # mm
analysisParameters["year"]="2018"

# build individual runs
analysisParameters["magFieldLabel"] = "0.0 T"
analysisParameters["monthday"]="0716"
analysisParameters["series"]="001"
listOfAnalysisParameters.append(deepcopy(analysisParameters))

# analysisParameters["magFieldLabel"] = "0.05 T"
# analysisParameters["monthday"]="0716"
# analysisParameters["series"]="007_008_merged"
# listOfAnalysisParameters.append(deepcopy(analysisParameters))
#
# analysisParameters["magFieldLabel"] = "0.5 T"
# analysisParameters["monthday"]="0716"
# analysisParameters["series"]="002"
# listOfAnalysisParameters.append(deepcopy(analysisParameters))
#
# analysisParameters["magFieldLabel"] = "1.5 T"
# analysisParameters["monthday"]="0716"
# analysisParameters["series"]="004"
# listOfAnalysisParameters.append(deepcopy(analysisParameters))
#
# analysisParameters["magFieldLabel"] = "4.0 T"
# analysisParameters["monthday"]="0716"
# analysisParameters["series"]="005"
# listOfAnalysisParameters.append(deepcopy(analysisParameters))


########## execution below ##########

def getOutputFilePrefix(analysisParameters):
    year          = analysisParameters["year"].replace(" ", "_")
    monthday      = analysisParameters["monthday"].replace(" ", "_")
    series        = analysisParameters["series"].replace(" ", "_")
    magFieldLabel = analysisParameters["magFieldLabel"].replace(" ", "_")
    my_cmap       = analysisParameters["cmap"].replace(" ", "_")
    contourLevels_PhiF = str(analysisParameters["contourLevels_PhiF"]).replace(" ", "_")
    contourLevels_Iisat = str(analysisParameters["contourLevels_Iisat"]).replace(" ", "_")
    return magFieldLabel+"_"+year+"_"+monthday+"_"+series+"_"+my_cmap+"_"+contourLevels_PhiF+"_"+contourLevels_Iisat

def plotElectrodes(plt):
    stainlessSteelColor="#333333"
    macorColor="#AAAAAA"

    # bottom electrode
    plt.fill([-30.0,  30.0,  30.0, -30.0], \
             [  0.0,   0.0, -10.0, -10.0], stainlessSteelColor)
    # bottom shield
    plt.fill([-33.5, -32.0, -32.0, -33.5], \
             [  0.0,   0.0, -10.0, -10.0], stainlessSteelColor)
    plt.fill([ 33.5,  32.0,  32.0,  33.5], \
             [  0.0,   0.0, -10.0, -10.0], stainlessSteelColor)
    # bottom Macor insulation
    plt.fill([-32.0, -30.0, -30.0, -32.0], \
             [- 5.0, - 5.0, -10.0, -10.0], macorColor)
    plt.fill([ 32.0,  30.0,  30.0,  32.0], \
             [- 5.0, - 5.0, -10.0, -10.0], macorColor)

    # top electrode
    plt.fill([-30.0,  30.0,  30.0, -30.0], \
             [ 30.0,  30.0,  40.0,  40.0], stainlessSteelColor)
    # top shield
    plt.fill([-33.5, -32.0, -32.0, -33.5], \
             [ 30.0,  30.0,  40.0,  40.0], stainlessSteelColor)
    plt.fill([ 33.5,  32.0,  32.0,  33.5], \
             [ 30.0,  30.0,  40.0,  40.0], stainlessSteelColor)
    # top Macor insulation
    plt.fill([-32.0, -30.0, -30.0, -32.0], \
             [ 35.0,  35.0,  40.0,  40.0], macorColor)
    plt.fill([ 32.0,  30.0,  30.0,  32.0], \
             [ 35.0,  35.0,  40.0,  40.0], macorColor)

def TeFun(x, I_ion, I_e, T_e):
    #return I_ion + S_ion*x + I_e*np.exp((x-Phi_P)/T_e)
    return I_ion + I_e*np.exp((x-80.0)/T_e)
Te_model = Model(TeFun)


def runAnalysis(listOfAnalysisParameters):
    positions_subdir="2D_Control_Piezo_Linearencoder"
    probechar_subdir="Langmuir_Kennlinie_SM2400"

    for analysisParameters in listOfAnalysisParameters:

        datapath            = analysisParameters["datapath"]
        year                = analysisParameters["year"]
        monthday            = analysisParameters["monthday"]
        series              = analysisParameters["series"]
        magFieldLabel       = analysisParameters["magFieldLabel"]
        outDir              = analysisParameters["outDir"]
        my_cmap             = analysisParameters["cmap"]
        contourLevels_PhiF  = analysisParameters["contourLevels_PhiF"]
        contourLevels_Iisat = analysisParameters["contourLevels_Iisat"]
        contourLevels_Te = analysisParameters["contourLevels_Te"]
        minPhiF             = analysisParameters["minPhiF"]
        maxPhiF             = analysisParameters["maxPhiF"]
        minIisat            = analysisParameters["minIisat"]
        maxIisat            = analysisParameters["maxIisat"]
        minTe            = analysisParameters["minTe"]
        maxTe            = analysisParameters["maxTe"]
        minR                = analysisParameters["minR"]
        maxR                = analysisParameters["maxR"]
        minZ                = analysisParameters["minZ"]
        maxZ                = analysisParameters["maxZ"]

        outputSuffix  = getOutputFilePrefix(analysisParameters)

        datadir = datapath / year / monthday / series

        print("looking for data in " + str(datadir))
        print("       => saving to "+ str(outDir / outputSuffix))

        # read positions
        posfilePiezo = datadir / positions_subdir / "path.dat"
        posPiezo = np.loadtxt(posfilePiezo)

        xOff = 28.5-30.0 # r
        #xScale = 0.9 # maybe camera viewing angle error?
        xScale = 1.0

        yOff = 54 # z
        #yScale = -0.9 # maybe camera viewing angle error?
        yScale = -1.0

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
        vacuum_meas_filename = datapath / year / "0716" / "500mT_no_plasma.dat"
        vacuum_char = np.loadtxt(vacuum_meas_filename)


#        r_profile_r = []
#        r_profile_Te = []


#        plt.figure()

        # read probe characteristics
        phif = np.zeros([num_z, num_r])
        iisat = np.zeros([num_z, num_r])
        te = np.zeros([num_z, num_r])
        i=0
        chars=[]
        for z in range(num_z):
            for ir,r in enumerate(range(num_r)):

                char_filename = datadir / probechar_subdir / f"{i:04d}-00.dat"
                #print(char_filename)
                char_data=np.loadtxt(char_filename)

                char_data[:,1] -= vacuum_char[:,1]

                # add 99.6V to voltage (HP 6516A power supply in series to SourceMeter)
        #        char_data[:,0]+=99.6

                chars.append(char_data)
                phif_here=np.interp(0, char_data[:,1],char_data[:,0])
#                phif_here=np.interp(1e-6,char_data[:,1],char_data[:,0])
                phif[z,r]=phif_here

#                iisat_here=np.interp(phif_here-50.0, char_data[:,0], char_data[:,1])
                iisat_here=np.interp(-30.0, char_data[:,0], char_data[:,1])
                iisat[z,r]=iisat_here


                #if z==num_z/2:
                    # radial profile, vertically centered between electrodes
                 #   r_profile_r.append(r)

                    # indices up to floating potential + 1V
                span = 12 # 8 works ok for most of the 0T data
                lastIdx = np.where(char_data[:,0] < phif_here)[0][-1]

                Uspan_n = 4.0
                Uspan_p = 2.0
                indices = np.where((char_data[:,0] > phif_here-Uspan_n)&(char_data[:,0] < phif_here+Uspan_p))[0]

                  #  plt.plot(char_data[lastIdx-span:lastIdx+span+2,0], char_data[lastIdx-span:lastIdx+span+2,1], ".", color="C"+str(ir%10), label="r="+str(r))

                    # fit linear model for ion saturation region plus exponential model for electron retardation region
                    # starting values and fit
                try:
                    #result = Te_model.fit(char_data[lastIdx-span:lastIdx+span+2,1], x=char_data[lastIdx-span:lastIdx+span+2,0], I_ion=-1e-5, I_e=1e-5, T_e=2.0)
                    result = Te_model.fit(char_data[indices,1], x=char_data[indices,0], I_ion=-1e-5, I_e=1e-5, T_e=2.0)
                    te[z,r] = result.best_values['T_e']
                except ValueError:
                    te[z,r] = np.nan

                    # print some debugging output of the fit
                   # print(result.fit_report())



                  #  r_profile_Te.append(result.best_values['T_e'])

                    # plot fit result
                   # subset_l = char_data[lastIdx-span:lastIdx+span+2,0][-1]
                   # subset_f = char_data[lastIdx-span:lastIdx+span+2,0][0]
                   # subset_w = subset_l-subset_f

                   # x_eval = np.linspace(subset_f-subset_w*0.05, subset_l+subset_w*0.05, 100)


                   # plt.plot(x_eval, result.eval(x=x_eval), color="C"+str(ir%10))






                i+=1

#
#        plt.xlabel("U / V")
#        plt.ylabel("I / A")
#        plt.legend(loc="upper left")
#        plt.grid(True)
#
#        plt.figure()
#        plt.plot(r_profile_r, r_profile_Te, ".-", label="Te")
#        plt.xlabel("r / m")
#        plt.ylabel("Te / eV")
#        plt.grid(True)
#        plt.legend(loc="best")

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
        plt.contourf(gridX, gridY, phif, cmap=plt.get_cmap(my_cmap), levels=contourLevels_PhiF+1, vmin=minPhiF, vmax=maxPhiF)
        m = cm.ScalarMappable(cmap=my_cmap) # from https://stackoverflow.com/questions/43150687/colorbar-limits-are-not-respecting-set-vmin-vmax-in-plt-contourf-how-can-i-more
        m.set_array(phif)
        m.set_clim(minPhiF, maxPhiF)
        cb=plt.colorbar(m, ax=plt.gca(), boundaries=np.linspace(minPhiF, maxPhiF, contourLevels_PhiF))
#         cb = clippedcolorbar(CS, extend='neither')
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-2,2))
        cb.update_ticks()
        cb.set_label("floating potential / V")
        plt.xlabel("r / mm")
        plt.ylabel("z / mm")
        plotElectrodes(plt)
        plt.axis("equal")
        plt.xlim([minR, maxR])
        plt.ylim([minZ, maxZ])
        plt.title(magFieldLabel)
        plt.tight_layout()
        plt.savefig(outDir / f"PhiF_{outputSuffix}.png")
#        cid = fig.canvas.mpl_connect('button_press_event', onclick)
#
        #fig2=plt.figure()
        #plt.imshow(np.flipud(iisat[z0:,:]*-1e6), interpolation="nearest", extent=[r_meshgrid[z0,0], r_meshgrid[z0,-1], z_meshgrid[z0,0], z_meshgrid[-1,0]])
        #plt.contourf(np.flipud(iisat[z0:,:]*-1e6), interpolation="nearest", extent=[gridX[z0,0]-1, gridX[z0,-1]+1, gridY[z0,0]+1, gridY[-1,0]-1], levels=20, cmp=my_cmap)
        plt.figure()
        #plt.pcolormesh(gridX, gridY, iisat*-1e6, cmap=plt.get_cmap(my_cmap), vmin=minIisat, vmax=maxIisat)
        plt.contourf(gridX, gridY, iisat*-1e6, cmap=plt.get_cmap(my_cmap), levels=contourLevels_Iisat+1, vmin=minIisat, vmax=maxIisat)
        m = cm.ScalarMappable(cmap=my_cmap)
        m.set_array(iisat*-1e6)
        m.set_clim(minIisat, maxIisat)
        cb=plt.colorbar(m, ax=plt.gca(), boundaries=np.linspace(minIisat, maxIisat, contourLevels_Iisat))
#         cb = clippedcolorbar(CS, extend='neither')
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-2,2))
        cb.update_ticks()
        cb.set_label("ion saturation current (at -30V) / uA")
        plt.xlabel("r / mm")
        plt.ylabel("z / mm")
        plotElectrodes(plt)
        plt.axis("equal")
        plt.xlim([minR, maxR])
        plt.ylim([minZ, maxZ])
        plt.title(magFieldLabel)
        plt.tight_layout()
        plt.savefig(outDir / f"Iisat_{outputSuffix}.png")
#        cid = fig2.canvas.mpl_connect('button_press_event', onclick)


        plt.figure()
        plt.contourf(gridX, gridY, te, cmap=plt.get_cmap(my_cmap), levels=contourLevels_Te+1, vmin=minTe, vmax=maxTe)
        m = cm.ScalarMappable(cmap=my_cmap)
        m.set_array(te)
        m.set_clim(minTe, maxTe)
        cb=plt.colorbar(m, ax=plt.gca(), boundaries=np.linspace(minTe, maxTe, contourLevels_Te))
#        cb = clippedcolorbar(CS, extend='neither')
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-2,2))
        cb.update_ticks()
        cb.set_label("T_e / eV")
        plt.xlabel("r / mm")
        plt.ylabel("z / mm")
        plotElectrodes(plt)
        plt.axis("equal")
        plt.xlim([minR, maxR])
        plt.ylim([minZ, maxZ])
        plt.title(magFieldLabel)
        plt.tight_layout()
        plt.savefig(outDir / f"Te_{outputSuffix}.png")

        plt.show()


runAnalysis(listOfAnalysisParameters)
