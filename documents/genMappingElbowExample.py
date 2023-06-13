# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 10:41:30 2023

@author: benoit.lepage
"""

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

####display parameters
displayScanPatch = True
scanLinesToDisplay = [9,]  #'all' to display all scan lines, list ex [0,2,3,] to display specific scan lines

mainAxisLineWidth = 5
borderLineWidth = 3
centerLineWidth = 1

#fixed parameter
arc = np.pi


def circleArc(middleAngleRad, arcRad, offCenterRadius, circleRadius, phi,nbSeg):
    theta = np.linspace(middleAngleRad-arcRad/2, middleAngleRad+arcRad/2, nbSeg)
    yinplane = circleRadius*np.cos(theta)+offCenterRadius
    x = yinplane*np.sin(phi)
    y = yinplane*np.cos(phi)
    z = circleRadius*np.sin(theta)
    return x,y,z

def add_patch(legend, patchList):
    from matplotlib.patches import Patch
    ax = legend.axes

    handles, labels = ax.get_legend_handles_labels()
    for i_patch in range(len(patchList)):
        handles.append(Patch(facecolor='w',edgecolor=patchList[i_patch]['patchColor']))
        labels.append(patchList[i_patch]['patchLabel'])

    legend._legend_box = None
    legend._init_legend_box(handles, labels)
    legend._set_loc(legend._loc)
    legend.set_title(legend.get_title().get_text())

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
if displayScanPatch:
    figde, axde = plt.subplots(nrows=1, ncols=1, figsize=(12, 12)) #This is a representation of the dataEncoding grid in gen Mapping

####this test parameter
scanResolution= 0.01  #v = 0 is the ref
indexResolution = 0.01
circleRadius = 0.2
offCenterRadius = 0.5
phiangle = np.pi/4
scanWidth = 0.060
scanOverlap = 0.005



#Nb U is evaluated on the extrados of the pipe to match the target scan resolution
ptsPerIndex = 1
angularStep = indexResolution/(circleRadius*ptsPerIndex)  
pathLenghtExtrados = (circleRadius+offCenterRadius)*phiangle
nbU = int(pathLenghtExtrados/scanResolution)     #Extend of the u axis
nbVscanlines =int(arc*circleRadius/(scanWidth-scanOverlap)) #scan lines orthogonal to the v axis
nbVgrid = int(arc*circleRadius/indexResolution) #Extend of the v axis


#Drawing a true u,v grid
#calculating and plotting the u,v grid
nbSeg = int(arc/angularStep)
xs = np.zeros([nbU,nbSeg])
ys = np.zeros([nbU,nbSeg])
zs = np.zeros([nbU,nbSeg])
refAngle = np.pi/2

for i_us in range(nbU):
    phi = i_us*phiangle/nbU
    xs[i_us,:],ys[i_us,:],zs[i_us,:] = circleArc(refAngle,arc,offCenterRadius,circleRadius,phi,nbSeg)
    if i_us == 0:
        ax.plot(xs[i_us,:],ys[i_us,:],zs[i_us,:],'r', linewidth = mainAxisLineWidth, label ='v axis ref line (u=0)')
    elif i_us == 1:
        ax.plot(xs[i_us,:],ys[i_us,:],zs[i_us,:],'k:', linewidth = centerLineWidth,  label ='regular u,v grid')        
    else:
        ax.plot(xs[i_us,:],ys[i_us,:],zs[i_us,:],'k:', linewidth = centerLineWidth)
    
for i_vsg in range(nbVgrid):
    segIdx = i_vsg*ptsPerIndex
    ax.plot(xs[:,segIdx],ys[:,segIdx],zs[:,segIdx],'k:')      

for i_vs in range(nbVscanlines):
    segIdx = int(i_vs*nbSeg/nbVscanlines)
    if i_vs == 0: 
        ax.plot(xs[:,segIdx],ys[:,segIdx],zs[:,segIdx],'r', linewidth = mainAxisLineWidth, label = 'u axis ref line (v=0) and first scan line')
    elif i_vs == 1:
        ax.plot(xs[:,segIdx],ys[:,segIdx],zs[:,segIdx],'k', linewidth = borderLineWidth, label = 'other scan lines')        
    else:
        ax.plot(xs[:,segIdx],ys[:,segIdx],zs[:,segIdx],'k', linewidth = borderLineWidth)      

if displayScanPatch:
    #calculating and plotting a subset of u',v' grid
    scanWidthRad = scanWidth/circleRadius
    nbSegp = int(scanWidthRad/angularStep)
    vAsVprime = []
    patchList = []
    for i_display in range(nbVscanlines):
        if scanLinesToDisplay == 'all' or i_display in scanLinesToDisplay:
            vAsVprime.append(i_display)
    
    for i_vp in range(len(vAsVprime)):
        #at us[0] == 0 phi = 0 so the angle for a specific v can be found
        thisVpColor = (float((vAsVprime[i_vp]%2)*0.6),float(vAsVprime[i_vp]/nbVscanlines),float(1 - (vAsVprime[i_vp]/nbVscanlines))) #this is only a complicated way to generate automatically different color for each inspection patch
        thisSegIdx = int(vAsVprime[i_vp]*nbSeg/nbVscanlines)
        vpRefAngle = np.arctan2(zs[0,thisSegIdx],ys[0,thisSegIdx]-offCenterRadius)
        
        thisPathLenght = ys[0,thisSegIdx]*phiangle
        nbUp = int(thisPathLenght/scanResolution)  #extend of this uprime map
        nbVp = int(nbSegp/ptsPerIndex)  #extend of this vprime map
        xps = np.zeros([nbUp,nbSegp])
        yps = np.zeros([nbUp,nbSegp])
        zps = np.zeros([nbUp,nbSegp])
        thisPath = {'patchColor':thisVpColor, 'patchLabel':" u'/v' local path# " + str(vAsVprime[i_vp])}
        patchList.append(thisPath)
        
        #drawing the scan patch on the elbow
        for i_ups in range(nbUp):
            phi = i_ups*phiangle/nbUp
            xps[i_ups,:],yps[i_ups,:],zps[i_ups,:] = circleArc(vpRefAngle,scanWidthRad,offCenterRadius,circleRadius,phi,nbSegp)
            if i_ups == 0 or i_ups == nbUp -1:
                linewidth = borderLineWidth
            else:
                linewidth = centerLineWidth
            ax.plot(xps[i_ups,:],yps[i_ups,:],zps[i_ups,:],color = thisVpColor, linewidth = linewidth)
        
        
        for i_vps in range(nbVp):
            segIdx = i_vps*ptsPerIndex
            if i_vps == 0 or i_vps == nbVp-1:
                linewidth = borderLineWidth
            else:
                linewidth = centerLineWidth            
            ax.plot(xps[:,segIdx],yps[:,segIdx],zps[:,segIdx],color = thisVpColor,linewidth = linewidth)   
        
        #drawing the scan path on the dataEncoding map
        usde = np.linspace(0,thisPathLenght,nbUp)
        vsde = np.linspace(i_vp*scanWidth-scanWidth/2,i_vp*scanWidth+scanWidth/2,nbVp)
        for i_vsde in range(nbUp):
            if i_vsde == 0 or i_ups == nbUp -1:
                linewidth = borderLineWidth
            else:
                linewidth = centerLineWidth                    
            axde.plot(np.full(nbVp,usde[i_vsde]),vsde,color = thisVpColor, linewidth = linewidth)
        
        
        for i_usde in range(nbVp):
            segIdx = i_vps*ptsPerIndex
            if i_usde == 0 or i_vsde == nbVp-1:
                linewidth = borderLineWidth
            else:
                linewidth = centerLineWidth            
            axde.plot(usde,np.full(nbUp,vsde[i_usde]),color = thisVpColor, linewidth = linewidth)
        
if displayScanPatch:
    add_patch(ax.legend(),patchList)
    figde.suptitle('Representation in the discrete Grid data Encoding')
    axde.set_xlabel('u axis (in the dataEncoding representation) (~m)')
    axde.set_ylabel('v axis (in the dataEncoding representation) (~m)')
else:
    ax.legend()


ax.set_box_aspect((np.ptp(xs), np.ptp(ys), np.ptp(zs)))
plt.show()

