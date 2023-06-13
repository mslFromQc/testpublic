# -*- coding: utf-8 -*-
"""
Created on Wed May  4 08:17:51 2022

@author: guillaume.p-april
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import  bisplrep, bisplev, griddata, splrep, splev
import GeometryClasses as geo




# ============================================================================
# single evaluation position
u_eval = -10e-3
v_eval = -5e-3
skewAngle =  1/4*np.pi
wedgeAngle = 30*np.pi/180# 1/4*np.pi
squintAngle = 0#15*np.pi/180# 0.16*np.pi
roofAngle = 0#10*np.pi/180# 0.1*np.pi
# specimen_name = 'cylinder_helix'
specimen_name = 'cylinder_loop'


# scan pattern in (u,v)
along_u = np.linspace(-0e-3, 15e-3, 10)
along_v = np.linspace(+3e-3, +13e-3, 10)
scan_name  = 'GRID_rasterUV_stepsAlongU'
# scan_name  = 'GRID_rasterUV_stepsAlongV'
# scan_name  = 'GRID_combUV_stepsAlongU'
# scan_name  = 'GRID_combUV_stepsAlongV'
#scan_name  = 'LINE_UV'
# ============================================================================


# ============================================================================
# PROPERTY STRUCTURE

# the wedge and probe
Nel = 10
prop = {'Wedge' : {
            'length' : 100e-3,
            'width' : 50e-3,
            'radiusAlongYw' : -400e-3, # the wedge is assumed to match exactly the surface
            'radiusAlongXw' : np.inf,
            'skewAngle' : skewAngle,
            'wedgeAngle' : wedgeAngle,
            'squintAngle' : squintAngle,
            'roofAngle' : roofAngle,
            'NormalPrimaryOffset' : 30e-3,
            'NormalSecondaryOffset' : 10e-3,
            'NormalTertiaryOffset' : 15e-3},
        'Probe' : {
            'Nel' : Nel, 
            'XpYpPlane_SecondaryAxisSide' : 10e-3,
            'XpYpPlane_PrimaryAxisSide' : 30e-3,
            'Ref_SecondaryAxis' : 0,
            'Ref_PrimaryAxis' : -5e-3,
            'ElPos_SecondaryAxis' : np.zeros( (Nel, ) ),
            'ElPos_PrimaryAxis' :   np.linspace(0, 9, Nel)*2.0e-3,
            'ElPos_NormalAxis'  :   np.zeros( (Nel, ) ),
            'ElNormal_SecondaryAxis' : np.zeros( (Nel, ) ),
            'ElNormal_PrimaryAxis' :   np.zeros( (Nel, ) ),
            'ElNormal_NormalAxis' :    np.ones( (Nel, ) ),
            'ElSideSecondary' : np.ones( (Nel, ) )*10e-3,
            'ElSidePrimary' : np.ones( (Nel, ) )*1.5e-3
            }
        }
    

# the specimen
if specimen_name == 'cylinder_loop':
    # cylinder parameters, closed loop reference line
    
    # physical 
    diam = 50e-3 # OD [m]
    R0 = diam/2 # radius [m]
    closedLineDiameter = 30e-3 # closed line diameter [m]
    refLineOrigin = [0, R0, closedLineDiameter/2] # origin on reference line, in general coordinates
    
    # display
    L = 60e-3 # length of cylinder specimen [m]
    NSamplesPhi = 100
    NSamplesAxis = 200
    NSamplesLength = 1001
    characteristic_length_for_display = 0.05*L
    
    # grid
    deltaU_and_deltaV = 1e-3 # grid step along u and v (equal)
    u_start = 0e-3 
    u_end = 55e-3
    v_start = 0e-3
    v_end = 10e-3
    
    prop['Specimen'] = {
        'name' :  'cylinder with closed loop line',
        'outerDiameter' : diam ,
        'closedLineDiameter' : closedLineDiameter ,
        'refLineOrigin' : refLineOrigin}
    prop['Display'] = { 
    'axisLength' : L, 
        'NSamplesPhi' : NSamplesPhi,
        'NSamplesAxis' : NSamplesAxis ,
        'NSamplesLine' : NSamplesLength,
        'characteristic_length_for_display' : characteristic_length_for_display
        }
    prop['Grid'] = {
        'deltaU_and_deltaV': deltaU_and_deltaV,
        'u_start': u_start, 
        'u_end': u_end, 
        'v_start': v_start, 
        'v_end': v_end
        }
        
    
    
    MySpecimen = geo.CylinderCLosedLine(prop)
    MyWedge = geo.SystemAxes(prop)
    #MyWedge.get_wedgePositionAndAxes( DebugFigs = True )

elif specimen_name == 'cylinder_helix':
    # cylinder parameters, helix reference line
    
    # physical 
    diam = 50e-3 # OD [m]
    R0 = diam/2 # radius [m]
    Theta =  0.45*np.pi # suface angle w.r.t. cylinder axis [rad]
    refLineOrigin = [0, R0, 0] # origin on reference line, in general coordinates
    
    # display
    L = 60e-3 # length of cylinder specimen [m]
    NSamplesPhi = 100
    NSamplesAxis = 200
    NSamplesLength = 2001
    characteristic_length_for_display = 0.05*L
    
    # grid
    deltaU_and_deltaV = 1e-3 # grid step along u and v (equal)
    u_start = 0e-3 
    u_end = 55e-3
    v_start = 0e-3
    v_end = 15e-3
    
    prop['Specimen'] = {
        'name' :  'cylinder with helix line',
        'outerDiameter' : diam ,
        'axisVsLineAngle' : Theta ,
        'refLineOrigin' : refLineOrigin}
    prop['Display'] = { 
        'axisLength' : L, 
        'NSamplesPhi' : NSamplesPhi,
        'NSamplesAxis' : NSamplesAxis ,
        'NSamplesLine' : NSamplesLength,
        'characteristic_length_for_display' : characteristic_length_for_display}
    prop['Grid'] = {
        'deltaU_and_deltaV': deltaU_and_deltaV,
        'u_start': u_start, 
        'u_end': u_end, 
        'v_start': v_start, 
        'v_end': v_end
        }
    
    MySpecimen = geo.CylinderHelixLine(prop)
    MyWedge = geo.SystemAxes(prop)
    #MyWedge.get_wedgePositionAndAxes( DebugFigs = True )
# ============================================================================



# ============================================================================
# PLOT THE SPECIMEN AND GRID COORDINATES

MySpecimen.get_analyticUVtoXYZgrid( prop, DebugFigs = True )

MySpecimen.get_surface()

MyUVandXYZ = geo.UVandXYZ(Specimen = MySpecimen)

MyUVandXYZ.drawSpecimenGrid( MySpecimen, prop )

# ============================================================================



# ============================================================================
# PLOT THE SPECIMEN AND LOCAL REFERENTIALS IN (U,V) COORDINATES

# initiate uv coordinates on specimen
MyUVandXYZ = geo.UVandXYZ(Specimen = MySpecimen)

# the position on the weld line
MyUVandXYZ.get_PointAndVectorsXYZFromUV_single(Specimen = MySpecimen, u = u_eval, v = 0)

# the position away from the weld line
MyUVandXYZ.get_PointAndVectorsXYZFromUV_single(Specimen = MySpecimen, u = u_eval, v = v_eval)

# skew angle rotated frame
MyUVandXYZ.rotateOnSkewAngle( prop['Wedge']['skewAngle'] )
 
# get the specimen surface
MySpecimen.get_surface()

# plot the specimen and the local referentials for the 2 positions defined above
#MyUVandXYZ.drawSpecimenLast2Positions( MySpecimen, prop )
# ============================================================================



# ============================================================================
# PLOT THE PLANE OF ELEMENTS 
#MyWedge.plot_probeElements( DebugFigs = True )

# ============================================================================




# ============================================================================
# PLOT THE SCAN PATTERN on the specimen

# initiate UV coordinates
MyUVandXYZ = geo.UVandXYZ(Specimen = MySpecimen)

# initiate the scan pattern
MyScanPattern = geo.scanPattern(along_u, along_v)

# evaluate the scan pattern based on scan modality
getattr(MyScanPattern, scan_name)()

# plot the scan pattern on specimen
MyScanPattern.drawSpecimenScan(MySpecimen, MyUVandXYZ, prop)
# ============================================================================

