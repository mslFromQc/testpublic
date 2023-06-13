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
skewAngle =  1/4*np.pi
wedgeAngle = 30*np.pi/180# 1/4*np.pi
squintAngle = 15*np.pi/180# 0.16*np.pi
roofAngle = 10*np.pi/180# 0.1*np.pi
# ============================================================================


# ============================================================================
# PROPERTY STRUCTURE

# the wedge and probe
Nel = 10
prop = {'Wedge' : {
            'length' : 40e-3,
            'width' : 50e-3,
            'radiusAlongYw' : +80e-3, # the wedge is assumed to match exactly the surface
            'radiusAlongXw' : np.inf,
            'skewAngle' : skewAngle,
            'wedgeAngle' : wedgeAngle,
            'squintAngle' : squintAngle,
            'roofAngle' : roofAngle,
            'NormalPrimaryOffset' : -30e-3,
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
    

# ============================================================================
# PLOT THE PLANE OF ELEMENTS 
MyWedge = geo.SystemAxes(prop)
    
MyWedge.plot_probeElements( figs = True )

MyWedge.get_wedgePositionAndAxes( figs = True)
# ============================================================================
