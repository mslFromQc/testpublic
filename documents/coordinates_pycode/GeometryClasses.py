# -*- coding: utf-8 -*-
"""
Created on Wed May 11 15:52:27 2022

@author: guillaume.p-april
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import  splrep, splev







# GEOMETRY PARAMETERS



# ============================================================================
# cylinder parameters, helix reference line

class CylinderHelixLine :
    
    def __init__(self, prop):
        
        self.prop = prop
        
        R0 = self.prop['Specimen']['outerDiameter']/2
        Theta = self.prop['Specimen']['axisVsLineAngle']
        refLineOrigin = self.prop['Specimen']['refLineOrigin']
        
        L =  self.prop['Display']['axisLength']
        phi = np.linspace(0, 2*np.pi, self.prop['Display']['NSamplesPhi']) #
        
        
        # line parameter
        Nsample = self.prop['Display']['NSamplesLine'] # odd to have P=0 at middle of symmetric display
        if Theta == np.pi/2:
            P = np.unique( np.sort( np.concatenate( (np.linspace(-np.pi*R0, +np.pi*R0, Nsample-1), np.asarray( [0] ) ) ) ) )
            ind = int( np.argmin( np.abs(P) ) )
        else:
            P = np.unique( np.sort( np.concatenate( (np.linspace((-L/2)/np.cos(Theta), (L/2)/np.cos(Theta), Nsample-1), np.asarray( [0] ) ) ) ) )
            ind = int( np.argmin( np.abs(P) ) )
            
        # position of line    
        phi0 = np.arctan2(refLineOrigin[1], refLineOrigin[0]) # origin
        phi = P*np.sin(Theta)/R0
        X_line = R0*np.cos(phi+phi0)
        Y_line = R0*np.sin(phi+phi0)
        Z_line = P*np.cos(Theta)
        
        Normal_X_line = -np.cos(phi+phi0)
        Normal_Y_line = -np.sin(phi+phi0)
        Normal_Z_line = np.zeros( (Nsample, ) )
 
        # origin of reference line
        self.OriginLine =  {
                    'X' : refLineOrigin[0],
                    'Y' : refLineOrigin[1],
                    'Z' : refLineOrigin[2]
                            }
        # reference line, along u        
        self.ParallelToLine = {
                    'Origin_ind' : ind,
                    'X' : X_line,
                    'Y' : Y_line,
                    'Z' : Z_line,
                    'Normal_X' : Normal_X_line,
                    'Normal_Y' : Normal_Y_line,
                    'Normal_Z' : Normal_Z_line
                    }




    def get_surface( self ) :
        L =  self.prop['Display']['axisLength']
        R0 = self.prop['Specimen']['outerDiameter']/2
        z = np.linspace(-L/2, L/2, self.prop['Display']['NSamplesAxis'])
        phi = np.linspace(0, 2*np.pi, self.prop['Display']['NSamplesPhi']) #
        
        phi_grid, Z_surf = np.meshgrid(phi, z)
        X_surf = R0*np.cos(phi_grid)
        Y_surf = R0*np.sin(phi_grid)
        
        # normal to surface
        X_surfNormal = -np.cos(phi_grid)
        Y_surfNormal = -np.sin(phi_grid)
        Z_surfNormal = np.zeros( np.shape(phi_grid) )
        
        self.Surface = {
                    'X' :  X_surf, 
                    'Y' :  Y_surf, 
                    'Z' :  Z_surf,
                    'Normal_X' :  X_surfNormal,
                    'Normal_Y' :  Y_surfNormal,
                    'Normal_Z' :  Z_surfNormal
                            }
    
    
    
    
    def get_LinePerpendicular(self, pos ):
        
        # line perpendicular to reference line, along v
        
        Theta = self.prop['Specimen']['axisVsLineAngle']
        R0 = self.prop['Specimen']['outerDiameter']/2
        L = self.prop['Display']['axisLength']
        
        Nsample = self.prop['Display']['NSamplesLine'] # odd to have P=0 at middle of symmetric display
        Theta_perp = Theta  - np.pi/2
        phi0 = np.arctan2(pos[1], pos[0])
        if np.abs(Theta_perp) == np.pi/2:
            P = np.unique( np.sort( np.concatenate( (np.linspace(-np.pi*R0, +np.pi*R0, Nsample-1), np.asarray( [0] ) ) ) ) )
            ind = int( np.argmin( np.abs(P) ) )
        else:
            P = np.unique( np.sort( np.concatenate( (np.linspace((-L/2-pos[2])/np.cos(Theta_perp), (L/2-pos[2])/np.cos(Theta_perp), Nsample-1), np.asarray( [0] ) ) ) ) )
            ind = int( np.argmin( np.abs(P) ) )
            
        phi = P*np.sin(Theta_perp)/R0
        X_line_perp = R0*np.cos(phi +phi0)
        Y_line_perp = R0*np.sin(phi +phi0)
        Z_line_perp = pos[2] + P*np.cos(Theta_perp)
        
        Normal_X_line_perp = -np.cos(phi +phi0)
        Normal_Y_line_perp = -np.sin(phi +phi0)
        Normal_Z_line_perp = np.zeros( (len(phi), ) )
    
        # line perpendicular to reference line, along v
        self.PerpendicularToLine = {
                    'Origin_ind' : ind, 
                    'X' : X_line_perp,
                    'Y' : Y_line_perp,
                    'Z' : Z_line_perp,
                    'Normal_X' : Normal_X_line_perp,
                    'Normal_Y' : Normal_Y_line_perp,
                    'Normal_Z' : Normal_Z_line_perp
                    }


    def get_analyticUVtoXYZgrid( self, prop , DebugFigs = False):
        
        Theta = self.prop['Specimen']['axisVsLineAngle']
        R0 = self.prop['Specimen']['outerDiameter']/2
        L = self.prop['Display']['axisLength']

        refLineOrigin = self.prop['Specimen']['refLineOrigin']
        X0 = refLineOrigin[0]
        Y0 = refLineOrigin[1]
        Z0 = refLineOrigin[2]

        deltaU_and_deltaV = self.prop['Grid']['deltaU_and_deltaV']
        u_start = self.prop['Grid']['u_start']
        u_end = self.prop['Grid']['u_end']
        v_start = self.prop['Grid']['v_start']
        v_end = self.prop['Grid']['v_end']

        uu = np.arange(u_start, u_end, deltaU_and_deltaV)
        vv = np.arange(v_start, v_end, deltaU_and_deltaV)
        
        U, V = np.meshgrid(uu, vv)
        
        phi = np.cos(Theta)/R0*V + np.sin(Theta)/R0*U
        phi0 = np.arctan2(Y0, X0)
        
        
        X = R0*np.cos( phi + phi0 )
        Y = R0*np.sin( phi + phi0 )
        Z = -np.sin(Theta)*V + np.cos(Theta)*U + Z0

        u_x = -np.sin(phi + phi0)*np.sin(Theta)
        u_y =  np.cos(phi + phi0)*np.sin(Theta)
        u_z =  np.cos(Theta)*np.ones( np.shape(u_x) )

        v_x = -np.sin(phi + phi0)*np.cos(Theta)
        v_y =  np.cos(phi + phi0)*np.cos(Theta)
        v_z = -np.sin(Theta)*np.ones( np.shape(u_x) )

        # all grid points correspond
        self.Grid = {
            'u' : U, # the grid u coordinate values
            'v' : V, # the grid v coordinate values
            'X' : X, # the grid coordinates in general coordinates X
            'Y' : Y, # the grid coordinates in general coordinates Y
            'Z' : Z, # the grid coordinates in general coordinates Z
            'u_X' : u_x, # the grid local unit vector along u in gen coordinates, X
            'u_Y' : u_y, # the grid local unit vector along u in gen coordinates, Y
            'u_Z' : u_z,  # the grid local unit vector along u in gen coordinates, Z
            'v_X' : v_x, # the grid local unit vector along v in gen coordinates, X
            'v_Y' : v_y, # the grid local unit vector along v in gen coordinates, Y
            'v_Z' : v_z  # the grid local unit vector along v in gen coordinates, Z
            }

# ============================================================================
# cylinder parameters, closed loop reference line

class CylinderCLosedLine :
    
    def __init__(self, prop):
        
        self.prop = prop
        
        R0 = self.prop['Specimen']['outerDiameter']/2
        rho = self.prop['Specimen']['closedLineDiameter']/2
        refLineOrigin = self.prop['Specimen']['refLineOrigin']
        
        L =  self.prop['Display']['axisLength']
        
        
        # line parameter
        Nsample = self.prop['Display']['NSamplesLine'] # odd to have P=0 at middle of symmetric display
        P = np.unique( np.sort( np.concatenate( (np.linspace(-1*np.pi, +1*np.pi, Nsample-1), np.asarray( [0] ) ) ) ) )
        ind = int( np.argmin( np.abs(P) ) )
            
        # position of line    
        phi = P
        phi0 = np.arctan2(refLineOrigin[0], refLineOrigin[2])
        X_line = rho*np.sin(phi + phi0)
        Y_line = np.sqrt(R0**2 - X_line**2)
        Z_line = rho*np.cos(phi + phi0)
        
        Normal_X_line = -X_line/np.sqrt(X_line**2 + Y_line**2)
        Normal_Y_line = -Y_line/np.sqrt(X_line**2 + Y_line**2)
        Normal_Z_line = np.zeros( (Nsample, ) )
 
        # origin of reference line
        self.OriginLine =  {
                    'X' : refLineOrigin[0],
                    'Y' : refLineOrigin[1],
                    'Z' : refLineOrigin[2]
                            }
        
        # reference line, along  u        
        self.ParallelToLine = {
                    'Origin_ind' : ind,
                    'X' : X_line,
                    'Y' : Y_line,
                    'Z' : Z_line,
                    'Normal_X' : Normal_X_line,
                    'Normal_Y' : Normal_Y_line,
                    'Normal_Z' : Normal_Z_line
                    }




    def get_surface( self ) :
        L =  self.prop['Display']['axisLength']
        z = np.linspace(-L/2, L/2, self.prop['Display']['NSamplesAxis'])
        phi = np.linspace(0, 2*np.pi, self.prop['Display']['NSamplesPhi']) #
        R0 = self.prop['Specimen']['outerDiameter']/2
        
        
        phi_grid, Z_surf = np.meshgrid(phi, z)
        X_surf = R0*np.cos(phi_grid)
        Y_surf = R0*np.sin(phi_grid)
        
        # normal to surface
        X_surfNormal = -np.cos(phi_grid)
        Y_surfNormal = -np.sin(phi_grid)
        Z_surfNormal = np.zeros( np.shape(phi_grid) )
        
        self.Surface = {
                    'X' :  X_surf, 
                    'Y' :  Y_surf, 
                    'Z' :  Z_surf,
                    'Normal_X' :  X_surfNormal,
                    'Normal_Y' :  Y_surfNormal,
                    'Normal_Z' :  Z_surfNormal
                            }
    
    
    
    
    def get_LinePerpendicular(self, pos, DebugFigs = False ):
        
        R0 = self.prop['Specimen']['outerDiameter']/2
        rho = self.prop['Specimen']['closedLineDiameter']/2
        
        # angular position of the point 
        phi0 = np.arctan2(pos[0], pos[2])
        
        # local tangent vector
        u_x = rho*np.cos(phi0)
        u_y = -rho**2/np.sqrt(R0**2 - rho**2*np.sin(phi0)**2)*np.cos(phi0)*np.sin(phi0)
        u_z = -rho*np.sin(phi0)

        Norm = np.sqrt( u_x**2 + u_y**2 + u_z**2 )
        u_x = u_x/Norm
        u_y = u_y/Norm
        u_z = u_z/Norm

        # local normal vector        
        phi0 = np.arctan2(pos[1], pos[0])
        Normal_x = -np.cos(phi0)
        Normal_y = -np.sin(phi0)
        Normal_z = 0
        
        v_x =  Normal_y*u_z - Normal_z*u_y
        v_y = -Normal_x*u_z + Normal_z*u_x
        v_z =  Normal_x*u_y - Normal_y*u_x


        # find the local angle with respect tp pipe weld
        if u_z < 0 and v_z > 0:
            # Theta_perp =  -np.arccos( v_z ) 
            Theta_perp =  -np.arccos( u_z ) - np.pi/2
        elif u_z < 0 and v_z < 0 :
            # Theta_perp =  -np.arccos( v_z )
            Theta_perp = +np.arccos( u_z ) - np.pi/2 #- np.pi/2 
        elif u_z > 0 and v_z < 0 :
            # Theta_perp = +np.arccos( v_z )
            Theta_perp = +np.arccos( u_z ) - np.pi/2
        elif u_z > 0 and v_z > 0 :
            # Theta_perp = +np.arccos( v_z )
            Theta_perp = -np.arccos( u_z ) - np.pi/2

        #Theta_perp  = +np.arccos( u_z ) - np.pi/2
        
        Nsample = self.prop['Display']['NSamplesLine'] # odd to have P=0 at middle of symmetric display

        # this is where the perpendicular line parameter is sampled
        P = np.unique( np.sort( np.concatenate( (np.linspace(-np.pi*R0, +np.pi*R0, Nsample-1), np.asarray( [0] ) ) ) ) )
        ind = int( np.argmin( np.abs(P) ) )
        
        # perpendicular line coordinates evaluation
        phi = P*np.sin(Theta_perp)/R0
        X_line_perp = R0*np.cos(phi +phi0)
        Y_line_perp = R0*np.sin(phi +phi0)
        Z_line_perp = pos[2] + P*np.cos(Theta_perp)
        
        # normal to perp line
        Normal_X_line_perp = -np.cos(phi +phi0)
        Normal_Y_line_perp = -np.sin(phi +phi0)
        Normal_Z_line_perp = np.zeros( (len(phi), ) )
        

        if DebugFigs == True:        
            fig = plt.figure(figsize = (15, 14))
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(self.ParallelToLine['X'], self.ParallelToLine['Y'], self.ParallelToLine['Z'], c = 'orange', linewidth = 2, label = 'parallel line')
            ax.plot(X_line_perp, Y_line_perp, Z_line_perp, c ='lightgreen', linewidth = 2, label = 'perpendicular line')
            max_lim = np.max([np.abs(X_line_perp), np.abs(Y_line_perp), np.abs(Z_line_perp)])
            ax.plot([pos[0], pos[0]+0.05*max_lim*u_x], 
                    [pos[1], pos[1]+0.05*max_lim*u_y], 
                    [pos[2], pos[2]+0.05*max_lim*u_z], c = 'red')
            ax.plot([pos[0], pos[0]+0.05*max_lim*v_x], 
                    [pos[1], pos[1]+0.05*max_lim*v_y], 
                    [pos[2], pos[2]+0.05*max_lim*v_z], c = 'green')
            
            ax.text( pos[0]+0.05*max_lim*u_x, pos[1]+0.05*max_lim*u_y, pos[2]+0.05*max_lim*u_z, 'u', c = 'red' )
            ax.text( pos[0]+0.05*max_lim*v_x, pos[1]+0.05*max_lim*v_y, pos[2]+0.05*max_lim*v_z, 'v', c = 'green' )

            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

            ax.legend()            

            ax.set_xlim(-max_lim, max_lim)
            ax.set_ylim(-max_lim, max_lim)
            ax.set_zlim(-max_lim, max_lim) 
            ax.set_box_aspect([1,1,1])
        
        # line perpendicular to reference line, along v
        self.PerpendicularToLine = {
                    'Origin_ind' : ind, 
                    'X' : X_line_perp,
                    'Y' : Y_line_perp,
                    'Z' : Z_line_perp,
                    'Normal_X' : Normal_X_line_perp,
                    'Normal_Y' : Normal_Y_line_perp,
                    'Normal_Z' : Normal_Z_line_perp
                    }


    def get_analyticUVtoXYZgrid( self, prop, DebugFigs = False):
        
        
        R0 = prop['Specimen']['outerDiameter']/2
        rho = prop['Specimen']['closedLineDiameter']/2
        L = prop['Display']['axisLength']
        
        
        X = self.ParallelToLine['X'] 
        Y = self.ParallelToLine['Y'] 
        Z = self.ParallelToLine['Z'] 
        refLineOrigin = self.prop['Specimen']['refLineOrigin']
        X0 = refLineOrigin[0]
        Y0 = refLineOrigin[1]
        Z0 = refLineOrigin[2]

         
        deltaU_and_deltaV = self.prop['Grid']['deltaU_and_deltaV']
    
        # local line u
        u = np.zeros( (len(X), ) )
        u0 = 0
        origin_found = False
        for p in np.arange(1, len(u)):
            du = np.sqrt(  (X[p]-X[p-1])**2 
                         + (Y[p]-Y[p-1])**2 
                         + (Z[p]-Z[p-1])**2 )
            du_test = np.sqrt(  (X0-X[p-1])**2 
                         + (Y0-Y[p-1])**2 
                         + (Z0-Z[p-1])**2 )
            if du_test < du and origin_found == False:
                u0 = u[p-1] + du_test
                
                origin_found = True
                
            u[p] = u[p-1] + du
                
        u = u - u0
        
        tck_u_x = splrep( u, X, k=3, s=0 )
        tck_u_y = splrep( u, Y, k=3, s=0 )
        tck_u_z = splrep( u, Z, k=3, s=0 )
         
                    
        u_start = self.prop['Grid']['u_start']
        u_end = self.prop['Grid']['u_end']
        v_start = self.prop['Grid']['v_start']
        v_end = self.prop['Grid']['v_end']

        uu = np.arange(u_start, u_end, deltaU_and_deltaV)
        vv = np.arange(v_start, v_end, deltaU_and_deltaV)
        
        X = splev(uu, tck_u_x, der = 0)
        Y = splev(uu, tck_u_y, der = 0)
        Z = splev(uu, tck_u_z, der = 0)
        
        theta = np.arctan2(X, Z)
        theta0 = np.arctan2(X0, Z0)
         
        if DebugFigs == True:
            fig = plt.figure(figsize = (13, 12))
            ax = fig.add_subplot(111)
            ax.scatter(uu, theta*180/np.pi)
            
            ax.tick_params(axis='both', which='major', labelsize=15)
            ax.set_xlabel('$u$ [m]', fontsize = 20)
            ax.set_ylabel(r'$\theta$ [deg]', fontsize = 20)
        
        U, V = np.meshgrid(uu, vv)
        THETA, V = np.meshgrid(theta, vv)
        XX, _ = np.meshgrid(X, vv)
        YY, _ = np.meshgrid(Y, vv)
        ZZ, _ = np.meshgrid(Z, vv)

        sinTheta = -R0*np.cos(THETA+theta0)/np.sqrt( R0**2-rho**2*np.sin(THETA+theta0)**4 )
        cosTheta = -np.sqrt( R0**2-rho**2*np.sin(THETA+theta0)**2 )/np.sqrt( R0**2-rho**2*np.sin(THETA+theta0)**4 )*np.sin(THETA+theta0)

        
        phi = cosTheta/R0*V 
        phi0 = np.arctan2(YY, XX)
        X = R0*np.cos( phi + phi0 )
        Y = R0*np.sin( phi + phi0 )
        Z = -sinTheta*V + ZZ

        u_x = -np.sin(phi + phi0)*sinTheta
        u_y =  np.cos(phi + phi0)*sinTheta
        u_z =  cosTheta*np.ones( np.shape(u_x) )

        v_x = -np.sin(phi + phi0)*cosTheta
        v_y =  np.cos(phi + phi0)*cosTheta
        v_z = -sinTheta*np.ones( np.shape(u_x) )

        # all grid points correspond
        self.Grid = {
            'u' : U, # the grid u coordinate values
            'v' : V, # the grid v coordinate values
            'X' : X, # the grid coordinates in general coordinates X
            'Y' : Y, # the grid coordinates in general coordinates Y
            'Z' : Z, # the grid coordinates in general coordinates Z
            'u_X' : u_x, # the grid local unit vector along u in gen coordinates, X
            'u_Y' : u_y, # the grid local unit vector along u in gen coordinates, Y
            'u_Z' : u_z,  # the grid local unit vector along u in gen coordinates, Z
            'v_X' : v_x, # the grid local unit vector along v in gen coordinates, X
            'v_Y' : v_y, # the grid local unit vector along v in gen coordinates, Y
            'v_Z' : v_z  # the grid local unit vector along v in gen coordinates, Z
            }

# ============================================================================


# ============================================================================
# GENERIC SURFACE POSITION SYSTEM (x,y,z) <-> (u,v)
class UVandXYZ :
    
    def __init__(self, Specimen):
        
        # retrieve Line position
        self.X_line = Specimen.ParallelToLine['X']
        self.Y_line = Specimen.ParallelToLine['Y']
        self.Z_line = Specimen.ParallelToLine['Z']
        self.Normal_X_line = Specimen.ParallelToLine['Normal_X']
        self.Normal_Y_line = Specimen.ParallelToLine['Normal_Y']
        self.Normal_Z_line = Specimen.ParallelToLine['Normal_Z']
        
        # local line parameters
        u = np.zeros( (len(self.X_line), ) )
        for p in np.arange(1, len(self.X_line)):
            u[p] = u[p-1] + np.sqrt( (self.X_line[p]-self.X_line[p-1])**2 
                                    + (self.Y_line[p]-self.Y_line[p-1])**2 
                                    + (self.Z_line[p]-self.Z_line[p-1])**2 )
        u = u - u[  Specimen.ParallelToLine['Origin_ind'] ]

        self.tck_u_x = splrep( u, self.X_line, k=3, s=0 )
        self.tck_u_y = splrep( u, self.Y_line, k=3, s=0 )
        self.tck_u_z = splrep( u, self.Z_line, k=3, s=0 )
        
        
        # the origin (u=0,v=0) 
        self.x_OriginPoint = splev(0, self.tck_u_x, der = 0)
        self.y_OriginPoint = splev(0, self.tck_u_y, der = 0)
        self.z_OriginPoint = splev(0, self.tck_u_z, der = 0)
     
        
    def get_PointAndVectorXYZFromU_single(self, u):

        self.single_x_alongU = float( splev(u, self.tck_u_x, der = 0) )
        self.single_y_alongU = float( splev(u, self.tck_u_y, der = 0) )
        self.single_z_alongU = float( splev(u, self.tck_u_z, der = 0) )        

        u_x = splev(u, self.tck_u_x, der = 1)
        u_y = splev(u, self.tck_u_y, der = 1)
        u_z = splev(u, self.tck_u_z, der = 1)
        
        Norm = np.sqrt( u_x**2 + u_y**2 + u_z**2 )
        self.single_u_x_alongU = float( u_x/Norm )
        self.single_u_y_alongU = float( u_y/Norm )
        self.single_u_z_alongU = float( u_z/Norm )


    def get_lineAlongV(self, Specimen, u):
        
        self.get_PointAndVectorXYZFromU_single(u)
        
        x0 = self.single_x_alongU
        y0 = self.single_y_alongU
        z0 = self.single_z_alongU
        u_x = self.single_u_x_alongU
        u_y = self.single_u_y_alongU
        u_z = self.single_u_z_alongU
        
        # check the orientation of the perp line data such that it is aligned in
        
        Specimen.get_LinePerpendicular( pos = [x0, y0, z0] )
        
        X_line_perp = Specimen.PerpendicularToLine['X']
        Y_line_perp = Specimen.PerpendicularToLine['Y']
        Z_line_perp = Specimen.PerpendicularToLine['Z']
        Normal_X_line_perp = Specimen.PerpendicularToLine['Normal_X']
        Normal_Y_line_perp = Specimen.PerpendicularToLine['Normal_Y']
        Normal_Z_line_perp = Specimen.PerpendicularToLine['Normal_Z']
        
        
        # find closest point in [X_line_perp, Y_line_perp, Z_line_perp] to [x_startPoint, y_startPoint, z_startPoint]
        # use the local approximation to arrange data along positive direction of v
        indv = np.argmin( (X_line_perp-x0)**2 + 
                         (Y_line_perp-y0)**2 + 
                         (Z_line_perp-z0)**2 )
        v_x_test = X_line_perp[indv+1]-X_line_perp[indv]
        v_y_test = Y_line_perp[indv+1]-Y_line_perp[indv]
        v_z_test = Z_line_perp[indv+1]-Z_line_perp[indv]
        
        # obtain an approximation for vector v based on closest point in parallel line data
        v_x_approx =  Normal_Y_line_perp[indv]*u_z - Normal_Z_line_perp[indv]*u_y
        v_y_approx = -Normal_X_line_perp[indv]*u_z + Normal_Z_line_perp[indv]*u_x
        v_z_approx =  Normal_X_line_perp[indv]*u_y - Normal_Y_line_perp[indv]*u_x
        
        
        if v_x_approx*v_x_test + v_y_approx*v_y_test + v_z_approx*v_z_test < 0 :
            self.X_line_perp = np.flip(X_line_perp)
            self.Y_line_perp = np.flip(Y_line_perp) 
            self.Z_line_perp = np.flip(Z_line_perp)
            self.Normal_X_line_perp = np.flip(Normal_X_line_perp)
            self.Normal_Y_line_perp = np.flip(Normal_Y_line_perp)
            self.Normal_Z_line_perp = np.flip(Normal_Z_line_perp)
            Specimen.PerpendicularToLine['Origin_ind'] = len(X_line_perp)-1 - Specimen.PerpendicularToLine['Origin_ind']
        else:
            self.X_line_perp = (X_line_perp)
            self.Y_line_perp = (Y_line_perp) 
            self.Z_line_perp = (Z_line_perp)
            self.Normal_X_line_perp = (Normal_X_line_perp)
            self.Normal_Y_line_perp = (Normal_Y_line_perp)
            self.Normal_Z_line_perp = (Normal_Z_line_perp)
        
        v_x_test = X_line_perp[indv+1]-X_line_perp[indv]
        v_y_test = Y_line_perp[indv+1]-Y_line_perp[indv]
        v_z_test = Z_line_perp[indv+1]-Z_line_perp[indv]
        Norm = np.sqrt( v_x_test**2 + v_y_test**2 + v_z_test**2 )
        
        v_x_test = v_x_test/Norm
        v_y_test = v_y_test/Norm
        v_z_test = v_z_test/Norm
        
    
        # arclength along v
        v = np.zeros( (len(X_line_perp), ) )
        for p in np.arange(1, len(X_line_perp)):
            v[p] = v[p-1] + np.sqrt( (self.X_line_perp[p]-self.X_line_perp[p-1])**2 
                                    + (self.Y_line_perp[p]-self.Y_line_perp[p-1])**2 
                                    + (self.Z_line_perp[p]-self.Z_line_perp[p-1])**2 )
        v = v - v[ Specimen.PerpendicularToLine['Origin_ind'] ]
        
        self.tck_Normal_v_x = splrep( v, self.Normal_X_line_perp, k=3, s=0 )
        self.tck_Normal_v_y = splrep( v, self.Normal_Y_line_perp, k=3, s=0 )
        self.tck_Normal_v_z = splrep( v, self.Normal_Z_line_perp, k=3, s=0 )
        
        self.tck_v_x = splrep( v, self.X_line_perp, k=3, s=0 )
        self.tck_v_y = splrep( v, self.Y_line_perp, k=3, s=0 )
        self.tck_v_z = splrep( v, self.Z_line_perp, k=3, s=0 )


    def get_PointAndVectorsXYZFromUV_single(self, Specimen, u, v):
        
        self.get_lineAlongV(Specimen, u)
        

        if hasattr(self, 'posUV_inXYZ'):
            
            # position at (u,v)
            self.posUV_inXYZ['x'].append( float( splev(v, self.tck_v_x, der = 0) ) )
            self.posUV_inXYZ['y'].append( float( splev(v, self.tck_v_y, der = 0) ) )
            self.posUV_inXYZ['z'].append( float( splev(v, self.tck_v_z, der = 0) ) )
            self.posUV_inXYZ['v'].append( v )
            self.posUV_inXYZ['u'].append( u )
            
            
            # v vector at (u, v)
            v_x = float( splev(v, self.tck_v_x, der = 1) )
            v_y = float( splev(v, self.tck_v_y, der = 1) )
            v_z = float( splev(v, self.tck_v_z, der = 1) )
            
            Norm = np.sqrt( v_x**2 + v_y**2 + v_z**2 )
            self.posUV_inXYZ['vec_v_x'].append( v_x/Norm )
            self.posUV_inXYZ['vec_v_y'].append( v_y/Norm )
            self.posUV_inXYZ['vec_v_z'].append( v_z/Norm )
            
            
            # normal vector at point (u, v)
            Normal_x = float( splev(v, self.tck_Normal_v_x, der = 0) )
            Normal_y = float( splev(v, self.tck_Normal_v_y, der = 0) )
            Normal_z = float( splev(v, self.tck_Normal_v_z, der = 0) )
            
            Norm = np.sqrt( Normal_x**2 + Normal_y**2 + Normal_z**2 )
            self.posUV_inXYZ['Normal_x'].append( Normal_x/Norm )
            self.posUV_inXYZ['Normal_y'].append( Normal_y/Norm )
            self.posUV_inXYZ['Normal_z'].append( Normal_z/Norm )
                
            
            # u vector at (u, v)
            self.posUV_inXYZ['vec_u_x'].append(  self.posUV_inXYZ['vec_v_y'][-1]*self.posUV_inXYZ['Normal_z'][-1] - self.posUV_inXYZ['vec_v_z'][-1]*self.posUV_inXYZ['Normal_y'][-1] )
            self.posUV_inXYZ['vec_u_y'].append( -self.posUV_inXYZ['vec_v_x'][-1]*self.posUV_inXYZ['Normal_z'][-1] + self.posUV_inXYZ['vec_v_z'][-1]*self.posUV_inXYZ['Normal_x'][-1] ) 
            self.posUV_inXYZ['vec_u_z'].append(  self.posUV_inXYZ['vec_v_x'][-1]*self.posUV_inXYZ['Normal_y'][-1] - self.posUV_inXYZ['vec_v_y'][-1]*self.posUV_inXYZ['Normal_x'][-1] )
            
        else:
            
            # position at (u,v)
            self.posUV_inXYZ = {'x' : [ float( splev(v, self.tck_v_x, der = 0) ) ],
                                'y' : [ float( splev(v, self.tck_v_y, der = 0) ) ],
                                'z' : [ float( splev(v, self.tck_v_z, der = 0) ) ],
                                'v' : [ v ],
                                'u' : [ u ] }
                        
            
            # v vector at (u, v)
            v_x = float( splev(v, self.tck_v_x, der = 1) )
            v_y = float( splev(v, self.tck_v_y, der = 1) )
            v_z = float( splev(v, self.tck_v_z, der = 1) )
            
            Norm = np.sqrt( v_x**2 + v_y**2 + v_z**2 )
            self.posUV_inXYZ['vec_v_x'] = [ v_x/Norm ]
            self.posUV_inXYZ['vec_v_y'] = [ v_y/Norm ]
            self.posUV_inXYZ['vec_v_z'] = [ v_z/Norm ]
            
            
            # normal vector at (u, v)
            Normal_x = float( splev(v, self.tck_Normal_v_x, der = 0) )
            Normal_y = float( splev(v, self.tck_Normal_v_y, der = 0) )
            Normal_z = float( splev(v, self.tck_Normal_v_z, der = 0) )
            
            Norm = np.sqrt( Normal_x**2 + Normal_y**2 + Normal_z**2 )
            self.posUV_inXYZ['Normal_x'] = [ Normal_x/Norm ]
            self.posUV_inXYZ['Normal_y'] = [ Normal_y/Norm ]
            self.posUV_inXYZ['Normal_z'] = [ Normal_z/Norm ]
            
            
            # u vector at (u, v)
            self.posUV_inXYZ['vec_u_x'] = [  self.posUV_inXYZ['vec_v_y'][-1]*self.posUV_inXYZ['Normal_z'][-1] - self.posUV_inXYZ['vec_v_z'][-1]*self.posUV_inXYZ['Normal_y'][-1] ]
            self.posUV_inXYZ['vec_u_y'] = [ -self.posUV_inXYZ['vec_v_x'][-1]*self.posUV_inXYZ['Normal_z'][-1] + self.posUV_inXYZ['vec_v_z'][-1]*self.posUV_inXYZ['Normal_x'][-1] ]
            self.posUV_inXYZ['vec_u_z'] = [  self.posUV_inXYZ['vec_v_x'][-1]*self.posUV_inXYZ['Normal_y'][-1] - self.posUV_inXYZ['vec_v_y'][-1]*self.posUV_inXYZ['Normal_x'][-1] ]
            
        
    def rotateOnSkewAngle( self, skewAngle ):


        self.posUV_inXYZ['RotSkew_vec_u_x'] = []
        self.posUV_inXYZ['RotSkew_vec_u_y'] = []
        self.posUV_inXYZ['RotSkew_vec_u_z'] = []
        
        self.posUV_inXYZ['RotSkew_vec_v_x'] = []
        self.posUV_inXYZ['RotSkew_vec_v_y'] = []
        self.posUV_inXYZ['RotSkew_vec_v_z'] = []
        
        
        if hasattr(skewAngle, "__len__"):
            self.posUV_inXYZ['SkewAngle'] = skewAngle
        else:
            skewAngle = [skewAngle]*len( self.posUV_inXYZ['x'] )
            self.posUV_inXYZ['SkewAngle'] = skewAngle
                        
            
        # skew rotated vector for 
        for p in range( len( self.posUV_inXYZ['x'] ) ):
            self.posUV_inXYZ['RotSkew_vec_u_x'].append(  np.cos(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_u_x'][p] + np.sin(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_v_x'][p] )
            self.posUV_inXYZ['RotSkew_vec_u_y'].append(  np.cos(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_u_y'][p] + np.sin(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_v_y'][p] )
            self.posUV_inXYZ['RotSkew_vec_u_z'].append(  np.cos(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_u_z'][p] + np.sin(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_v_z'][p] )
            
            self.posUV_inXYZ['RotSkew_vec_v_x'].append( -np.sin(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_u_x'][p] + np.cos(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_v_x'][p] )
            self.posUV_inXYZ['RotSkew_vec_v_y'].append( -np.sin(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_u_y'][p] + np.cos(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_v_y'][p] )
            self.posUV_inXYZ['RotSkew_vec_v_z'].append( -np.sin(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_u_z'][p] + np.cos(-np.pi/2 + skewAngle[p])*self.posUV_inXYZ['vec_v_z'][p] )
        
    
    
    
    def drawSpecimenLast2Positions( self, MySpecimen, prop ):
        
        fig = plt.figure(figsize = (15, 14))
        ax = fig.add_subplot(111, projection='3d')
        
        
        dispScale = np.max([np.max(np.abs(MySpecimen.Surface['X'])), np.max(np.abs(MySpecimen.Surface['Y'])), np.max(np.abs(MySpecimen.Surface['Z']))])
        ldisp = prop['Display']['characteristic_length_for_display']
        
        ax.plot_surface(MySpecimen.Surface['X'], MySpecimen.Surface['Y'], MySpecimen.Surface['Z'], color = 'gray', alpha=0.3)
        ax.scatter( MySpecimen.OriginLine['X'], MySpecimen.OriginLine['Y'], MySpecimen.OriginLine['Z'], c = 'orange', linewidth = 2)
        
        
        # plot the local u v vectors on the previous to last (u,v) position
        ax.plot(MySpecimen.ParallelToLine['X'], MySpecimen.ParallelToLine['Y'], MySpecimen.ParallelToLine['Z'], c = 'orange', linewidth = 2)
        ax.plot(self.X_line_perp, self.Y_line_perp, self.Z_line_perp, c ='lightgreen', linewidth = 2)
        ax.plot([self.posUV_inXYZ['x'][-2], self.posUV_inXYZ['x'][-2]+ldisp*self.posUV_inXYZ['vec_u_x'][-2]], 
                [self.posUV_inXYZ['y'][-2], self.posUV_inXYZ['y'][-2]+ldisp*self.posUV_inXYZ['vec_u_y'][-2]], 
                [self.posUV_inXYZ['z'][-2], self.posUV_inXYZ['z'][-2]+ldisp*self.posUV_inXYZ['vec_u_z'][-2]], c = 'red')
        ax.plot([self.posUV_inXYZ['x'][-2], self.posUV_inXYZ['x'][-2]+ldisp*self.posUV_inXYZ['vec_v_x'][-2]], 
                [self.posUV_inXYZ['y'][-2], self.posUV_inXYZ['y'][-2]+ldisp*self.posUV_inXYZ['vec_v_y'][-2]], 
                [self.posUV_inXYZ['z'][-2], self.posUV_inXYZ['z'][-2]+ldisp*self.posUV_inXYZ['vec_v_z'][-2]], c = 'green')
        ax.text( self.posUV_inXYZ['x'][-2]+ldisp*self.posUV_inXYZ['vec_u_x'][-2], 
                 self.posUV_inXYZ['y'][-2]+ldisp*self.posUV_inXYZ['vec_u_y'][-2], 
                 self.posUV_inXYZ['z'][-2]+ldisp*self.posUV_inXYZ['vec_u_z'][-2], 'u', c = 'red' )
        ax.text( self.posUV_inXYZ['x'][-2]+ldisp*self.posUV_inXYZ['vec_v_x'][-2],
                 self.posUV_inXYZ['y'][-2]+ldisp*self.posUV_inXYZ['vec_v_y'][-2], 
                 self.posUV_inXYZ['z'][-2]+ldisp*self.posUV_inXYZ['vec_v_z'][-2], 'v', c = 'green' )
                
        
        # plot the rotated frame at the last (u,v) position
        if 'RotSkew_vec_u_x' in self.posUV_inXYZ.keys():
            
            if self.posUV_inXYZ['SkewAngle'][-1] != 0:
                Ntheta = 100
                thetas = np.linspace(0, self.posUV_inXYZ['SkewAngle'][-1], Ntheta)
                s_x = self.posUV_inXYZ['x'][-1] + ldisp*(np.cos(thetas)*self.posUV_inXYZ['vec_u_x'][-1] + np.sin(thetas)*self.posUV_inXYZ['vec_v_x'][-1])
                s_y = self.posUV_inXYZ['y'][-1] + ldisp*(np.cos(thetas)*self.posUV_inXYZ['vec_u_y'][-1] + np.sin(thetas)*self.posUV_inXYZ['vec_v_y'][-1])
                s_z = self.posUV_inXYZ['z'][-1] + ldisp*(np.cos(thetas)*self.posUV_inXYZ['vec_u_z'][-1] + np.sin(thetas)*self.posUV_inXYZ['vec_v_z'][-1])
                ax.plot(s_x, s_y, s_z, c = 'blue', alpha = 0.5)
                my_str = 'skew angle : {angle:.0f}deg'
                ax.text(s_x[int(Ntheta/2)], s_y[int(Ntheta/2)], s_z[int(Ntheta/2)], my_str.format(angle = self.posUV_inXYZ['SkewAngle'][-1]/np.pi*180))
            
            ax.plot([self.posUV_inXYZ['x'][-1], self.posUV_inXYZ['x'][-1]+ldisp*self.posUV_inXYZ['vec_u_x'][-1]], 
                    [self.posUV_inXYZ['y'][-1], self.posUV_inXYZ['y'][-1]+ldisp*self.posUV_inXYZ['vec_u_y'][-1]], 
                    [self.posUV_inXYZ['z'][-1], self.posUV_inXYZ['z'][-1]+ldisp*self.posUV_inXYZ['vec_u_z'][-1]], c = 'red')
            
            ax.plot([self.posUV_inXYZ['x'][-1], self.posUV_inXYZ['x'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_u_x'][-1]], 
                    [self.posUV_inXYZ['y'][-1], self.posUV_inXYZ['y'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_u_y'][-1]], 
                    [self.posUV_inXYZ['z'][-1], self.posUV_inXYZ['z'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_u_z'][-1]], c = 'magenta')
            ax.plot([self.posUV_inXYZ['x'][-1], self.posUV_inXYZ['x'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_v_x'][-1]], 
                    [self.posUV_inXYZ['y'][-1], self.posUV_inXYZ['y'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_v_y'][-1]], 
                    [self.posUV_inXYZ['z'][-1], self.posUV_inXYZ['z'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_v_z'][-1]], c = 'blue')
            
            ax.text( self.posUV_inXYZ['x'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_u_x'][-1],  
                     self.posUV_inXYZ['y'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_u_y'][-1], 
                     self.posUV_inXYZ['z'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_u_z'][-1], "x'", c = 'magenta' )
            ax.text( self.posUV_inXYZ['x'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_v_x'][-1], 
                     self.posUV_inXYZ['y'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_v_y'][-1],
                     self.posUV_inXYZ['z'][-1] + ldisp*self.posUV_inXYZ['RotSkew_vec_v_z'][-1], "y'", c = 'blue' )
          
        
        my_str = '({u:.4f}, {v:.4f})'
        ax.text( self.posUV_inXYZ['x'][-2], self.posUV_inXYZ['y'][-2], self.posUV_inXYZ['z'][-2], my_str.format(u = self.posUV_inXYZ['u'][-1], v = 0) )
        ax.text( self.posUV_inXYZ['x'][-1], self.posUV_inXYZ['y'][-1], self.posUV_inXYZ['z'][-1], my_str.format(u = self.posUV_inXYZ['u'][-1], v = self.posUV_inXYZ['v'][-1]) )
        ax.text( MySpecimen.OriginLine['X'], MySpecimen.OriginLine['Y'], MySpecimen.OriginLine['Z'], '(0,0)' )


        ax.set_xlabel('x', fontsize = 20)
        ax.set_ylabel('y', fontsize = 20)
        ax.set_zlabel('z', fontsize = 20)

        ax.set_xlim(-dispScale, dispScale)
        ax.set_ylim(-dispScale, dispScale)
        ax.set_zlim(-dispScale, dispScale) 
        ax.set_box_aspect([1,1,1])
    
    
    def drawSpecimenGrid( self, MySpecimen, prop , DebugFigs = False):
        
        fig = plt.figure(figsize = (15, 14))
        ax = fig.add_subplot(111, projection='3d')
        ax.tick_params(axis='both', which='major', labelsize=15)
        
        dispScale = np.max([np.max(np.abs(MySpecimen.Surface['X'])), np.max(np.abs(MySpecimen.Surface['Y'])), np.max(np.abs(MySpecimen.Surface['Z']))])
        ldisp = prop['Display']['characteristic_length_for_display']
        
        ax.plot_surface(MySpecimen.Surface['X'], MySpecimen.Surface['Y'], MySpecimen.Surface['Z'], color = 'lightgray', alpha=0.3)
        ax.plot(MySpecimen.ParallelToLine['X'], MySpecimen.ParallelToLine['Y'], MySpecimen.ParallelToLine['Z'], c = 'orange', linewidth = 2)
        ax.scatter( MySpecimen.OriginLine['X'], MySpecimen.OriginLine['Y'], MySpecimen.OriginLine['Z'], c = 'orange', linewidth = 2)
        
        
        # plot the local u v grid
        for p in range(MySpecimen.Grid['u'].shape[0]):
            ax.plot(MySpecimen.Grid['X'][ p ], MySpecimen.Grid['Y'][ p ], MySpecimen.Grid['Z'][ p ], c = 'red') 

        for p in range(MySpecimen.Grid['u'].shape[1]):
            ax.plot(MySpecimen.Grid['X'][ :, p ], MySpecimen.Grid['Y'][ :, p ], MySpecimen.Grid['Z'][ :, p ], c = 'green')
        
        
        ax.set_xlabel('x', fontsize = 20)
        ax.set_ylabel('y', fontsize = 20)
        ax.set_zlabel('z', fontsize = 20)
        
        Xlims = ax.get_xlim()
        Ylims = ax.get_ylim()
        Zlims = ax.get_zlim()
        
        DX = Xlims[1]-Xlims[0]
        DY = Ylims[1]-Ylims[0]
        DZ = Zlims[1]-Zlims[0]
        MaxSize = np.max([DX, DY, DZ])
        
        ax.set_box_aspect([DX/MaxSize, DY/MaxSize, DZ/MaxSize])
        fig.tight_layout()
        
        
        
        
        fig = plt.figure(figsize = (13, 12))
        ax = fig.add_subplot( 111 )
        ax.set_xlabel("$u$", fontsize = 20)
        ax.set_ylabel("$v$", fontsize = 20)
        ax.tick_params(axis='both', which='major', labelsize=15)
        
        
        ind = np.where( np.abs(MySpecimen.Grid['v']) < 1e-14 )[0]
        u_ref = np.sort(np.unique(MySpecimen.Grid['u'][ ind ]))
        
        ax.plot( u_ref, np.zeros( (len(u_ref), )),  linewidth = 4, c = 'orange', alpha = 0.5)
        ax.scatter( MySpecimen.Grid['u'], MySpecimen.Grid['v'],  c = 'cornflowerblue', s = 8)
        

        ax.set_aspect('equal', adjustable='box')
        fig.tight_layout()
        
        

class scanPattern:

    def __init__(self, along_u, along_v):
        self.along_u = along_u
        self.along_v = along_v
        
        
    def GRID_rasterUV_stepsAlongU( self ):
        scan_u = np.asarray([])
        scan_v = np.asarray([])
        for p in range(len(self.along_u)):
            if p/2 == np.round(p/2): 
                scan_v = np.concatenate( (scan_v, self.along_v) )
            else:
                scan_v = np.concatenate( (scan_v, np.flip( self.along_v ) ) )
            scan_u = np.concatenate( (scan_u, np.full( (len(self.along_v), ), self.along_u[p])) )
        self.scan_u = scan_u
        self.scan_v = scan_v
        
    
    def GRID_rasterUV_stepsAlongV( self ):
        scan_u = np.asarray([])
        scan_v = np.asarray([])
        for p in range(len(self.along_v)):
            if p/2 == np.round(p/2): 
                scan_u = np.concatenate( (scan_u, self.along_u) )
            else:
                scan_u = np.concatenate( (scan_u, np.flip( self.along_u ) ) )
            scan_v = np.concatenate( (scan_v, np.full( (len(self.along_u), ), self.along_v[p])) )
        self.scan_u = scan_u
        self.scan_v = scan_v
        
    
    def GRID_combUV_stepsAlongV( self ):
        scan_u = np.asarray([])
        scan_v = np.asarray([])
        for p in range(len(self.along_v)): 
            scan_u = np.concatenate( (scan_u, self.along_u) )
            scan_u = np.concatenate( (scan_u, np.flip( self.along_u )) )
            scan_v = np.concatenate( (scan_v, np.full( (len(self.along_u), ), self.along_v[p])) )
            scan_v = np.concatenate( (scan_v, np.full( (len(self.along_u), ), self.along_v[p])) )
        self.scan_u = scan_u
        self.scan_v = scan_v
        
    
    def GRID_combUV_stepsAlongU( self ):
        scan_u = np.asarray([])
        scan_v = np.asarray([])
        for p in range(len(self.along_u)): 
            scan_v = np.concatenate( (scan_v, self.along_v) )
            scan_v = np.concatenate( (scan_v, np.flip( self.along_v )) )
            scan_u = np.concatenate( (scan_u, np.full( (len(self.along_v), ), self.along_u[p])) )
            scan_u = np.concatenate( (scan_u, np.full( (len(self.along_v), ), self.along_u[p])) )
        self.scan_u = scan_u
        self.scan_v = scan_v
        
        
    def LINE_UV( self ):
        self.scan_u = self.along_u
        self.scan_v = self.along_v
        

    def drawSpecimenScan( self, MySpecimen, MyUVandXYZ, prop ):
        
        for p in range(len(self.scan_u)):
            
            MyUVandXYZ.get_PointAndVectorsXYZFromUV_single(Specimen = MySpecimen, u = self.scan_u[p], v = self.scan_v[p])

        
        fig = plt.figure(figsize = (15, 14))
        ax = fig.add_subplot(111, projection='3d')
        
        
        dispScale = np.max([np.max(np.abs(MySpecimen.Surface['X'])), np.max(np.abs(MySpecimen.Surface['Y'])), np.max(np.abs(MySpecimen.Surface['Z']))])
        ldisp = prop['Display']['characteristic_length_for_display']
        
        
        # plot the surface and origin of u line
        ax.plot_surface(MySpecimen.Surface['X'], MySpecimen.Surface['Y'], MySpecimen.Surface['Z'], color = 'gray', alpha=0.3)
        ax.scatter( MySpecimen.OriginLine['X'], MySpecimen.OriginLine['Y'], MySpecimen.OriginLine['Z'], c = 'orange', linewidth = 2)
        
        
        # plot the local u line
        ax.plot(MySpecimen.ParallelToLine['X'], MySpecimen.ParallelToLine['Y'], MySpecimen.ParallelToLine['Z'], c = 'orange', linewidth = 2)
        
        # plot the local referential at the first point on the scan
        ax.plot([MyUVandXYZ.posUV_inXYZ['x'][0], MyUVandXYZ.posUV_inXYZ['x'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_u_x'][0]], 
                [MyUVandXYZ.posUV_inXYZ['y'][0], MyUVandXYZ.posUV_inXYZ['y'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_u_y'][0]], 
                [MyUVandXYZ.posUV_inXYZ['z'][0], MyUVandXYZ.posUV_inXYZ['z'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_u_z'][0]], c = 'red')
        ax.plot([MyUVandXYZ.posUV_inXYZ['x'][0], MyUVandXYZ.posUV_inXYZ['x'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_v_x'][0]], 
                [MyUVandXYZ.posUV_inXYZ['y'][0], MyUVandXYZ.posUV_inXYZ['y'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_v_y'][0]], 
                [MyUVandXYZ.posUV_inXYZ['z'][0], MyUVandXYZ.posUV_inXYZ['z'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_v_z'][0]], c = 'green')
        ax.text( MyUVandXYZ.posUV_inXYZ['x'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_u_x'][0], 
                 MyUVandXYZ.posUV_inXYZ['y'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_u_y'][0], 
                 MyUVandXYZ.posUV_inXYZ['z'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_u_z'][0], 'u', c = 'red' )
        ax.text( MyUVandXYZ.posUV_inXYZ['x'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_v_x'][0],
                 MyUVandXYZ.posUV_inXYZ['y'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_v_y'][0], 
                 MyUVandXYZ.posUV_inXYZ['z'][0]+ldisp*MyUVandXYZ.posUV_inXYZ['vec_v_z'][0], 'v', c = 'green' )

        # plot the trajectory
        ax.plot(MyUVandXYZ.posUV_inXYZ['x'], MyUVandXYZ.posUV_inXYZ['y'], MyUVandXYZ.posUV_inXYZ['z'], c = 'indianred', marker = '+')
    
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        ax.set_xlim(-dispScale, dispScale)
        ax.set_ylim(-dispScale, dispScale)
        ax.set_zlim(-dispScale, dispScale) 
        ax.set_box_aspect([1,1,1])
        
        
    
class SystemAxes:
    
    def __init__(self, prop):
        
        self.prop = prop
        
    
    def plot_probeElements(self, figs = False):
        # COD case
        if np.abs( self.prop['Wedge']['radiusAlongYw'] ) < np.inf :
            
            R0 = self.prop['Wedge']['radiusAlongYw']
            l0 = self.prop['Wedge']['length']
            w0 = self.prop['Wedge']['width']

            
            # Origin of wedge system
            Ox_w = 0
            Oy_w = 0
            Oz_w = 0
       
            
            # z_sw offset contact to specimen
            if np.abs( self.prop['Wedge']['radiusAlongYw'] ) > 10e32 :
                z_w0 = 0
            elif self.prop['Wedge']['radiusAlongYw'] > 0:
                z_w0 = -( np.abs(R0) - np.sqrt( R0**2 - (l0/2)**2 ) )
            else:
                z_w0 = 0
        
            alpha = self.prop['Wedge']['squintAngle']
            beta  = self.prop['Wedge']['wedgeAngle']
            gamma = self.prop['Wedge']['roofAngle']
            
            # the ref offsets
            x_ref = self.prop['Wedge']['NormalSecondaryOffset']
            y_ref = self.prop['Wedge']['NormalPrimaryOffset']
            z_ref = self.prop['Wedge']['NormalTertiaryOffset'] + z_w0
            
            
            
            
            # nominal primary axis unit vector 
            hat_yp_xw =  np.sin(alpha)*np.cos(beta)
            hat_yp_yw =  np.cos(alpha)*np.cos(beta)
            hat_yp_zw =  np.sin(beta)
            
            # nominal secondary unit vector, choose 
            hat_xp_xw = -np.cos(alpha)*np.cos(gamma) - np.sin(alpha)*np.sin(beta)*np.sin(gamma)
            hat_xp_yw = +np.sin(alpha)*np.cos(gamma) - np.cos(alpha)*np.sin(beta)*np.sin(gamma)
            hat_xp_zw =  np.cos(beta)*np.sin(gamma)
            
            # nominal normal vector
            hat_zp_xw = -np.cos(alpha)*np.sin(gamma) + np.sin(alpha)*np.sin(beta)*np.cos(gamma)
            hat_zp_yw = +np.sin(alpha)*np.sin(gamma) + np.cos(alpha)*np.sin(beta)*np.cos(gamma)
            hat_zp_zw = -np.cos(beta)*np.cos(gamma)
            
            
            
            
            if figs == True:
                
                dispScale = np.max([l0, w0])
                ldisp = 0.1*dispScale
                dispScale = 1.05*dispScale
                
                # ============================================================
                # plot system of axes in globality
                fig = plt.figure(figsize = (15, 14))
                ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho')
        
                # the primary, secondary, tertiary axes
                ax.plot([0, 0], [Oy_w, y_ref], [0, 0], c = 'green')
                ax.plot([Ox_w, x_ref], [y_ref, y_ref], [0, 0], c = 'green')
                ax.plot([x_ref, x_ref], [y_ref, y_ref], [Oz_w, z_ref], c = 'green')
                ax.text( 0, Oy_w + 0.5*(y_ref-Oy_w), 0 , "$y_{ref}$", c = 'green' , fontsize = 15)
                ax.text( Ox_w + 0.5*(x_ref-Ox_w), y_ref, 0 , "$x_{ref}$", c = 'green' , fontsize = 15)
                ax.text( x_ref, y_ref, Oz_w + 0.5*(z_ref-Oz_w) , "$z_{ref}$", c = 'green' , fontsize = 15)                
                
                # wedge referential
                ax.plot([Ox_w, Ox_w + ldisp*1], [0, 0], [0, 0], c = 'k')
                ax.plot([0, 0], [Oy_w, Oy_w + ldisp*1], [0, 0], c = 'k')
                ax.plot([0, 0], [0, 0], [Oz_w, Oz_w + ldisp*1], c = 'k')
                ax.text( Ox_w + ldisp*1,  0, 0, "$x_w$", c = 'k' , fontsize = 15)
                ax.text( 0,  Oy_w + ldisp*1, 0, "$y_w$", c = 'k' , fontsize = 15)
                ax.text( 0,  0, Oz_w + ldisp*1, "$z_w$", c = 'k' , fontsize = 15)
                ax.text( Ox_w, Oy_w, Oz_w, "$O_w$", c = 'k' , fontsize = 15) 
        
                # probe placement
                for p in range( self.prop['Probe']['Nel'] ):
                
                    Ly_p, Lx_p = np.meshgrid(np.linspace(-0.5*self.prop['Probe']['ElSidePrimary'][p]   , +0.5*self.prop['Probe']['ElSidePrimary'][p], 3), 
                                             np.linspace(-0.5*self.prop['Probe']['ElSideSecondary'][p] , +0.5*self.prop['Probe']['ElSideSecondary'][p], 3)) 
                    XX = ( self.prop['Probe']['ElPos_SecondaryAxis'][p] - self.prop['Probe']['Ref_SecondaryAxis'] + Lx_p )*hat_xp_xw + ( self.prop['Probe']['ElPos_PrimaryAxis'][p] - self.prop['Probe']['Ref_PrimaryAxis'] + Ly_p)*hat_yp_xw + x_ref
                    YY = ( self.prop['Probe']['ElPos_SecondaryAxis'][p] - self.prop['Probe']['Ref_SecondaryAxis'] + Lx_p )*hat_xp_yw + ( self.prop['Probe']['ElPos_PrimaryAxis'][p] - self.prop['Probe']['Ref_PrimaryAxis'] + Ly_p)*hat_yp_yw + y_ref
                    ZZ = ( self.prop['Probe']['ElPos_SecondaryAxis'][p] - self.prop['Probe']['Ref_SecondaryAxis'] + Lx_p )*hat_xp_zw + ( self.prop['Probe']['ElPos_PrimaryAxis'][p] - self.prop['Probe']['Ref_PrimaryAxis'] + Ly_p)*hat_yp_zw + z_ref
                    ax.plot_surface(XX, YY, ZZ, color = "cornflowerblue", alpha = 0.5, antialiased = True)


                # plot the yp axis 
                ldisp_tmp = np.max([ldisp, ( 1.05*np.max(self.prop['Probe']['ElPos_PrimaryAxis'] + 0.5*self.prop['Probe']['ElSidePrimary']) - self.prop['Probe']['Ref_PrimaryAxis'] ) ])
                
                ax.plot( [x_ref, x_ref + ldisp_tmp*hat_yp_xw ],
                         [y_ref, y_ref + ldisp_tmp*hat_yp_yw ], 
                         [z_ref, z_ref + ldisp_tmp*hat_yp_zw ],
                         c = 'gray')
                ax.text(x_ref + ldisp_tmp*hat_yp_xw, 
                        y_ref + ldisp_tmp*hat_yp_yw, 
                        z_ref + ldisp_tmp*hat_yp_zw,
                        "$y_p$", c = 'gray', fontsize = 15)
                
                # plot the xp axis 
                ldisp_tmp = np.max([ldisp, ( 1.05*np.max(self.prop['Probe']['ElPos_SecondaryAxis'] + 0.5*self.prop['Probe']['ElSideSecondary']) - self.prop['Probe']['Ref_SecondaryAxis'] ) ])
                
                ax.plot( [x_ref, x_ref + ldisp_tmp*hat_xp_xw ],
                         [y_ref, y_ref + ldisp_tmp*hat_xp_yw ], 
                         [z_ref, z_ref + ldisp_tmp*hat_xp_zw ],
                         c = 'gray')
                ax.text(x_ref +  ldisp_tmp *hat_xp_xw, 
                        y_ref +  ldisp_tmp *hat_xp_yw, 
                        z_ref +  ldisp_tmp *hat_xp_zw,
                        "$x_p$", c = 'gray', fontsize = 15)
                


                # the projection of the hat_y_p vector in (x_w, y_w) plane and the wedge angle
                hat_yp_proj_xw = hat_yp_xw
                hat_yp_proj_yw = hat_yp_yw
                hat_yp_proj_zw = 0
                
                Norm = np.sqrt( hat_yp_proj_xw**2 + hat_yp_proj_yw**2 + hat_yp_proj_zw**2 )
                hat_yp_proj_xw = hat_yp_proj_xw/Norm
                hat_yp_proj_yw = hat_yp_proj_yw/Norm
                hat_yp_proj_zw = hat_yp_proj_zw/Norm
                
                Ntheta = 100
                
                # plot the arc of the wedge angle
                if beta != 0:
                    thetas = np.linspace(0, beta, Ntheta)
                        
                    s_xw = ldisp*np.cos(thetas)*hat_yp_proj_xw
                    s_yw = ldisp*np.cos(thetas)*hat_yp_proj_yw
                    s_zw = ldisp*np.sin(thetas)*1
                    ax.plot([x_ref, x_ref + ldisp*hat_yp_proj_xw], [y_ref, y_ref + ldisp*hat_yp_proj_yw], [z_ref, z_ref + hat_yp_proj_zw], c = 'indianred', alpha = 0.5)
                    ax.plot(x_ref + s_xw, y_ref + s_yw, z_ref + s_zw, c='indianred', alpha = 0.5, label = '$\\beta$ : wedge angle ; $\\sin\\beta = \\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{z}}_w$')
                    ax.text(x_ref + s_xw[int(Ntheta/2)], y_ref + s_yw[int(Ntheta/2)], z_ref + s_zw[int(Ntheta/2)], "$\\beta$", c = 'indianred', fontsize = 15)
                    
                
                # plot the arc of the squint angle
                if alpha != 0:
                    thetas = np.linspace(0, alpha, Ntheta)
                        
                    s_xw = ldisp*np.sin(thetas)
                    s_yw = ldisp*np.cos(thetas)
                    s_zw = ldisp*np.sin(thetas)*0
                    ax.plot([x_ref, x_ref ], [y_ref, y_ref + ldisp], [z_ref, z_ref ], c = 'limegreen', alpha=0.5)
                    ax.plot(x_ref + s_xw, y_ref + s_yw, z_ref + s_zw, c='limegreen', alpha= 0.5, label = '$\\alpha$ : squint angle ; $\\sin\\alpha = \\frac{ \\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{x}}_w }{\\sqrt{ (\\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{x}}_w)^2 + (\\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{y}}_w)^2 }}$')
                    ax.text(x_ref + s_xw[int(Ntheta/2)], y_ref + s_yw[int(Ntheta/2)], z_ref + s_zw[int(Ntheta/2)], "$\\alpha$", c = 'limegreen', fontsize = 15)
                    
                
                # the projection of the hat_x_p vector in (x_w, y_w) plane and the roof angle
                hat_xp_proj_xw = np.cos(gamma)*hat_xp_xw + np.sin(gamma)*hat_zp_xw
                hat_xp_proj_yw = np.cos(gamma)*hat_xp_yw + np.sin(gamma)*hat_zp_yw
                hat_xp_proj_zw = 0
                
                Norm = np.sqrt( hat_xp_proj_xw**2 + hat_xp_proj_yw**2 + hat_xp_proj_zw**2 )
                hat_xp_proj_xw = hat_xp_proj_xw/Norm
                hat_xp_proj_yw = hat_xp_proj_yw/Norm
                hat_xp_proj_zw = hat_xp_proj_zw/Norm
                
                # plot the arc of the roof angle
                if gamma != 0:
                    thetas = np.linspace(0, gamma, Ntheta)
                        
                    s_xw = ldisp*(np.cos(thetas)*hat_xp_xw + np.sin(thetas)*hat_zp_xw)
                    s_yw = ldisp*(np.cos(thetas)*hat_xp_yw + np.sin(thetas)*hat_zp_yw)
                    s_zw = ldisp*(np.cos(thetas)*hat_xp_zw + np.sin(thetas)*hat_zp_zw)
                    ax.plot([x_ref, x_ref + ldisp*hat_xp_proj_xw], [y_ref, y_ref + ldisp*hat_xp_proj_yw], [z_ref, z_ref + hat_xp_proj_zw], c = 'mediumpurple', alpha=0.5)
                    ax.plot(x_ref + s_xw, y_ref + s_yw, z_ref + s_zw, c='mediumpurple', alpha= 0.5, label = '$\\gamma$ : roof angle ; $\\sin\\gamma = \\frac{\\hat{\mathbf{x}}_p \\cdot \\hat{\\mathbf{z}}_w}{\\cos\\beta}$')
                    ax.text(x_ref + s_xw[int(Ntheta/2)], y_ref + s_yw[int(Ntheta/2)], z_ref + s_zw[int(Ntheta/2)], "$\\gamma$", c = 'mediumpurple', fontsize = 15)
                    



                Xlims = ax.get_xlim()
                Ylims = ax.get_ylim()
                Zlims = ax.get_zlim()
                
                DX = Xlims[1]-Xlims[0]
                DY = Ylims[1]-Ylims[0]
                DZ = Zlims[1]-Zlims[0]
                MaxSize = np.max([DX, DY, DZ])
                
                ax.set_box_aspect([DX/MaxSize, DY/MaxSize, DZ/MaxSize])
                fig.tight_layout()
                ax.view_init(azim=-90, elev=90)
        
    def get_wedgePositionAndAxes(self, figs = False):
        
        # COD case
        if np.abs( self.prop['Wedge']['radiusAlongYw'] ) < np.inf :
            
            R0 = self.prop['Wedge']['radiusAlongYw']
            l0 = self.prop['Wedge']['length']
            w0 = self.prop['Wedge']['width']

            if R0**2-(l0/2)**2 < 0:
                print("Exit early R0**2-(l0/2)**2 < 0")
                return
            

            # Origin of wedge system
            x_w0 = 0
            y_w0 = 0
            z_w0 = np.max( [0, R0-np.sqrt(R0**2-(l0/2)**2)  ] )
       
            # angular opening of wedge bottom
            theta = np.arctan( (l0/2)/np.sqrt( R0**2 - (l0/2)**2 ) )
            
            # z_sw offset contact to specimen
            if np.abs( self.prop['Wedge']['radiusAlongYw'] ) > 10e32 :
                z_sw = 0
            elif self.prop['Wedge']['radiusAlongYw'] > 0:
                z_sw = -( np.abs(R0) - np.sqrt( R0**2 - (l0/2)**2 ) )
            else:
                z_sw = 0
        
            alpha = self.prop['Wedge']['squintAngle']
            beta  = self.prop['Wedge']['wedgeAngle']
            gamma = self.prop['Wedge']['roofAngle']


            # the nominal offsets
            x_pw = self.prop['Wedge']['NormalSecondaryOffset']
            y_pw = self.prop['Wedge']['NormalPrimaryOffset']
            z_pw = self.prop['Wedge']['NormalTertiaryOffset']
            
            
            # nominal primary axis unit vector 
            hat_yp_xw =  np.sin(alpha)*np.cos(beta)
            hat_yp_yw =  np.cos(alpha)*np.cos(beta)
            hat_yp_zw =  np.sin(beta)
            
            # nominal secondary unit vector, choose 
            hat_xp_xw = -np.cos(alpha)*np.cos(gamma) - np.sin(alpha)*np.sin(beta)*np.sin(gamma)
            hat_xp_yw = +np.sin(alpha)*np.cos(gamma) - np.cos(alpha)*np.sin(beta)*np.sin(gamma)
            hat_xp_zw =  np.cos(beta)*np.sin(gamma)
            
            # nominal normal vector
            hat_zp_xw = -np.cos(alpha)*np.sin(gamma) + np.sin(alpha)*np.sin(beta)*np.cos(gamma)
            hat_zp_yw = +np.sin(alpha)*np.sin(gamma) + np.cos(alpha)*np.sin(beta)*np.cos(gamma)
            hat_zp_zw = -np.cos(beta)*np.cos(gamma)
            
            
            # center of curvature, and z' vector    
            radCenter_xw = x_w0
            radCenter_yw = -l0/2
            radCenter_zw = -R0 + ( R0-np.sign(R0)*np.sqrt(R0**2-(l0/2)**2) )
            
            # vector along part tangent x_prime
            x_prime_xw = 1
            x_prime_yw = 0
            x_prime_zw = 0

            
            # vector along part normal z_prime
            z_prime_xw = radCenter_xw - x_w0
            z_prime_yw = radCenter_yw - y_w0
            z_prime_zw = -(radCenter_zw) - z_w0

            Norm = np.sqrt(z_prime_xw**2 + z_prime_yw**2 + z_prime_zw**2)
            z_prime_xw = z_prime_xw/Norm
            z_prime_yw = z_prime_yw/Norm
            z_prime_zw = z_prime_zw/Norm
            
            
            # vector tangent to part surface y_prime
            y_prime_xw =  z_prime_yw*x_prime_zw - z_prime_zw*x_prime_yw
            y_prime_yw = -z_prime_xw*x_prime_zw + z_prime_zw*x_prime_xw
            y_prime_zw =  z_prime_xw*x_prime_yw - z_prime_yw*x_prime_xw
            
            
            if figs == True:
                
                dispScale = np.max([l0, w0])
                ldisp = 0.1*dispScale
                dispScale = 1.05*dispScale
                
                # ============================================================
                # plot system of axes in globality
                fig = plt.figure(figsize = (15, 14))
                ax = fig.add_subplot(111, projection='3d')
                
                
                
                
                
                # the curvature along wedge contact with part
                Ntheta = 100
                thetas = np.linspace(-theta, theta, Ntheta)
                rad_xw = np.full( (Ntheta,), 0)
                rad_yw = (R0)*np.sin(thetas) + radCenter_yw
                rad_zw = (R0)*np.cos(thetas) + radCenter_zw
                ax.plot(np.concatenate( (rad_xw, np.asarray([0, 0, 0])) ), 
                        np.concatenate( (rad_yw, np.asarray([rad_yw[-1], rad_yw[0], rad_yw[0]]) ) ), 
                        np.concatenate( (rad_zw, np.asarray([0, 0, rad_zw[0]])) ) , c = 'gold', alpha = 0.3)
                
                # the primary, secondary, tertiary axes
                ax.plot([0, 0], [y_w0, y_pw], [0, 0], c = 'orange')
                ax.plot([x_w0, x_pw], [y_pw, y_pw], [0, 0], c = 'orange')
                ax.plot([x_pw, x_pw], [y_pw, y_pw], [z_w0, z_pw], c = 'orange')
                ax.text( 0, y_w0 + 0.5*(y_pw-y_w0), 0 , "primary offset", c = 'orange' )
                ax.text( x_w0 + 0.5*(x_pw-x_w0), y_pw, 0 , "secondary offset", c = 'orange' )
                ax.text( x_pw, y_pw, z_w0 + 0.5*(z_pw-z_w0) , "tertiary offset", c = 'orange' )                
                
                # wedge referential
                ax.plot([x_w0, x_w0 + ldisp*1], [0, 0], [0, 0], c = 'k')
                ax.plot([0, 0], [y_w0, y_w0 + ldisp*1], [0, 0], c = 'k')
                ax.plot([0, 0], [0, 0], [0, 0 + ldisp*1], c = 'k')
                ax.text( x_w0 + ldisp*1,  0, 0, "$x_w$", c = 'k' )
                ax.text( 0,  y_w0 + ldisp*1, 0, "$y_w$", c = 'k' )
                ax.text( 0,  0,  ldisp*1, "$z_w$", c = 'k' )  
                
                # the local surface referential (rotated part referential)
                # ax.plot([x_w0, x_w0 + ldisp*x_prime_xw], [y_w0, y_w0 + ldisp*x_prime_yw], [z_w0 + z_sw, z_w0 + z_sw + ldisp*x_prime_zw], c = 'gray')
                # ax.plot([x_w0, x_w0 + ldisp*y_prime_xw], [y_w0, y_w0 + ldisp*y_prime_yw], [z_w0 + z_sw, z_w0 + z_sw + ldisp*y_prime_zw], c = 'gray')
                # ax.plot([x_w0, x_w0 + ldisp*z_prime_xw], [y_w0, y_w0 + ldisp*z_prime_yw], [z_w0 + z_sw, z_w0 + z_sw + ldisp*z_prime_zw], c = 'gray')
                # ax.text( x_w0 + ldisp*x_prime_xw,  y_w0 + ldisp*x_prime_yw, z_w0 + z_sw + ldisp*x_prime_zw, "$x'$", c = 'gray' )                  
                # ax.text( x_w0 + ldisp*y_prime_xw,  y_w0 + ldisp*y_prime_yw, z_w0 + z_sw + ldisp*y_prime_zw, "$y'$", c = 'gray' )
                # ax.text( x_w0 + ldisp*z_prime_xw,  y_w0 + ldisp*z_prime_yw, z_w0 + z_sw + ldisp*z_prime_zw, "$z'$", c = 'gray' ) 
                
                # the probe normal placement
                ax.plot([x_pw, x_pw + ldisp*hat_xp_xw], [y_pw, y_pw + ldisp*hat_xp_yw], [z_pw, z_pw + ldisp*hat_xp_zw], c = 'cornflowerblue')
                ax.plot([x_pw, x_pw + ldisp*hat_yp_xw], [y_pw, y_pw + ldisp*hat_yp_yw], [z_pw, z_pw + ldisp*hat_yp_zw], c = 'cornflowerblue')
                ax.plot([x_pw, x_pw + ldisp*hat_zp_xw], [y_pw, y_pw + ldisp*hat_zp_yw], [z_pw, z_pw + ldisp*hat_zp_zw], c = 'cornflowerblue')
                ax.text( x_pw + ldisp*hat_xp_xw,  y_pw + ldisp*hat_xp_yw, z_pw + ldisp*hat_xp_zw, "$x_p$", c = 'cornflowerblue' ) 
                ax.text( x_pw + ldisp*hat_yp_xw,  y_pw + ldisp*hat_yp_yw, z_pw + ldisp*hat_yp_zw, "$y_p$", c = 'cornflowerblue' )
                ax.text( x_pw + ldisp*hat_zp_xw,  y_pw + ldisp*hat_zp_yw, z_pw + ldisp*hat_zp_zw, "$z_p$", c = 'cornflowerblue' )                
                
                
                
                # the projection of the hat_y_p vector in (x_w, y_w) plane and the wedge angle
                hat_yp_proj_xw = hat_yp_xw
                hat_yp_proj_yw = hat_yp_yw
                hat_yp_proj_zw = 0
                
                Norm = np.sqrt( hat_yp_proj_xw**2 + hat_yp_proj_yw**2 + hat_yp_proj_zw**2 )
                hat_yp_proj_xw = hat_yp_proj_xw/Norm
                hat_yp_proj_yw = hat_yp_proj_yw/Norm
                hat_yp_proj_zw = hat_yp_proj_zw/Norm
                
                if beta != 0:
                    thetas = np.linspace(0, beta, Ntheta)
                        
                    s_xw = ldisp*np.cos(thetas)*hat_yp_proj_xw
                    s_yw = ldisp*np.cos(thetas)*hat_yp_proj_yw
                    s_zw = ldisp*np.sin(thetas)*1
                    ax.plot([x_pw, x_pw + ldisp*hat_yp_proj_xw], [y_pw, y_pw + ldisp*hat_yp_proj_yw], [z_pw, z_pw + hat_yp_proj_zw], c = 'indianred', alpha = 0.5)
                    ax.plot(x_pw + s_xw, y_pw + s_yw, z_pw + s_zw, c='indianred', alpha = 0.5, label = '$\\beta$ : wedge angle ; $\\sin\\beta = \\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{z}}_w$')
                    ax.text(x_pw + s_xw[int(Ntheta/2)], y_pw + s_yw[int(Ntheta/2)], z_pw + s_zw[int(Ntheta/2)], "$\\beta$", c = 'indianred')
                    print('wedge angle:')
                    print('   sin beta_nominal   = ' + str(np.sin(beta)))
                    print('   sin beta_empirical = ' + str( hat_yp_zw*1 ))
                
                
                if alpha != 0:
                    thetas = np.linspace(0, alpha, Ntheta)
                        
                    s_xw = ldisp*np.sin(thetas)
                    s_yw = ldisp*np.cos(thetas)
                    s_zw = ldisp*np.sin(thetas)*0
                    ax.plot([x_pw, x_pw ], [y_pw, y_pw + ldisp], [z_pw, z_pw ], c = 'limegreen', alpha=0.5)
                    ax.plot(x_pw + s_xw, y_pw + s_yw, z_pw + s_zw, c='limegreen', alpha= 0.5, label = '$\\alpha$ : squint angle ; $\\sin\\alpha = \\frac{ \\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{x}}_w }{\\sqrt{ (\\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{x}}_w)^2 + (\\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{y}}_w)^2 }}$')
                    ax.text(x_pw + s_xw[int(Ntheta/2)], y_pw + s_yw[int(Ntheta/2)], z_pw + s_zw[int(Ntheta/2)], "$\\alpha$", c = 'limegreen')
                    print('squint angle:')
                    print('   sin alpha_nominal   = ' + str(np.sin(alpha)))
                    print('   sin alpha_empirical = ' + str( hat_yp_xw*1/np.sqrt(hat_yp_xw**2+hat_yp_yw**2) ))
                
                
                # the projection of the hat_x_p vector in (x_w, y_w) plane and the roof angle
                hat_xp_proj_xw = np.cos(gamma)*hat_xp_xw + np.sin(gamma)*hat_zp_xw
                hat_xp_proj_yw = np.cos(gamma)*hat_xp_yw + np.sin(gamma)*hat_zp_yw
                hat_xp_proj_zw = 0
                
                Norm = np.sqrt( hat_xp_proj_xw**2 + hat_xp_proj_yw**2 + hat_xp_proj_zw**2 )
                hat_xp_proj_xw = hat_xp_proj_xw/Norm
                hat_xp_proj_yw = hat_xp_proj_yw/Norm
                hat_xp_proj_zw = hat_xp_proj_zw/Norm
                
                if gamma != 0:
                    thetas = np.linspace(0, gamma, Ntheta)
                        
                    s_xw = ldisp*(np.cos(thetas)*hat_xp_xw + np.sin(thetas)*hat_zp_xw)
                    s_yw = ldisp*(np.cos(thetas)*hat_xp_yw + np.sin(thetas)*hat_zp_yw)
                    s_zw = ldisp*(np.cos(thetas)*hat_xp_zw + np.sin(thetas)*hat_zp_zw)
                    ax.plot([x_pw, x_pw + ldisp*hat_xp_proj_xw], [y_pw, y_pw + ldisp*hat_xp_proj_yw], [z_pw, z_pw + hat_xp_proj_zw], c = 'mediumpurple', alpha=0.5)
                    ax.plot(x_pw + s_xw, y_pw + s_yw, z_pw + s_zw, c='mediumpurple', alpha= 0.5, label = '$\\gamma$ : roof angle ; $\\sin\\gamma = \\frac{\\hat{\mathbf{x}}_p \\cdot \\hat{\\mathbf{z}}_w}{\\cos\\beta}$')
                    ax.text(x_pw + s_xw[int(Ntheta/2)], y_pw + s_yw[int(Ntheta/2)], z_pw + s_zw[int(Ntheta/2)], "$\\gamma$", c = 'mediumpurple')
                    print('roof angle:')
                    print('   sin gamma_nominal   = ' + str(np.sin(gamma)))
                    print('   sin gamma_empirical = ' + str( hat_xp_zw*1/np.cos(beta) ))
                
                
                print("test orthogonality:")
                print("   hatx_p.haty_p = " + str(hat_xp_xw*hat_yp_xw + hat_xp_yw*hat_yp_yw + hat_xp_zw*hat_yp_zw))
                print("   hatx_p.hatz_p = " + str(hat_xp_xw*hat_zp_xw + hat_xp_yw*hat_zp_yw + hat_xp_zw*hat_zp_zw))
                print("   haty_p.hatz_p = " + str(hat_yp_xw*hat_zp_xw + hat_yp_yw*hat_zp_yw + hat_yp_zw*hat_zp_zw))
                
                
                
                ax.set_xlabel("$x_w'$")
                ax.set_ylabel("$y_w'$")
                ax.set_zlabel("$z_w'$")
                
                ax.legend()
                
                ax.set_xlim(-dispScale, dispScale)
                ax.set_ylim(-dispScale, dispScale)
                ax.set_zlim(-dispScale, dispScale) 
                ax.set_box_aspect([1,1,1])
                # ============================================================                
                                
                
                # ============================================================
                # zoom in nominal probe axes in wedge reference frame 
                fig = plt.figure(figsize = (15, 14))
                ax = fig.add_subplot(111, projection='3d')
                
                
                # the plane perpendicular to y_p
                Lx_p, Lz_p = np.meshgrid(np.linspace(-ldisp, ldisp, 5), np.linspace(-ldisp, ldisp, 5)) 
                XX = Lx_p*hat_xp_xw + Lz_p*hat_zp_xw + x_pw
                YY = Lx_p*hat_xp_yw + Lz_p*hat_zp_yw + y_pw
                ZZ = Lx_p*hat_xp_zw + Lz_p*hat_zp_zw + z_pw
                ax.plot_surface(XX, YY, ZZ, color = "cornflowerblue", alpha = 0.1)
                ax.plot_wireframe(XX, YY, ZZ, color = "cornflowerblue", alpha = 0.2)
                
                
                
                # wedge referential at reference point
                ax.plot([x_pw, x_pw + ldisp*1], [y_pw, y_pw], [z_pw, z_pw], c = 'gold')
                ax.plot([x_pw, x_pw], [y_pw, y_pw + ldisp*1], [z_pw, z_pw], c = 'gold')
                ax.plot([x_pw, x_pw], [y_pw, y_pw], [z_pw, z_pw + ldisp*1], c = 'gold')
                ax.text( x_pw + ldisp*1,  y_pw, z_pw, "$x_w'$", c = 'gold' , fontsize = 20)
                ax.text( x_pw,  y_pw + ldisp*1, z_pw, "$y_w'$", c = 'gold' , fontsize = 20)
                ax.text( x_pw,  y_pw, z_pw + ldisp*1, "$z_w'$", c = 'gold' , fontsize = 20)  
                
                
                # the probe nominal normal placement
                ax.plot([x_pw, x_pw + ldisp*hat_xp_xw], [y_pw, y_pw + ldisp*hat_xp_yw], [z_pw, z_pw + ldisp*hat_xp_zw], c = 'cornflowerblue', label = '$x_p$ : nominal secondary axis')
                ax.plot([x_pw, x_pw + ldisp*hat_yp_xw], [y_pw, y_pw + ldisp*hat_yp_yw], [z_pw, z_pw + ldisp*hat_yp_zw], c = 'cornflowerblue', label = '$y_p$ : nominal primary axis')
                ax.plot([x_pw, x_pw + ldisp*hat_zp_xw], [y_pw, y_pw + ldisp*hat_zp_yw], [z_pw, z_pw + ldisp*hat_zp_zw], c = 'cornflowerblue', label = '$z_p$ : normal axis')
                ax.text( x_pw + ldisp*hat_xp_xw,  y_pw + ldisp*hat_xp_yw, z_pw + ldisp*hat_xp_zw, "$x_p$", c = 'cornflowerblue' , fontsize = 20) 
                ax.text( x_pw + ldisp*hat_yp_xw,  y_pw + ldisp*hat_yp_yw, z_pw + ldisp*hat_yp_zw, "$y_p$", c = 'cornflowerblue' , fontsize = 20)
                ax.text( x_pw + ldisp*hat_zp_xw,  y_pw + ldisp*hat_zp_yw, z_pw + ldisp*hat_zp_zw, "$z_p$", c = 'cornflowerblue' , fontsize = 20) 
                
                
                # the projection of the hat_y_p vector in (x_w, y_w) plane and the wedge angle
                hat_yp_proj_xw = hat_yp_xw
                hat_yp_proj_yw = hat_yp_yw
                hat_yp_proj_zw = 0
                
                Norm = np.sqrt( hat_yp_proj_xw**2 + hat_yp_proj_yw**2 + hat_yp_proj_zw**2 )
                hat_yp_proj_xw = hat_yp_proj_xw/Norm
                hat_yp_proj_yw = hat_yp_proj_yw/Norm
                hat_yp_proj_zw = hat_yp_proj_zw/Norm
                
                if beta != 0:
                    thetas = np.linspace(0, beta, Ntheta)
                        
                    s_xw = ldisp*np.cos(thetas)*hat_yp_proj_xw
                    s_yw = ldisp*np.cos(thetas)*hat_yp_proj_yw
                    s_zw = ldisp*np.sin(thetas)*1
                    ax.plot([x_pw, x_pw + ldisp*hat_yp_proj_xw], [y_pw, y_pw + ldisp*hat_yp_proj_yw], [z_pw, z_pw + hat_yp_proj_zw], c = 'indianred', alpha = 0.5)
                    ax.plot(x_pw + s_xw, y_pw + s_yw, z_pw + s_zw, c='indianred', alpha = 0.5, label = '$\\beta$ : wedge angle ; $\\sin\\beta = \\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{z}}_w$')
                    ax.text(x_pw + s_xw[int(Ntheta/2)], y_pw + s_yw[int(Ntheta/2)], z_pw + s_zw[int(Ntheta/2)], "$\\beta$", c = 'indianred', fontsize = 20)
                    
                
                
                # the projection of the hat_y_p vector in (x_w, y_w) plane, then its projection along the y_w axis, and the squint angle
                hat_yp_projproj_xw = 0
                hat_yp_projproj_yw = hat_yp_yw
                hat_yp_projproj_zw = 0
                
                Norm = np.sqrt( hat_yp_projproj_xw**2 + hat_yp_projproj_yw**2 + hat_yp_projproj_zw**2 )
                hat_yp_projproj_xw = hat_yp_projproj_xw/Norm
                hat_yp_projproj_yw = hat_yp_projproj_yw/Norm
                hat_yp_projproj_zw = hat_yp_projproj_zw/Norm
                
                if alpha != 0:
                    thetas = np.linspace(0, alpha, Ntheta)
                        
                    s_xw = ldisp*np.sin(thetas)
                    s_yw = ldisp*np.cos(thetas)
                    s_zw = ldisp*np.sin(thetas)*0
                    ax.plot([x_pw, x_pw ], [y_pw, y_pw + ldisp], [z_pw, z_pw ], c = 'limegreen', alpha=0.5)
                    ax.plot(x_pw + s_xw, y_pw + s_yw, z_pw + s_zw, c='limegreen', alpha= 0.5, label = '$\\alpha$ : squint angle ; $\\sin\\alpha = \\frac{ \\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{x}}_w }{\\sqrt{ (\\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{x}}_w)^2 + (\\hat{\mathbf{y}}_p \\cdot \hat{\\mathbf{y}}_w )^2}}$')
                    ax.text(x_pw + s_xw[int(Ntheta/2)], y_pw + s_yw[int(Ntheta/2)], z_pw + s_zw[int(Ntheta/2)], "$\\alpha$", c = 'limegreen', fontsize = 20)
                
                
                
                
                # the projection of the hat_x_p vector in (x_w, y_w) plane and the roof angle
                hat_xp_proj_xw = np.cos(gamma)*hat_xp_xw + np.sin(gamma)*hat_zp_xw
                hat_xp_proj_yw = np.cos(gamma)*hat_xp_yw + np.sin(gamma)*hat_zp_yw
                hat_xp_proj_zw = 0
                
                Norm = np.sqrt( hat_xp_proj_xw**2 + hat_xp_proj_yw**2 + hat_xp_proj_zw**2 )
                hat_xp_proj_xw = hat_xp_proj_xw/Norm
                hat_xp_proj_yw = hat_xp_proj_yw/Norm
                hat_xp_proj_zw = hat_xp_proj_zw/Norm
                
                if gamma != 0:
                    thetas = np.linspace(0, gamma, Ntheta)
                        
                    s_xw = ldisp*(np.cos(thetas)*hat_xp_xw + np.sin(thetas)*hat_zp_xw)
                    s_yw = ldisp*(np.cos(thetas)*hat_xp_yw + np.sin(thetas)*hat_zp_yw)
                    s_zw = ldisp*(np.cos(thetas)*hat_xp_zw + np.sin(thetas)*hat_zp_zw)
                    ax.plot([x_pw, x_pw + ldisp*hat_xp_proj_xw], [y_pw, y_pw + ldisp*hat_xp_proj_yw], [z_pw, z_pw + hat_xp_proj_zw], c = 'mediumpurple', alpha=0.5)
                    ax.plot(x_pw + s_xw, y_pw + s_yw, z_pw + s_zw, c='mediumpurple', alpha= 0.5, label = '$\\gamma$ : roof angle ; $\\sin\\gamma = \\frac{\\hat{\mathbf{x}}_p \\cdot \\hat{\\mathbf{z}}_w}{\\cos\\beta}$')
                    ax.text(x_pw + s_xw[int(Ntheta/2)], y_pw + s_yw[int(Ntheta/2)], z_pw + s_zw[int(Ntheta/2)], "$\\gamma$", c = 'mediumpurple', fontsize = 20)
                
                                
                ax.set_xlabel("$x_w'$", fontsize = 20)
                ax.set_ylabel("$y_w'$", fontsize = 20)
                ax.set_zlabel("$z_w'$", fontsize = 20)
                
                ax.legend()
                
                ax.set_xlim(x_pw+-ldisp, x_pw+ldisp)
                ax.set_ylim(y_pw+-ldisp, y_pw+ldisp)
                ax.set_zlim(z_pw+-ldisp, z_pw+ldisp) 
                ax.set_box_aspect([1,1,1])
                
                                
# ============================================================================


# ============================================================================
# mechanical scan patterns


# ============================================================================
