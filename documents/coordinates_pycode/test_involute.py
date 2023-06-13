# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:31:32 2022

@author: guillaume.p-april
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import  splrep, splev, LinearNDInterpolator
from scipy.spatial import Delaunay

dp = 0.01011
deltav_deltau = 0.010021


def MySampledParabolicSpecimen(dp):
    
    W = 1
    L = 1
    a = -0.8
    
    NL = int(L/dp)
    NW = int(W/dp)
    
    x = np.linspace(-W/2, W/2, NW)
    y = np.linspace(-L/2, L/2, NL)
    
    X , Y = np.meshgrid(x, y)
    Z = a*X**2

    NormalX = 2*a*X
    NormalY = 0
    NormalZ = -1
    
    Norm = np.sqrt( NormalX**2 + NormalY**2 + NormalZ**2 )
    
    NormalX = NormalX/Norm
    NormalY = NormalY/Norm
    NormalZ = NormalZ/Norm
    
    return X, Y, Z, NormalX, NormalY, NormalZ 



def MySampledCylinderSpecimen(dp):
    Dx = 0.5
    Dy = 0.5
    R0 = 0.75
    
    Dphi = np.arcsin(Dx/R0)
    phi0 = np.pi/2 
    
    dl = R0*np.sin(Dphi)
    y = np.linspace(-dl, dl, int(2*dl/dp))
    phi = np.linspace( phi0 - Dphi, phi0 + Dphi, int(2*dl/dp) )
    
    phi_grid, Y = np.meshgrid(phi, y)
    
    X = R0*np.cos(phi_grid)
    Z = R0*np.sin(phi_grid)
    
    # normal to surface
    NormalX = -np.cos(phi_grid)
    NormalY =  np.zeros( np.shape(phi_grid) )
    NormalZ = -np.sin(phi_grid)
    
    return X, Y, Z, NormalX, NormalY, NormalZ 
    

def MyParabolicTrajectory(X, Y, Z, NormalX, NormalY, NormalZ, deltau):
    
    interp_Z, interp_NormalX, interp_NormalY, interp_NormalZ = get_ZAndNormInterp(X, Y, Z, NormalX, NormalY, NormalZ)
    
    a = 1.0
    
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()
    
    NormalX = NormalX.flatten()
    NormalY = NormalY.flatten()
    NormalZ = NormalZ.flatten()
    
    x = np.linspace( np.min(X), np.max(X), 100 )
    y = a*x**2
    
    y_log = (y >= np.min(Y) )*(y <= np.max(Y) ) == 1
    x = x[y_log]
    y = y[y_log]
    z = np.zeros( (len(x), ) )

    normalx = np.zeros( (len(x), ) )
    normaly = np.zeros( (len(x), ) )
    normalz = np.zeros( (len(x), ) )
    
    for p in range(len(x)):
        
        
        z[p] = interp_Z( np.asarray( [ x[p], y[p] ] ) )
        
        normalx[p] = interp_NormalX( np.asarray( [ x[p], y[p] ] ) )    
        normaly[p] = interp_NormalY( np.asarray( [ x[p], y[p] ] ) )
        normalz[p] = interp_NormalZ( np.asarray( [ x[p], y[p] ] ) )
    
        
    u = get_u(x, y, z)
    
    u_uniform = np.arange(np.min(u), np.max(u), deltav_deltau )    
    
    tck_u_x = splrep( u, x, k=3, s=0 )
    tck_u_y = splrep( u, y, k=3, s=0 )
    tck_u_z = splrep( u, z, k=3, s=0 )
            
    tck_u_normalx = splrep( u, normalx, k=3, s=0 )
    tck_u_normaly = splrep( u, normaly, k=3, s=0 )
    tck_u_normalz = splrep( u, normalz, k=3, s=0 )
    
    x = splev(u_uniform, tck_u_x, der = 0) 
    y = splev(u_uniform, tck_u_y, der = 0) 
    z = splev(u_uniform, tck_u_z, der = 0)
    
    normalx = splev(u_uniform, tck_u_normalx, der = 0) 
    normaly = splev(u_uniform, tck_u_normaly, der = 0) 
    normalz = splev(u_uniform, tck_u_normalz, der = 0)
    

    return x, y, z, normalx, normaly, normalz


def MyLineTrajectory(X, Y, Z, NormalX, NormalY, NormalZ, deltau):
    
    interp_Z, interp_NormalX, interp_NormalY, interp_NormalZ = get_ZAndNormInterp(X, Y, Z, NormalX, NormalY, NormalZ)
    
    Theta = 0.65
    R0 = 0.75
    
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()
    
    NormalX = NormalX.flatten()
    NormalY = NormalY.flatten()
    NormalZ = NormalZ.flatten()
    
    if np.abs( Theta-np.pi/2 ) > 1e-3: 
        r = np.linspace( np.min(Y)/np.cos(Theta), np.max(Y)/np.cos(Theta), 100 )
    else:
        r = np.linspace( np.min(X)/np.sin(Theta), np.max(X)/np.sin(Theta), 100 )
    
        
    u = np.linspace(r[0], r[-1], 100)
    x = R0*np.cos(1/R0*np.sin(Theta)*u + np.pi/2) 
    y = u*np.cos(Theta) 
    
    y_log = (y >= np.min(Y) )*(y <= np.max(Y) ) == 1
    x = x[y_log]
    y = y[y_log]
    x_log = (x >= np.min(X) )*(x <= np.max(X) ) == 1
    x = x[x_log]
    y = y[x_log]
    
    z = np.zeros( (len(x), ) )

    normalx = np.zeros( (len(x), ) )
    normaly = np.zeros( (len(x), ) )
    normalz = np.zeros( (len(x), ) )
    
    for p in range(len(x)):
        
        
        z[p] = interp_Z( np.asarray( [ x[p], y[p] ] ) )
        
        normalx[p] = interp_NormalX( np.asarray( [ x[p], y[p] ] ) )    
        normaly[p] = interp_NormalY( np.asarray( [ x[p], y[p] ] ) )
        normalz[p] = interp_NormalZ( np.asarray( [ x[p], y[p] ] ) )
    
        
    u = get_u(x, y, z)
    
    u_uniform = np.arange(np.min(u), np.max(u), deltav_deltau )    
    
    tck_u_x = splrep( u, x, k=3, s=0 )
    tck_u_y = splrep( u, y, k=3, s=0 )
    tck_u_z = splrep( u, z, k=3, s=0 )
            
    tck_u_normalx = splrep( u, normalx, k=3, s=0 )
    tck_u_normaly = splrep( u, normaly, k=3, s=0 )
    tck_u_normalz = splrep( u, normalz, k=3, s=0 )
    
    x = splev(u_uniform, tck_u_x, der = 0) 
    y = splev(u_uniform, tck_u_y, der = 0) 
    z = splev(u_uniform, tck_u_z, der = 0)
    
    normalx = splev(u_uniform, tck_u_normalx, der = 0) 
    normaly = splev(u_uniform, tck_u_normaly, der = 0) 
    normalz = splev(u_uniform, tck_u_normalz, der = 0)
    

    return x, y, z, normalx, normaly, normalz



def MyLineTrajectory_Analytical(u, v, X0, Y0, Z0):
     
    R0 = 0.75
    Theta = 0.65
    
    x = R0*np.cos(-np.cos(Theta)/R0*v + np.sin(Theta)/R0*u + np.arctan2( Z0, X0 ))
    y = v*np.sin(Theta) + u*np.cos(Theta) + Y0
    z = R0*np.sin(-np.cos(Theta)/R0*v + np.sin(Theta)/R0*u + np.arctan2( Z0, X0 ))

    return x, y, z



def get_ZAndNormInterp(X, Y, Z, NormalX, NormalY, NormalZ):
    
    X_flat = X.flatten()
    Y_flat = Y.flatten()
    Z_flat = Z.flatten()
            
    points = np.zeros( (len(X_flat), 2))
    points[:,0] = X_flat
    points[:,1] = Y_flat
    
    tri = Delaunay( points )  # Compute the triangulation

    # interpolaters
    interp_Z = LinearNDInterpolator(tri, Z_flat)
    interp_NormalX = LinearNDInterpolator( tri, NormalX.flatten() )
    interp_NormalY = LinearNDInterpolator( tri, NormalY.flatten() )
    interp_NormalZ = LinearNDInterpolator( tri, NormalZ.flatten() )
    
    
    return interp_Z, interp_NormalX, interp_NormalY, interp_NormalZ
    
    


def get_u(x, y, z):
    
    # local line u
    u = np.zeros( (len(x), ) )
    for p in np.arange(1, len(u)):
        u[p] = u[p-1] + np.sqrt(  (x[p]-x[p-1])**2 
                                + (y[p]-y[p-1])**2 
                                + (z[p]-z[p-1])**2 )
        
    return u
    


def get_NextInvolute(X_flat, Y_flat, Z_flat, interp_Z, interp_NormalX, interp_NormalY, interp_NormalZ, x, y, z, normalx, normaly, normalz, deltav, u):
    
    Xmin = np.min(X_flat)
    Xmax = np.max(X_flat)
    Ymin = np.min(Y_flat)
    Ymax = np.max(Y_flat)
    Zmin = np.min(Z_flat)
    Zmax = np.max(Z_flat)
    
    
    tck_u_x = splrep( u, x, k=3, s=0 )
    tck_u_y = splrep( u, y, k=3, s=0 )
    tck_u_z = splrep( u, z, k=3, s=0 )

    du_x = splev(u, tck_u_x, der = 1) 
    du_y = splev(u, tck_u_y, der = 1) 
    du_z = splev(u, tck_u_z, der = 1)
    
    Norm = np.sqrt( du_x**2 + du_y**2 + du_z**2 )    
    du_x = du_x/Norm
    du_y = du_y/Norm
    du_z = du_z/Norm
    
    dv_x =  normaly*du_z - normalz*du_y
    dv_y = -normalx*du_z + normalz*du_x
    dv_z =  normalx*du_y - normaly*du_x
    
    Norm = np.sqrt( dv_x**2 + dv_y**2 + dv_z**2 )
    dv_x = dv_x/Norm
    dv_y = dv_y/Norm
    dv_z = dv_z/Norm
    
    x_next = []
    y_next = []
    z_next = []
    
    normalx_next = []
    normaly_next = []
    normalz_next = []
    
    u_next = []
    
    
    for p in range( len(x) ):
        
        x_test = x[p] + deltav*dv_x[p]
        y_test = y[p] + deltav*dv_y[p]
        z_test = interp_Z( np.asarray( [x_test, y_test] ) )[0]
        
        
        if x_test <= Xmin or x_test >= Xmax:
            bleh = 1
        elif y_test <= Ymin or y_test >= Ymax:
            bleh = 1
        else:
            x_next.append(x_test)
            y_next.append(y_test)
            z_next.append(z_test)
            
            normalx_next.append( interp_NormalX( np.asarray( [x_test, y_test] ) )[0] )
            normaly_next.append( interp_NormalY( np.asarray( [x_test, y_test] ) )[0] )
            normalz_next.append( interp_NormalZ( np.asarray( [x_test, y_test] ) )[0] )
            
            u_next.append( u[p] )
            

            
    x_next = np.asarray(x_next)
    y_next = np.asarray(y_next)    
    z_next = np.asarray(z_next)
    
    normalx_next = np.asarray(normalx_next)
    normaly_next = np.asarray(normaly_next)    
    normalz_next = np.asarray(normalz_next)
    
    u_next = np.asarray(u_next)
    
    return x_next, y_next, z_next, normalx_next, normaly_next, normalz_next, u_next
    


def get_involutes(x_ref, y_ref, z_ref, normalx_ref, normaly_ref, normalz_ref, X, Y, Z, NormalX, NormalY, NormalZ, deltav):
    
    Deltav_over_deltav = 0
    
    interp_Z, interp_NormalX, interp_NormalY, interp_NormalZ = get_ZAndNormInterp(X, Y, Z, NormalX, NormalY, NormalZ)
    
    Xmin = np.min(X)
    Xmax = np.max(X)
    Ymin = np.min(Y)
    Ymax = np.max(Y)
    Zmin = np.min(Z)
    Zmax = np.max(Z)
    
    X_flat = X.flatten()
    Y_flat = Y.flatten()
    Z_flat = Z.flatten()
    
    NormalX_flat = NormalX.flatten()
    NormalY_flat = NormalY.flatten()
    NormalZ_flat = NormalZ.flatten()
    
    x_cur = x_ref
    y_cur = y_ref
    z_cur = z_ref
    
    normalx_cur = normalx_ref
    normaly_cur = normaly_ref
    normalz_cur = normalz_ref
    
    u = get_u(x_ref, y_ref, z_ref)
    u_cur = u
    
    Involutes_FWD_x = [x_ref]
    Involutes_FWD_y = [y_ref]
    Involutes_FWD_z = [z_ref]

    p = 0
    v_FWD = [0]
    u_FWD = [u_cur]
    while len(u_cur) > 10:
        x_next, y_next, z_next, normalx_next, normaly_next, normalz_next, u_next = get_NextInvolute(X_flat, Y_flat, Z_flat, 
                                                                                            interp_Z, interp_NormalX, interp_NormalY, interp_NormalZ, 
                                                                                            x_cur, y_cur, z_cur, 
                                                                                            normalx_cur, normaly_cur, normalz_cur, 
                                                                                            deltav, u_cur )
        
        Involutes_FWD_x.append( x_next )
        Involutes_FWD_y.append( y_next )
        Involutes_FWD_z.append( z_next )
        
        v_FWD.append( v_FWD[-1] + 1*deltav )
        u_FWD.append( u_next )
        print('v = ' + str(v_FWD[-1]))

        # if p == Deltav_over_deltav:
        #     Involutes_FWD_x.append( x_next )
        #     Involutes_FWD_y.append( y_next )
        #     Involutes_FWD_z.append( z_next )
            
        #     v_FWD.append( v_FWD[-1] + Deltav_over_deltav*deltav )
        #     u_FWD.append( u_next )
        #     p = 0
        #     print('v = ' + str(v_FWD[-1]))
        # else:
        #     p += 1
                
        x_cur = x_next
        y_cur = y_next
        z_cur = z_next
        
        normalx_cur = normalx_next
        normaly_cur = normaly_next
        normalz_cur = normalz_next
        
        u_cur = u_next
        
    v_FWD = np.asarray( v_FWD )
    
    
    
    u_cur = u
    x_cur = x_ref
    y_cur = y_ref
    z_cur = z_ref
    
    normalx_cur = normalx_ref
    normaly_cur = normaly_ref
    normalz_cur = normalz_ref
    
    
    Involutes_BWD_x = []
    Involutes_BWD_y = []
    Involutes_BWD_z = []

    p = 0
    v_BWD = []
    u_BWD = []
    while len(u_cur) > 10:
        x_next, y_next, z_next, normalx_next, normaly_next, normalz_next, u_next = get_NextInvolute(X_flat, Y_flat, Z_flat, 
                                                                                            interp_Z, interp_NormalX, interp_NormalY, interp_NormalZ, 
                                                                                            x_cur, y_cur, z_cur, 
                                                                                            normalx_cur, normaly_cur, normalz_cur, 
                                                                                            -deltav, u_cur )

        Involutes_BWD_x.append( x_next )
        Involutes_BWD_y.append( y_next )
        Involutes_BWD_z.append( z_next )
        
        if len(v_BWD) == 0:
            v_BWD.append( 1*-deltav )
        else:
            v_BWD.append( v_BWD[-1] + 1*-deltav )
        u_BWD.append( u_next )
        print('v = ' + str(v_BWD[-1]))
        
        # if p == Deltav_over_deltav:
        #     Involutes_BWD_x.append( x_next )
        #     Involutes_BWD_y.append( y_next )
        #     Involutes_BWD_z.append( z_next )
            
        #     if len(v_BWD) == 0:
        #         v_BWD.append( Deltav_over_deltav*-deltav )
        #     else:
        #         v_BWD.append( v_BWD[-1] + Deltav_over_deltav*-deltav )
        #     u_BWD.append( u_next )
        #     p = 0
        #     print('v = ' + str(v_BWD[-1]))
        # else:
        #     p += 1
                
        x_cur = x_next
        y_cur = y_next
        z_cur = z_next
        
        normalx_cur = normalx_next
        normaly_cur = normaly_next
        normalz_cur = normalz_next
        
        u_cur = u_next
        
    v_BWD = np.asarray( v_BWD )

    
    v_tmp = np.concatenate( (v_BWD, v_FWD) )
    u_BWD.extend( u_FWD )
    Involutes_BWD_x.extend( Involutes_FWD_x )
    Involutes_BWD_y.extend( Involutes_FWD_y )
    Involutes_BWD_z.extend( Involutes_FWD_z )
    
    ind = np.argsort( v_tmp )

    Involutes_x = []
    Involutes_y = []
    Involutes_z = []
    v2 = []
    u = []
    for p in range( len(v_tmp) ):
        Involutes_x.append( Involutes_BWD_x[ind[p]] )
        Involutes_y.append( Involutes_BWD_y[ind[p]] )
        Involutes_z.append( Involutes_BWD_z[ind[p]] )
        
        v2.append( v_tmp[ind[p]] )
        u.append( u_BWD[ind[p]] )
        
        
    PerpInvolutes_x = []
    PerpInvolutes_y = []
    PerpInvolutes_z = []
    v = []

    len_u = 0    
    for p in range(len(u)):
        if len_u < len(u[p]):
            u_ref = u[p]
            len_u = len(u[p])

    for p in range(len(u_ref)):
        u0 = u_ref[p]
        
        PerpInvolutes_x_tmp = []
        PerpInvolutes_y_tmp = []
        PerpInvolutes_z_tmp = []
        v_tmp = []
        
        for q in range(len(u)):
            u_log = ( np.abs(u0 - u[q]) < 10e-11 )
            
            if sum(u_log) > 0:
                PerpInvolutes_x_tmp.append( Involutes_x[q][u_log][0] )
                PerpInvolutes_y_tmp.append( Involutes_y[q][u_log][0] )
                PerpInvolutes_z_tmp.append( Involutes_z[q][u_log][0] )
            
                v_tmp.append( v2[q] )
                
        PerpInvolutes_x.append( PerpInvolutes_x_tmp )
        PerpInvolutes_y.append( PerpInvolutes_y_tmp )
        PerpInvolutes_z.append( PerpInvolutes_z_tmp )
        v.append( np.asarray(v_tmp) )
        
        
    v_alongPerpInvolutes = v
    v_alongInvolutes = v2
    u_alongInvolutes = u
    
    return u_alongInvolutes, v_alongInvolutes, v_alongPerpInvolutes, Involutes_x, Involutes_y, Involutes_z, PerpInvolutes_x, PerpInvolutes_y, PerpInvolutes_z
    


# specimen shape and reference line, sampled on tight grid     

# X, Y, Z, NormalX, NormalY, NormalZ  = MySampledParabolicSpecimen(dp)
X, Y, Z, NormalX, NormalY, NormalZ  = MySampledCylinderSpecimen(dp)

# x_ref, y_ref, z_ref, normalx_ref, normaly_ref, normalz_ref = MyParabolicTrajectory(X, Y, Z, NormalX, NormalY, NormalZ, deltau = deltav_deltau)
x_ref, y_ref, z_ref, normalx_ref, normaly_ref, normalz_ref = MyLineTrajectory(X, Y, Z, NormalX, NormalY, NormalZ, deltau = deltav_deltau)


#  
u, v_alongInvolutes, v, Involutes_x, Involutes_y, Involutes_z, PerpInvolutes_x, PerpInvolutes_y, PerpInvolutes_z = get_involutes(x_ref, y_ref, z_ref, normalx_ref, normaly_ref, normalz_ref, X, Y, Z, NormalX, NormalY, NormalZ, deltav = deltav_deltau)


u_anal = np.linspace(0, 0.5, 25)
v_anal = np.zeros( (len(u_anal), ) ) 
x_anal0, y_anal0, z_anal0 = MyLineTrajectory_Analytical(u_anal, v_anal, x_ref[0], y_ref[0], z_ref[0])

u_anal = np.full( (25, ) , u_anal[-1] )
v_anal =  np.linspace(0, 0.5, len(u_anal))
x_anal1, y_anal1, z_anal1 = MyLineTrajectory_Analytical(u_anal, v_anal, x_ref[0], y_ref[0], z_ref[0])



fig = plt.figure(figsize = (13, 12))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, alpha = 0.3, color = 'gray')
ax.plot(x_ref, y_ref, z_ref, linewidth = 4, c = 'orange')
ax.plot(x_anal0, y_anal0, z_anal0, linewidth = 4, c = 'blue')
ax.plot(x_anal1, y_anal1, z_anal1, linewidth = 4, c = 'blue')
ldisp = 0.05*np.max( [ np.max(X)-np.min(X), np.max(Y)-np.min(Y), np.max(Z)-np.min(Z) ] )
# for p in range(len(normalx_ref)):
#     plt.plot([x_ref[p], x_ref[p]+ldisp*normalx_ref[p] ], 
#          [y_ref[p], y_ref[p]+ldisp*normaly_ref[p] ], 
#          [z_ref[p], z_ref[p]+ldisp*normalz_ref[p] ], c = 'gray')

for p in range(len(u)):
    
    if len(u[p]) > 0:
        ind0 = 0
        ind1 = len(u[p])
        while ind0 != len(u[p]):
            for q in range(ind0, len(u[p])-1):
                if u[p][q+1]-u[p][q] > 1.1*deltav_deltau:
                    ind1 = q + 1
                    break
            
            if ind0 != ind1:
                plt.plot(Involutes_x[p][ind0:ind1], Involutes_y[p][ind0:ind1], Involutes_z[p][ind0:ind1], c = 'red', alpha = 0.5)
            
            ind0 = ind1 
            ind1 = len(u[p])
    
for p in range(len(v)):
    
    if len(v[p]) > 0:
        ind0 = 0
        ind1 = len(v[p])
        while ind0 != len(v[p]):
            for q in range(ind0, len(v[p])-1):
                if v[p][q+1]-v[p][q] > 1.1*deltav_deltau:
                    ind1 = q + 1
                    break
            
            if ind0 != ind1:
                plt.plot(PerpInvolutes_x[p][ind0:ind1], PerpInvolutes_y[p][ind0:ind1], PerpInvolutes_z[p][ind0:ind1], c = 'green', alpha = 0.5)
            
            ind0 = ind1 
            ind1 = len(v[p])



    

ax.set_xlabel("$x$", fontsize = 20)
ax.set_ylabel("$y$", fontsize = 20)
ax.set_zlabel("$z$", fontsize = 20)


ax.set_xlim(np.min(X), np.max(X))
ax.set_ylim(np.min(Y), np.max(Y))
ax.set_zlim(np.min(Z), np.max(Z))

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

ind = np.where( np.abs(v_alongInvolutes) < 1e-14 )[0][0]
ax.plot( u[ind], v_alongInvolutes[ind]*np.ones( (len(u[ind], )) ),  linewidth = 4, c = 'orange', alpha = 0.5)
for p in range(len(u)):
    if len(u[p]) > 0:
        ax.scatter(u[p], v_alongInvolutes[p]*np.ones( (len(u[p]), )), c = 'cornflowerblue', s = 8)
ax.set_aspect('equal', adjustable='box')
fig.tight_layout()

print('test on diff value of u on ref line, target: deltau = ' + str(deltav_deltau) )
for p in range(len(u[ind])-1):
    print('   ' + str( np.sqrt( ( Involutes_x[ind][p+1]-Involutes_x[ind][p] )**2 +
                       ( Involutes_y[ind][p+1]-Involutes_y[ind][p] )**2 +
                       ( Involutes_z[ind][p+1]-Involutes_z[ind][p] )**2 ) ) )

print('test on diff value of v near ref line, target: deltav = ' + str(deltav_deltau) )
for p in range(len(v)):
    
    ind = np.where( np.abs(v[p]) < 1e-14 )[0]

    if len(ind) > 0:
        ind = ind[0]
        
        if ind == len(v[p])-1:
            print('   ' + str( np.sqrt( ( PerpInvolutes_x[p][ind-1]-PerpInvolutes_x[p][ind] )**2 +
               ( PerpInvolutes_y[p][ind-1]-PerpInvolutes_y[p][ind] )**2 +
               ( PerpInvolutes_z[p][ind-1]-PerpInvolutes_z[p][ind] )**2 ) ) )
        else:
            print('   ' + str( np.sqrt( ( PerpInvolutes_x[p][ind+1]-PerpInvolutes_x[p][ind] )**2 +
                           ( PerpInvolutes_y[p][ind+1]-PerpInvolutes_y[p][ind] )**2 +
                           ( PerpInvolutes_z[p][ind+1]-PerpInvolutes_z[p][ind] )**2 ) ) )
        
        
    