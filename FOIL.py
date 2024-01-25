# FUNCTION - FOIL
# Written by: Jordan Nash
# Started   : 11/01/2024 - Transferred from MATLAB to Python
# Updated   :
#
# PURPOSE
# - Compute the integral expression for constant strength source panels
#
# %% IMPORT LIBRARIES AND FUNCTIONS
#
#
import numpy as np
import math as math
# %% COMPUTE FOIL GEOMETRY
#
#

def FOIL(nx, c, m, p, t, alpha, Xgap, Zgap, offset, deltaB, deltaC):
   
    x = (1-np.cos(np.linspace(0, np.pi, nx)))/2
   
    # CAMBER DEFINITION
    Zc_1 = ((m*c) / (p**2)) * (2 * p *x[0:int(np.round(nx * p))] - x[0:int(np.round(nx * p))]**2)
    Zc_2 = ((m*c) / ((1-p)**2)) * (1-2 * p + 2 * p *x[int(np.round(nx * p)):nx] - x[int(np.round(nx * p)):nx]**2)
    Zcamb = np.concatenate((Zc_1, Zc_2))
    Xcamb = x*c
   
    # GRADIENT DEIFNITION
    dZc_1 = ((2*m) / (p**2)) * (p - x[0:int(np.round(nx * p))])
    dZc_2 = ((2*m) / ((1-p)**2)) * (p - x[int(np.round(nx * p)):nx])
    dZc = np.concatenate((dZc_1, dZc_2))
   
    # THICKNESS DEFINITION
    Zt = ((c*t)/(0.2))*(0.2969*x**0.5 - 0.126*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)
   
    theta = np.arctan(dZc)
   
    # UPPER AND LOWER SURFACE DEFINITION
    Xu = x - Zt*np.sin(theta)
    Xl = x + Zt*np.sin(theta)
    XB = c*np.concatenate((np.flip(Xl[1:len(Xl)]), Xu))
   
    Zu = Zcamb + Zt*np.cos(theta)
    Zl = Zcamb - Zt*np.cos(theta)
    ZB = np.concatenate((np.flip(Zl[1:len(Zl)]), Zu)) 
    
    # PREALLOCATION
    XC = np.zeros(2*nx-2)
    ZC = np.zeros(2*nx-2)
    beta = np.zeros(2*nx-2)
    phi = np.zeros(2*nx-2)
    S = np.zeros(2*nx-2)
    dx = np.zeros(2*nx-2)
    dz = np.zeros(2*nx-2)
    
    # COLLOCATION POINTS AND PANEL ORIENTATIONS
    for i in range(2*nx-2):
        XC[i] = (XB[i] + XB[i+1])/2
        ZC[i] = (ZB[i] + ZB[i+1])/2
        dx[i]      = XB[i+1]-XB[i]                                                     # Change in X between boundary points
        dz[i]      = ZB[i+1]-ZB[i]                                                     # Change in Y between boundary points
        S[i]    = (dx[i]**2 + dz[i]**2)**0.5                                              # Length of the panel
        phi[i]  = math.atan2(dz[i],dx[i])                                                 # Angle of panel (positive X-axis to inside face)
        if (phi[i] < 0):                                                            # Make all panel angles positive [rad]
            phi[i] = phi[i] + 2*np.pi
    
    beta = ((phi + (np.pi/2)) - alpha)                                          # Angle of panel normal [rad]
    beta[beta > 2*np.pi] = beta[beta > 2*np.pi] - 2*np.pi                           # Make all panel angles between 0 and 2pi [rad]
    
    # Compute the boundary layer panel boundary points and panel collocation points
    DXC = np.diff(XC)
    DZC = np.diff(ZC)
    
    DXB = np.diff(XB)
    DZB = np.diff(ZB)
    
    normC = np.sqrt(DXC**2 + DZC**2)
    normB = np.sqrt(DXB**2 + DZB**2)
    
    NXC = np.concatenate(([DZC[0]], DZC / normC))
    NZC = np.concatenate(([-DXC[0]], DXC / normC)) 
    NXB = np.concatenate(([DZB[0]], DZB / normB))
    NZB = np.concatenate(([-DXB[0]], DXB / normB)) 
    
    XC = XC - deltaC*NXC + offset + Xgap
    ZC = ZC + deltaC*NZC + Zgap
    XB = XB - deltaB*NXB + offset + Xgap
    ZB = ZB + deltaB*NZB + Zgap
    
    # Compute the chord 
    Xchord = np.linspace(0, c, 2*nx-2) + offset + Xgap
    Zchord = np.zeros(2*nx-2) + Zgap
    
    # Compute the camber
    Xcamb = Xcamb + offset + Xgap
    Zcamb = Zcamb + Zgap

    # Compute the offset for next foil
    offset = XB[0]
    
    return XB, ZB, XC, ZC, Xchord, Zchord, Xcamb, Zcamb, offset, alpha, beta, phi, S


