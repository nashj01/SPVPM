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

# %% COMPUTE FOIL BOUNDARY POINTS
#
#

def BOUNDPTS(alpha, nx, c, m, p, t, Xgap, Zgap, offset):
   
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
    Zu = Zcamb + Zt*np.cos(theta)
    Zl = Zcamb - Zt*np.cos(theta)
    
    # Evaluation of the boundary points
    XB = c*np.concatenate((np.flip(Xl[1:len(Xl)]), Xu))
    ZB = np.concatenate((np.flip(Zl[1:len(Zl)]), Zu)) 
    XB = ((XB)*np.cos(-alpha) - (ZB)*np.sin(-alpha))
    ZB = ((XB)*np.sin(-alpha) + (ZB)*np.cos(-alpha))
    XB = XB + offset + Xgap
    ZB = ZB + Zgap
    
    # Compute the chord 
    Xchord = np.linspace(0, c, 2*nx-1) 
    Zchord = np.zeros(2*nx-1)
    Xchord = ((Xchord)*np.cos(-alpha) - (Zchord)*np.sin(-alpha))
    Zchord = ((Xchord)*np.sin(-alpha) + (Zchord)*np.cos(-alpha))
    Xchord = Xchord + offset + Xgap
    Zchord = Zchord + Zgap

    
    # Compute the camber
    
    Xcamb = ((Xcamb)*np.cos(-alpha) - (Zcamb)*np.sin(-alpha))
    Zcamb = ((Xcamb)*np.sin(-alpha) + (Zcamb)*np.cos(-alpha))
    Xcamb= Xcamb + offset + Xgap
    Zcamb = Zcamb + Zgap

    # Compute the offset for next foil
    offset = XB[0]
    
    return XB, ZB, Xchord, Zchord, Xcamb, Zcamb, offset

# %% COMPUTE FOIL CONTROL POINTS
#
#

def COLLOPTS(multi):
    
    # Preallocation
    XB = multi["XB"]
    ZB = multi["ZB"]
    XC = np.zeros(len(multi["XB"])-1)
    ZC = np.zeros(len(multi["XB"])-1)
    beta = np.zeros(len(multi["XB"])-1)
    phi = np.zeros(len(multi["XB"])-1)
    S = np.zeros(len(multi["XB"])-1)
    dx = np.zeros(len(multi["XB"])-1)
    dz = np.zeros(len(multi["XB"])-1)
    
    
    # Collocation points and panel orientations
    for j in range(len(XB)-1):
        XC[j] = (XB[j] + XB[j+1])/2
        ZC[j] = (ZB[j] + ZB[j+1])/2
        dx[j] = XB[j+1]-XB[j]                                                     
        dz[j] = ZB[j+1]-ZB[j]                                                    
        S[j] = (dx[j]**2 + dz[j]**2)**0.5                                              
        phi[j] = math.atan2(dz[j],dx[j])                                              
        if (phi[j] < 0):                                                        
            phi[j] = phi[j] + 2*np.pi
    
    
        beta = ((phi + (np.pi/2)) - 2*multi["alpha"])                                          
        beta[beta > 2*np.pi] = beta[beta > 2*np.pi] - 2*np.pi  
    
    return XC, ZC, S, beta, phi















































































    
#     # Compute the boundary layer panel boundary points and panel collocation points
#     DXC = np.diff(XC)
#     DZC = np.diff(ZC)
    
#     DXB = np.diff(XB)
#     DZB = np.diff(ZB)
    
#     normC = np.sqrt(DXC**2 + DZC**2)
#     normB = np.sqrt(DXB**2 + DZB**2)
    
#     NXC = np.concatenate(([DZC[0]], DZC / normC))
#     NZC = np.concatenate(([-DXC[0]], DXC / normC)) 
#     NXB = np.concatenate(([DZB[0]], DZB / normB))
#     NZB = np.concatenate(([-DXB[0]], DXB / normB)) 
    
#     XC = XC - deltaC*NXC + offset + Xgap
#     ZC = ZC + deltaC*NZC + Zgap
#     XB = XB - deltaB*NXB + offset + Xgap
#     ZB = ZB + deltaB*NZB + Zgap
    
#     # Compute the chord 
#     Xchord = np.linspace(0, c, 2*nx-2) + offset + Xgap
#     Zchord = np.zeros(2*nx-2) + Zgap
    
#     # Compute the camber
#     Xcamb = Xcamb + offset + Xgap
#     Zcamb = Zcamb + Zgap

#     # Compute the offset for next foil
#     offset = XB[0]
    
#     return XB, ZB, XC, ZC, Xchord, Zchord, Xcamb, Zcamb, offset, alpha, beta, phi, S