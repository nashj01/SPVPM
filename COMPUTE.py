# FUNCTION - COMPUTE
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

np.seterr('raise')

# %% COMPUTE THE IJ SOURCE GEOMETRIC INTEGRALS
#
#

def COMPUTE_IJ_SPM(XC,ZC,XB,ZB,phi,S):
    
    # Number of panels
    numPan = len(XC)                                                                # Number of panels/control points
    
    # Initialize arrays
    I = np.zeros([numPan,numPan])                                                   # Initialize I integral matrix
    J = np.zeros([numPan,numPan])                                                   # Initialize J integral matrix
    
    # Compute integral
    for i in range(numPan):                                                         # Loop over i panels
        for j in range(numPan):                                                     # Loop over j panels
            if (j != i):                                                            # If the i and j panels are not the same
                # Compute intermediate values
                A  = -(XC[i]-XB[j])*np.cos(phi[j])-(ZC[i]-ZB[j])*np.sin(phi[j])     # A term
                B  = (XC[i]-XB[j])**2 + (ZC[i]-ZB[j])**2                            # B term
                Cn = np.sin(phi[i]-phi[j])                                          # C term (normal)
                Dn = -(XC[i]-XB[j])*np.sin(phi[i])+(ZC[i]-ZB[j])*np.cos(phi[i])     # D term (normal)
                Ct = -np.cos(phi[i]-phi[j])                                         # C term (tangential)
                Dt = (XC[i]-XB[j])*np.cos(phi[i])+(ZC[i]-ZB[j])*np.sin(phi[i])      # D term (tangential)
                E  = np.sqrt(B-A**2)                                                # E term
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):       # If E term is 0 or complex or a NAN or an INF
                    I[i,j] = 0                                                      # Set I value equal to zero
                    J[i,j] = 0                                                      # Set J value equal to zero
                else:
                    # Compute I (needed for normal velocity), Ref [1]
                    term1  = 0.5*Cn*np.log((S[j]**2 + 2*A*S[j] + B)/B)              # First term in I equation
                    term2  = ((Dn-A*Cn)/E)*(math.atan2((S[j]+A),E)-math.atan2(A,E)) # Second term in I equation
                    I[i,j] = term1 + term2                                          # Compute I integral
                    
                    # Compute J (needed for tangential velocity), Ref [2]
                    term1  = 0.5*Ct*np.log((S[j]**2 + 2*A*S[j] + B)/B)              # First term in I equation
                    term2  = ((Dt-A*Ct)/E)*(math.atan2((S[j]+A),E)-math.atan2(A,E)) # Second term in I equation
                    J[i,j] = term1 + term2                                          # Compute J integral
                
            # Zero out any problem values
            if (np.iscomplex(I[i,j]) or np.isnan(I[i,j]) or np.isinf(I[i,j])):      # If I term is complex or a NAN or an INF
                I[i,j] = 0                                                          # Set I value equal to zero
            if (np.iscomplex(J[i,j]) or np.isnan(J[i,j]) or np.isinf(J[i,j])):      # If J term is complex or a NAN or an INF
                J[i,j] = 0                                                          # Set J value equal to zero
    
    return I, J                                                                     # Return both I and J matrices

# %% COMPUTE THE KL VORTEX GEOMETRIC INTEGRALS
#
#

def COMPUTE_KL_VPM(XC,ZC,XB,ZB,phi,S):
    
    # Number of panels
    numPan = len(XC)                                                                # Number of panels
    
    # Initialize arrays
    K = np.zeros([numPan,numPan])                                                   # Initialize K integral matrix
    L = np.zeros([numPan,numPan])                                                   # Initialize L integral matrix
    
    # Compute integral
    for i in range(numPan):                                                         # Loop over i panels
        for j in range(numPan):                                                     # Loop over j panels
            if (j != i):                                                            # If panel j is not the same as panel i
                # Compute intermediate values
                A  = -(XC[i]-XB[j])*np.cos(phi[j])-(ZC[i]-ZB[j])*np.sin(phi[j])     # A term
                B  = (XC[i]-XB[j])**2 + (ZC[i]-ZB[j])**2                            # B term
                Cn = -np.cos(phi[i]-phi[j])                                         # C term (normal)
                Dn = (XC[i]-XB[j])*np.cos(phi[i])+(ZC[i]-ZB[j])*np.sin(phi[i])      # D term (normal)
                Ct = np.sin(phi[j]-phi[i])                                          # C term (tangential)
                Dt = (XC[i]-XB[j])*np.sin(phi[i])-(ZC[i]-ZB[j])*np.cos(phi[i])      # D term (tangential)
                E  = np.sqrt(B-A**2)                                                # E term
                if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):       # If E term is 0 or complex or a NAN or an INF
                    K[i,j] = 0                                                      # Set K value equal to zero
                    L[i,j] = 0                                                      # Set L value equal to zero
                else:
                    # Compute K
                    term1  = 0.5*Cn*np.log((S[j]**2 + 2*A*S[j] + B)/B)              # First term in K equation
                    term2  = ((Dn-A*Cn)/E)*(math.atan2((S[j]+A),E)-math.atan2(A,E)) # Second term in K equation
                    K[i,j] = term1 + term2                                          # Compute K integral
                    
                    # Compute L
                    term1  = 0.5*Ct*np.log((S[j]**2 + 2*A*S[j] + B)/B)              # First term in L equation
                    term2  = ((Dt-A*Ct)/E)*(math.atan2((S[j]+A),E)-math.atan2(A,E)) # Second term in L equation
                    L[i,j] = term1 + term2                                          # Compute L integral
            
            # Zero out any problem values
            if (np.iscomplex(K[i,j]) or np.isnan(K[i,j]) or np.isinf(K[i,j])):      # If K term is complex or a NAN or an INF
                K[i,j] = 0                                                          # Set K value equal to zero
            if (np.iscomplex(L[i,j]) or np.isnan(L[i,j]) or np.isinf(L[i,j])):      # If L term is complex or a NAN or an INF
                L[i,j] = 0                                                          # Set L value equal to zero
    
    return K, L                             