# FUNCTION - COMPUTE X AND Y GEOMETRIC INTEGRALS FOR SOURCE PANEL METHOD
# Written by: JoshTheEngineer
# YouTube   : www.youtube.com/joshtheengineer
# Website   : www.joshtheengineer.com
# Started   : 01/23/19
# Updated   : 01/23/19 - Started code in MATLAB
#                      - Works as expected
#           : 02/03/19 - Transferred to Python
#                      - Works as expected
#           : 04/28/20 - Fixed E value error handling
#
# PURPOSE
# - Compute the geometric integral at point P due to source panels
# - Source panel strengths are constant, but can change from panel to panel
# - Geometric integral for X-direction: Mx(pj)
# - Geometric integral for Y-direction: My(pj)
#
# REFERENCE
# - [1]: Streamline Geometric Integral SPM, Mx(pj) and My(pj)
#           Link: https://www.youtube.com/watch?v=BnPZjGCatcg
#
# INPUTS
# - XP  : X-coordinate of computation point, P
# - YP  : Y-coordinate of computation point, P
# - XB  : X-coordinate of boundary points
# - YB  : Y-coordinate of boundary points
# - phi : Angle between positive X-axis and interior of panel
# - S   : Length of panel
# 
# OUTPUTS
# - Mx  : Value of X-direction geometric integral (Ref [1])
# - My  : Value of Y-direction geometric integral (Ref [1])

import numpy as np
import math as math

np.seterr('raise')

def STREAMLINE_SPM(XP,ZP,XB,ZB,phi,S):
    
    # Number of panels
    numPan = len(XB)-1                                                          # Number of panels
    
    # Initialize arrays
    Mx = np.zeros(numPan)                                                       # Initialize Ix integral array
    Mz = np.zeros(numPan)                                                       # Initialize Iy integral array
    
    # Compute integral
    for j in range(numPan):                                                     # Loop over all panels
        # Compute intermediate values
        A = -(XP-XB[j])*np.cos(phi[j]) - (ZP-ZB[j])*np.sin(phi[j])              # A term
        B  = (XP-XB[j])**2 + (ZP-ZB[j])**2;                                     # B term
        Cx = -np.cos(phi[j]);                                                   # C term (X-direction)
        Dx = XP - XB[j];                                                        # D term (X-direction)
        Cz = -np.sin(phi[j]);                                                   # C term (Y-direction)
        Dz = ZP - ZB[j];                                                        # D term (Y-direction)
        E  = math.sqrt(B-A**2);                                                 # E term
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):           # If E term is 0 or complex or a NAN or an INF
            Mx[j] = 0                                                           # Set Mx value equal to zero
            Mz[j] = 0                                                           # Set My value equal to zero
        else:
            # Compute Mx, Ref [1]
            term1 = 0.5*Cx*np.log((S[j]**2 + 2*A*S[j]+B)/B);                    # First term in Mx equation
            term2 = ((Dx-A*Cx)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));   # Second term in Mx equation
            Mx[j] = term1 + term2;                                              # Compute Mx integral
            
            # Compute My, Ref [1]
            term1 = 0.5*Cz*np.log((S[j]**2 + 2*A*S[j]+B)/B);                    # First term in My equation
            term2 = ((Dz-A*Cz)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));   # Second term in My equation
            Mz[j] = term1 + term2;                                              # Compute My integral

        # Zero out any problem values
        if (np.iscomplex(Mx[j]) or np.isnan(Mx[j]) or np.isinf(Mx[j])):         # If Mx term is complex or a NAN or an INF
            Mx[j] = 0                                                           # Set Mx value equal to zero
        if (np.iscomplex(Mz[j]) or np.isnan(Mz[j]) or np.isinf(Mz[j])):         # If My term is complex or a NAN or an INF
            Mz[j] = 0                                                           # Set My value equal to zero
    
    return Mx, Mz                                                              # Return both Mx and My matrices