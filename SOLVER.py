# FUNCTION - SOLVER
# Written by: Jordan Nash
# Started   : 24/01/2024
# Updated   :
# 
#
# PURPOSE
# - Executes the Source Panel Vortex Panel method
#
# %% IMPORT LIBRARIES AND FUNCTIONS
#
#

import numpy as np

from COMPUTE import COMPUTE_IJ_SPM
from COMPUTE import COMPUTE_KL_VPM


# %% IMPORT LIBRARIES AND FUNCTIONS
#
#

def A_MAT(flow_data, multi):
    
    I, J = COMPUTE_IJ_SPM(XC,ZC,XB,ZB,phi,S)  
    K, L = COMPUTE_KL_VPM(XC,ZC,XB,ZB,phi,S) 
    
    
    return