# FUNCTION - SPVP
# Written by: Jordan Nash
# Started   : 24/01/2024
# Updated   : 27/01/2024
# 
#
# PURPOSE
# - Frontend for panel method code
#
# %% IMPORT LIBRARIES
#
#

import time
import os

# %% IMPORT FUNCTIONS
#
#

from PREPROCESSOR import XML
from PREPROCESSOR import GEOMETRY
from PREPROCESSOR import PLTGEOMETRY
from PREPROCESSOR import PRECONSOLE

#from SOLVER import A_MAT

# %% PREPROCESSOR
# - Imports a pre-configured xml file that defines a single/multi-element airfoil
# - Generates the individual and system discretised geometry
# - Plots the geometry
# - Prints results

os.system('cls' if os.name == 'nt' else 'clear')

status = True

while status == True:
    # Start the timer
    start = time.time()
    
    # Import scripted geometry data
    geometry_data, flow_data = XML("NACA2412_5_ELEMENT.xml")
    
    # Compute geometry
    individual, multi = GEOMETRY(geometry_data)
    
    # Plot geometry
    PLTGEOMETRY(individual, multi)
    
    # End the timer
    end = time.time()
    
    # Calculate the duration
    duration = end-start
    
    # Console plot
    status = PRECONSOLE(flow_data, multi, duration)
         
# %% SOLVER
#
#
    # Start solver
    #A_MAT(flow_data, multi)
    




# %% POSTPROCESSOR
#
#

