# FUNCTION - SPVP
# Written by: Jordan Nash
# Started   : 24/01/2024
# Updated   :
# 
#
# PURPOSE
# - Frontend for panel method code
#
# %% IMPORT LIBRARIES AND FUNCTIONS
#
#

from PREPROCESS import XML
from PREPROCESS import GEOMETRY
from PREPROCESS import PLTGEOMETRY
from PREPROCESS import PRECON
# %% PREPROCESSOR
#
#

# Import scripted geometry data
file = XML("NACA2412.xml")

# Compute geometry
geometry = GEOMETRY(file)

# Plot geometry
PLTGEOMETRY(geometry)

# Console plot
PRECON(file, geometry)
# %% SOLVER
#
#





# %% POSTPROCESSOR
#
#

