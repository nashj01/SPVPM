# Multi-element Source Panel/Vortex Panel Method

![image](https://github.com/nashj01/SPVPM/assets/138644901/819a4e28-bf0c-4090-a833-34971db78092)

This repository contains a simplified linear alternative to CFD for analysing the aerodynamic properties of 2D and soon to be pseudo 3D aerodynamic bodies. Thus far it contains the following subroutines:

 - _PREPROCESSOR_
 - _SOLVER_
 - _POSTPROCESSOR_
 - _SPVP_


**PREPROCESSOR**

Handles the import of pre-configured xml data files and generates the geometry of single/multi-element foil systems

**SOLVER**

Implements a Source Panel/Vortex Panel (SPVP) method to evaluate the influence coefficients and primary matrix of linear problems

**POSTPROCESSOR**

Evaluates aerodynamic coefficients, results output, and visualisation through contour and streamline plots

**SPVP**

Handles the frontend user requests
