# FUNCTION - PLT2D
# Written by: Jordan Nash
# Started   : 11/01/2024 - Transferred from MATLAB to Python
#
# PURPOSE
# - Compute the integral expression for constant strength source panels
# - Source panel strengths are constant, but can change from panel to panel
# - Geometric integral for panel-normal    : I(ij)
# - Geometric integral for panel-tangential: J(ij)
#
# INPUTS
# - XC  : X-coordinate of control points
# - YC  : Y-coordinate of control points
# - XB  : X-coordinate of boundary points
# - YB  : Y-coordinate of boundary points
# - phi : Angle between positive X-axis and interior of panel
# - S   : Length of panel
# 
# OUTPUTS
# - I   : Value of panel-normal integral (Eq. 3.163 in Anderson or Ref [1])
# - J   : Value of panel-tangential integral (Eq. 3.165 in Anderson or Ref [2])

import numpy as np
import matplotlib.pyplot as plt

def PLT2D(alpha, XB, ZB, XC, ZC, XS, ZS, Xc, Zc, CP, nx):
    
    fig, ax1 = plt.subplots(1, 1, figsize=(10, 10))

    XBP = XB*np.cos(-alpha) - ZB*np.sin(-alpha)
    ZBP = XB*np.sin(-alpha) + ZB*np.cos(-alpha)
    XCP = XC*np.cos(-alpha) - ZC*np.sin(-alpha)
    ZCP = XC*np.sin(-alpha) + ZC*np.cos(-alpha)
    XSP = XS*np.cos(-alpha) - ZS*np.sin(-alpha)
    ZSP = XS*np.sin(-alpha) + ZS*np.cos(-alpha)
    XcP = Xc*np.cos(-alpha) - Zc*np.sin(-alpha)
    ZcP = Xc*np.sin(-alpha) + Zc*np.cos(-alpha)

    ax1.plot(XBP, ZBP, color='black', label='Panel')
    ax1.scatter(XCP, ZCP, color='black', label='Collocation')
    ax1.plot(XSP, ZSP, color='black', label='Chord')
    ax1.plot(XcP, ZcP, color='black', linestyle='dashed', label='Camber')

    ax2 = ax1.twinx()
    ax2.scatter(XC[nx:2*nx], CP[nx:2*nx], marker='^', color='red', label='Upper Surface')
    ax2.scatter(XC[0:nx], CP[0:nx], marker='^', color='blue', label='Lower Surface')
    ax2.invert_yaxis()

    ax1.set_title('AIRFOIL')
    ax1.legend(loc='upper right')
    ax2.legend(loc='upper center')
    ax1.axis("equal")

    plt.show()