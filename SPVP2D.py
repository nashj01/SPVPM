# FUNCTION - SPVP2D
# Written by: Jordan Nash
# Started   : 11/01/2024
# Updated   : 24/01/2024 - "Added boundary layer effects"
# 
#
# PURPOSE
# - Compute the integral expression for constant strength source panels
# - Source panel strengths are constant, but can change from panel to panel
# - Geometric integral for panel-normal    : I(ij)
# - Geometric integral for panel-tangential: J(ij)
# _____________________________________________________________________________
#
# _____________________________________________________________________________
import matplotlib.pyplot as plt
from matplotlib import path
import numpy as np
from GEO2D import GEO2D
from COMPUTE_IJ_SPM import COMPUTE_IJ_SPM
from COMPUTE_KL_VPM import COMPUTE_KL_VPM
from STREAMLINE_SPM import STREAMLINE_SPM
from STREAMLINE_VPM import STREAMLINE_VPM
# _____________________________________________________________________________

# Initialise geometry
nx = 50
c = 1
p = 0.4
M = 0.02
t = 0.12

# Initialise flow
Vinf = 55
rho = 1.225
mu = (1.4207)*10**-5
alpha = np.radians(10)

# Pre-allocation
numPan = 2*nx
deltaC = np.zeros(numPan)
deltaB = np.zeros(numPan+1)
REX = np.zeros(numPan)
REGC = np.zeros(numPan)
REGB = np.zeros(numPan+1)
RE = (rho*Vinf*c)/(mu)

# Initialise the figure
fig, ax1 = plt.subplots(1, 1, figsize=(7,7))

for k in range(2):
    # %% COMPUTE GEOMETRY
    XB, ZB, XC, ZC, XS, ZS, dx, dz, S, beta, phi, Xc, Zc = GEO2D(alpha, deltaB, deltaC, nx+1, c, p, M, t)

    # %% COMPUTE SOURCE PANEL STRENGTHS
    I, J = COMPUTE_IJ_SPM(XC,ZC,XB,ZB,phi,S)                                        # Call COMPUTE_IJ_SPM function (Refs [2] and [3])
    
    # %% COMPUTE VORTEX PANEL STRENGTHS
    K, L = COMPUTE_KL_VPM(XC,ZC,XB,ZB,phi,S)                                        # Call COMPUTE_KL_VPM function (Refs [6] and [7])
    
    # %% COMPUTE A MATRIX
    A = np.zeros([numPan,numPan])                                                   # Initialize the A matrix
    for i in range(numPan):                                                         # Loop over all i panels
        for j in range(numPan):                                                     # Loop over all j panels
            if (i == j):                                                            # If the panels are the same
                A[i,j] = np.pi                                                      # Set A equal to pi
            else:                                                                   # If panels are not the same
                A[i,j] = I[i,j]                                                     # Set A equal to I
    
    # Right column of A matrix
    newAV = np.zeros((numPan,1))                                                    # Used to enlarge the A matrix to account for gamma column
    A     = np.hstack((A,newAV))                                                    # Horizontally stack the A matrix with newAV to get enlarged matrix
    for i in range(numPan):                                                         # Loop over all i panels (rows)
        A[i,numPan] = -sum(K[i,:])                                                  # Add gamma term to right-most column of A matrix
    
    # Bottom row of A matrix
    newAH = np.zeros((1,numPan+1))                                                  # Used to enlarge the A matrix to account for Kutta condition equation
    A     = np.vstack((A,newAH))                                                    # Vertically stack the A matrix with newAH to get enlarged matrix
    for j in range(numPan):                                                         # Loop over all j panels (columns)
        A[numPan,j] = J[0,j] + J[numPan-1,j]                                        # Source contribution of Kutta condition equation
    A[numPan,numPan] = -(sum(L[0,:] + L[numPan-1,:])) + 2*np.pi                     # Vortex contribution of Kutta condition equation 
    
    # %% COMPUTE B ARRAY
    b = np.zeros(numPan)                                                            # Initialize the b array
    for i in range(numPan):                                                         # Loop over all i panels (rows)
        b[i] = -Vinf*2*np.pi*np.cos(beta[i])                                        # Compute RHS array
    
    # Last element of b array (Kutta condition)
    b = np.append(b,-Vinf*2*np.pi*(np.sin(beta[0]) + np.sin(beta[numPan-1])))       # Add Kutta condition equation RHS to b array
    
    # %% COMPUTE RESULT ARRAY
    resArr = np.linalg.solve(A,b)                                                   # Solve system of equation for all source strengths and single vortex strength
    
    # Separate lam and gamma values from result 
    lam   = resArr[0:len(resArr)-1]                                                 # All panel source strengths
    gamma = resArr[len(resArr)-1]                                                   # Constant vortex strength
    
    # %% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS AND REGIMES
    Vt = np.zeros(numPan)                                                           # Initialize tangential velocity
    CP = np.zeros(numPan)                                                           # Initialize pressure coefficient
    
    for i in range(numPan):                                                         # Loop over all panels
        term1 = Vinf*np.sin(beta[i])                                                # Uniform flow term
        term2 = (1/(2*np.pi))*sum(lam*J[i,:])                                       # Source panel terms when j is not equal to i
        term3 = gamma/2                                                             # Vortex panel term when j is equal to i
        term4 = -(gamma/(2*np.pi))*sum(L[i,:])                                      # Vortex panel terms when j is not equal to i
        
        Vt[i] = term1 + term2 + term3 + term4                                       # Compute tangential velocity on panel i
        CP[i] = 1-(Vt[i]/Vinf)**2
                                                                                    # Compute pressure coefficient on panel i
        if k == 0:  
            if i < round(0.05*numPan):
                S[i] = S[round(0.05*numPan)]
            if i > round(0.95*numPan):
                S[i] = S[round(0.95*numPan)]
    
            REX[i] = (rho*Vt[i]*S[i])/(mu)
            
            if abs(REX[i]) >= 0.2*10**6:
                REGB[i] = 1 
                REGC[i] = 1

    # %% COMPUTE THE BOUNDARY LAYER GROWTH
    
    # Account for clockwise order of operation when dealing with boundary layer transition
    if alpha >= 0:
        for i in range(numPan):
            if i > 1:
                if REGB[i-1] == 1:
                    REGB[i] = 1
                if REGC[i-1] == 1:
                    REGC[i] = 1      
    if alpha < 0:
        for i in reversed(range(numPan)):
            if i < numPan-1:
                if REGB[i] == 1:
                    REGB[i-1] = 1
                if REGC[i] == 1:
                    REGC[i-1] = 1
            
    # Compute the boundary layer thickness growth and the skin friction coefficient            
    for i in range(numPan):     
        if REGB[i] == 0:
            deltaB[i] = (5.836/(np.sqrt(abs(REX[i]))))*XB[i]
        if REGB[i] == 1:
            deltaB[i] = (0.37/(np.power(abs(REX[i]), 1/5)))*XB[i]
            
        if REGC[i] == 0 :
            CF = (0.664/(np.sqrt(abs(REX[i]))))
            deltaC[i] = (5.836/(np.sqrt(abs(REX[i]))))*XC[i]
        if REGC[i] == 1:
            CF = (0.0576/(np.power(abs(REX[i]), 1/5)))
            deltaC[i] = (0.37/(np.power(abs(REX[i]), 1/5)))*XC[i]
                
        
# %% PLOT ORIGINAL AIRFOIL AND BOUNDARY LAYER GROWTH
    
    # Apply the 2D rotation matrix to plot rotated airfoil
    XBP = XB*np.cos(-alpha) - ZB*np.sin(-alpha)
    ZBP = XB*np.sin(-alpha) + ZB*np.cos(-alpha)
    XCP = XC*np.cos(-alpha) - ZC*np.sin(-alpha)
    ZCP = XC*np.sin(-alpha) + ZC*np.cos(-alpha)
    XSP = XS*np.cos(-alpha) - ZS*np.sin(-alpha)
    ZSP = XS*np.sin(-alpha) + ZS*np.cos(-alpha)
    XcP = Xc*np.cos(-alpha) - Zc*np.sin(-alpha)
    ZcP = Xc*np.sin(-alpha) + Zc*np.cos(-alpha)
  
    # Original geometry plot
    if k == 0:
        ax1.plot(XBP, ZBP, color='black', label='Panel')
        ax1.scatter(XCP, ZCP, color='black', label='Collocation')
        ax1.plot(XSP, ZSP, color='black', label='Chord')
        ax1.plot(XcP, ZcP, color='black', linestyle='dashed', label='Camber')
        
    # Modified geometry plot visualised as a boundary layer
    if k == 1:
        colors = ['blue' if val == 0 else 'red' for val in REGC]        
        for j in range(len(XCP) - 1):
            if j > 1:
                plt.plot([XCP[j], XCP[j + 1]], [ZCP[j], ZCP[j + 1]], color=colors[j], linewidth=4)
        ax2 = ax1.twinx()
        ax2.scatter(XC[nx:2*nx], CP[nx:2*nx], marker='^', color='red', label='Upper Surface')
        ax2.plot(XC[nx:2*nx], CP[nx:2*nx], color='red')
        ax2.scatter(XC[0:nx], CP[0:nx], marker='^', color='blue', label='Lower Surface')
        ax2.plot(XC[0:nx], CP[0:nx], color='blue')
        ax2.invert_yaxis()
        ax1.set_title('SOURCE PANEL VORTEX PANEL [SPVP]')
        ax1.legend(loc='upper right')
        ax2.legend(loc='upper center')
        ax1.axis("equal")
        plt.show()
            
# %% COMPUTE LIFT AND MOMENT COEFFICIENTS

# Compute normal and axial force coefficients
CN = -CP*S*np.sin(beta)                                                         # Normal force coefficient 
CA = -CP*S*np.cos(beta)                                                         # Axial force coefficient

# Compute lift, drag, and moment coefficients
CL = (sum(CN*np.cos(alpha)) - sum(CA*np.sin(alpha)))/c                          # Lift coefficient via the Source Panel Vortex Panel formulation
KJ = ((2*sum(gamma*S))/(Vinf))/c                                                # Lift coefficient via the Kutta Joukouski formulation
CDI = (sum(CN*np.sin(alpha)) + sum(CA*np.cos(alpha)))/(Vinf)                    # Inviscid drag coefficient via the Source Panel Vortex Panel formulation
CDV = CDI + sum(0.5*rho*Vt*S*CF)                                                # Inviscid + Viscous drag coefficient via the inclusion of a skin friction coeffcient
CM = sum(CP*(XC-0.25)*S*np.cos(phi))                                            # Moment coefficient via the Source Panel Vortex Panel formulation

# %% OUTPUT RESULTS

print("")
print("")
print("")
print("================= GEOMETRY  =================")
print("_____________________________________________")
print(" CHORD PANELS [nx]        : %2.8f" % (2*nx))
print(" CHORD [c]                : %2.8f" % c)
print(" THICKNESS [t/c]          : %2.8f" % t)
print(" CAMBER [M]               : %2.8f" % M)
print("_____________________________________________")
print(" CONTROL [NPTS]           : %2.8f" % (len(XC)))                                                  
print(" BOUNDARY [NPTS]          : %2.8f" % (len(XB))) 
print("_____________________________________________")
print("")
print("")
print("")
print("=============== SPVP RESULTS ================")
print("_____________________________________________")
print(" LIFT COEF [SPVP]         : %2.8f" % CL)                                                  
print(" LIFT COEF [K-J]          : %2.8f" % KJ) 
print("_____________________________________________")
print(" DRAG COEF [SPVP(INVCD)]  : %2.8f" % CDI)  
print(" DRAG COEF [SPVP(VISC)]   : %2.8f" % CDV)                                                  
print("_____________________________________________")
print(" MOMENT COEF [SPVP]       : %2.8f" % CM)   
print("_____________________________________________")
print(" CENTER OF PRESSURE [Xcp] : %2.8f" % ((CM/CL)+(0.25*c)))                                       
print("_____________________________________________")
print(" REYNOLDS [RE]            : %2.8f" % RE)                                       
print("_____________________________________________")



# %% COMPUTE STREAMLINES - REFS [4] and [8]


# Grid parameters
nGridX = 100                                                                # X-grid for streamlines and contours
nGridZ = 100                                                                # Y-grid for streamlines and contours
xVals  = [min(XB)-0.5, max(XB)+0.5]                                         # X-grid extents [min, max]
zVals  = [min(ZB)-0.3, max(ZB)+0.3]                                         # Y-grid extents [min, max]

# Streamline parameters
slPct  = 25                                                                 # Percentage of streamlines of the grid
Zsl    = np.linspace(zVals[0],zVals[1],int((slPct/100)*nGridZ))             # Create array of Y streamline starting points
Xsl    = xVals[0]*np.ones(len(Zsl))                                         # Create array of X streamline starting points
XYsl   = np.vstack((Xsl.T,Zsl.T)).T                                         # Concatenate X and Y streamline starting points

# Generate the grid points
Xgrid  = np.linspace(xVals[0],xVals[1],nGridX)                              # X-values in evenly spaced grid
Zgrid  = np.linspace(zVals[0],zVals[1],nGridZ)                              # Y-values in evenly spaced grid
XX, ZZ = np.meshgrid(Xgrid,Zgrid)                                           # Create meshgrid from X and Y grid arrays

# Initialize velocities
Vx     = np.zeros([nGridX,nGridZ])                                          # Initialize X velocity matrix
Vz     = np.zeros([nGridX,nGridZ])                                          # Initialize Y velocity matrix

# Path to figure out if grid point is inside polygon or not
AF     = np.vstack((XB.T,ZB.T)).T                                           # Concatenate XB and YB geometry points
afPath = path.Path(AF)                                                      # Create a path for the geometry

# Solve for grid point X and Y velocities
for m in range(nGridX):                                                     # Loop over X-grid points
    print("m: %i" % m)
    for n in range(nGridZ):                                                 # Loop over Y-grid points
        XP     = XX[m,n]                                                    # Current iteration's X grid point
        ZP     = ZZ[m,n]                                                    # Current iteration's Y grid point
        Mx, Mz = STREAMLINE_SPM(XP,ZP,XB,ZB,phi,S)                          # Compute streamline Mx and My values
        Nx, Nz = STREAMLINE_VPM(XP,ZP,XB,ZB,phi,S)                          # Compute streamline Nx and Ny values
        
        # Check if grid points are in object
        # - If they are, assign a velocity of zero
        if afPath.contains_points([(XP,ZP)]):                               # If (XP,YP) is in the body
            Vx[m,n] = 0                                                     # Set X-velocity equal to zero
            Vz[m,n] = 0                                                     # Set Y-velocity equal to zero
        else:
            Vx[m,n] = (Vinf*np.cos(alpha) + sum(lam*Mx/(2*np.pi))            # Compute X-velocity
                                         + sum(-gamma*Nx/(2*np.pi)))
            Vz[m,n] = (Vinf*np.sin(alpha) + sum(lam*Mz/(2*np.pi))            # Compute Y-velocity
                                         + sum(-gamma*Nz/(2*np.pi)))

# Compute grid point velocity magnitude and pressure coefficient
Vxz  = np.sqrt(Vx**2 + Vz**2)                                               # Compute magnitude of velocity vector []
CpXZ = 1 - (Vxz/Vinf)**2                                                    # Pressure coefficient []

            

# FIGURE: Airfoil streamlines                                                        # Create figure
plt.cla()                                                                   # Get ready for plotting
np.seterr(under="ignore")                                                   # Ignore underflow error message
plt.streamplot(XX,ZZ,Vx,Vz, linewidth=0.5, density=40, color='r',           # Plot streamlines
               arrowstyle='-', start_points=XYsl)
plt.clim(vmin=0, vmax=2)
plt.fill(XB,ZB,'k')                                                         # Plot airfoil as black polygon
plt.xlabel('X Units')                                                       # Set X-label
plt.ylabel('Z Units')                                                       # Set Y-label
plt.gca().set_aspect('equal')                                               # Set axes equal
plt.xlim(xVals)                                                             # Set X-limits
plt.ylim(zVals)                                                             # Set Y-limits
plt.show()                                                                  # Display plot

plt.cla()                                                                   # Get ready for plotting
plt.contourf(XX,ZZ,CpXZ,500,cmap='jet')                                     # Plot contour
plt.fill(XB,ZB,'k')                                                         # Plot airfoil as black polygon
plt.xlabel('X Units')                                                       # Set X-label
plt.ylabel('Z Units')                                                       # Set Y-label
plt.gca().set_aspect('equal')                                               # Set axes equal
plt.xlim(xVals)                                                             # Set X-limits
plt.ylim(zVals)                                                             # Set Y-limits
plt.show()                                                                  # Display plot