# FUNCTION - PREPROCESS
# Written by: Jordan Nash
# Started   : 24/01/2024
# Updated   :
# 
#
# PURPOSE
# - Imports geometry and flow data from XML format
#
# %% IMPORT LIBRARIES AND FUNCTIONS
#
#

import os
import matplotlib.pyplot as plt
import numpy as np
from FOIL import FOIL
# %% READ XML DATA FILE
#
#

def XML(file):
    base_path = r"F:\WORK\DSTG\Panel Methods\2D\Source Vortex\VER3\XML"
    file_path = os.path.join(base_path, file)
    
    with open(file_path) as file:
        content = file.read()

    # Extract wing parameters
    wing = content.split("<wing>")[1:]
    wing_data = []
    for wing_content in wing:
        wing_data.append({
            "foil": wing_content.split("<foil>")[1].split("</foil>")[0],
            "panels": int(wing_content.split("<panels>")[1].split("</panels>")[0]),
            "chord": float(wing_content.split("<chord>")[1].split("</chord>")[0]),
            "camber": float(wing_content.split("<camber>")[1].split("</camber>")[0]),
            "camberpos": float(wing_content.split("<camberpos>")[1].split("</camberpos>")[0]),
            "thick": float(wing_content.split("<thick>")[1].split("</thick>")[0]),
            "angle": np.radians(float(wing_content.split("<angle>")[1].split("</angle>")[0])),
            "Xgap": float(wing_content.split("<Xgap>")[1].split("</Xgap>")[0]),
            "Zgap": float(wing_content.split("<Zgap>")[1].split("</Zgap>")[0]),
            "wing_name": "wing"
        })

    # Extract flap parameters
    flaps = content.split("<flap>")[1:]
    flap_data = []
    for i, flap_content in enumerate(flaps, start=1):
        flap_data.append({
            "foil": flap_content.split("<foil>")[1].split("</foil>")[0],
            "panels": int(flap_content.split("<panels>")[1].split("</panels>")[0]),
            "chord": float(flap_content.split("<chord>")[1].split("</chord>")[0]),
            "camber": float(flap_content.split("<camber>")[1].split("</camber>")[0]),
            "camberpos": float(flap_content.split("<camberpos>")[1].split("</camberpos>")[0]),
            "thick": float(flap_content.split("<thick>")[1].split("</thick>")[0]),
            "angle": np.radians(float(flap_content.split("<angle>")[1].split("</angle>")[0])),
            "Xgap": float(flap_content.split("<Xgap>")[1].split("</Xgap>")[0]),
            "Zgap": float(flap_content.split("<Zgap>")[1].split("</Zgap>")[0]),
            "flap_name": f"flap{i}"
        })

    # Extract slat parameters
    slats = content.split("<slat>")[1:]
    slat_data = []
    for i, slat_content in enumerate(slats, start=1):
        slat_data.append({
            "foil": slat_content.split("<foil>")[1].split("</foil>")[0],
            "panels": int(slat_content.split("<panels>")[1].split("</panels>")[0]),
            "chord": float(slat_content.split("<chord>")[1].split("</chord>")[0]),
            "camber": float(slat_content.split("<camber>")[1].split("</camber>")[0]),
            "camberpos": float(slat_content.split("<camberpos>")[1].split("</camberpos>")[0]),
            "thick": float(slat_content.split("<thick>")[1].split("</thick>")[0]),
            "angle": np.radians(float(slat_content.split("<angle>")[1].split("</angle>")[0])),
            "Xgap": float(slat_content.split("<Xgap>")[1].split("</Xgap>")[0]),
            "Zgap": float(slat_content.split("<Zgap>")[1].split("</Zgap>")[0]),
            "slat_name": f"slat{i}"
        })
        
    # Extract flow parameters
    flow = content.split("<flow>")[1:]
    flow_data = []
    for i, flow_content in enumerate(flow, start=1):
        flow_data.append({
            "velocity": float(slat_content.split("<velocity>")[1].split("</velocity>")[0]),
            "density": float(slat_content.split("<density>")[1].split("</density>")[0]),
            "viscosity": float(slat_content.split("<viscosity>")[1].split("</viscosity>")[0]),
        })
    
    return {"wing": wing_data, "flap": flap_data, "slat": slat_data, "flow": flow_data}
# %% COMPUTE GEOMETRY
#
#

def GEOMETRY(file):    
    slat_data = []
    for i in range(len(file["slat"])):
        if i == 0:
            offset = 0
        XB, ZB, XC, ZC, Xchord, Zchord, Xcamb, Zcamb, offset, alpha, beta, phi, S = FOIL(file["slat"][i]["panels"], file["slat"][i]["chord"], file["slat"][i]["camber"], file["slat"][i]["camberpos"], file["slat"][i]["thick"], file["slat"][i]["angle"], file["slat"][i]["Xgap"], file["slat"][i]["Zgap"], offset, 0, 0)
        slat_data.append({
                 "XB": XB,
                 "ZB": ZB,
                 "XC": XC,
                 "ZC": ZC,
                 "Xchord": Xchord,
                 "Zchord": Zchord,
                 "Xcamb": Xcamb,
                 "Zcamb": Zcamb,
                 "alpha": alpha,
                 "beta": beta,
                 "phi": phi,
                 "slat_name": f"slat{i}"
             })
    wing_data = []
    for i in range(len(file["wing"])):
        XB, ZB, XC, ZC, Xchord, Zchord, Xcamb, Zcamb, offset, alpha, beta, phi, S = FOIL(file["wing"][i]["panels"], file["wing"][i]["chord"], file["wing"][i]["camber"], file["wing"][i]["camberpos"], file["wing"][i]["thick"], file["wing"][i]["angle"], file["wing"][i]["Xgap"], file["wing"][i]["Zgap"], offset, 0, 0)
        wing_data.append({
                 "XB": XB,
                 "ZB": ZB,
                 "XC": XC,
                 "ZC": ZC,
                 "Xchord": Xchord,
                 "Zchord": Zchord,
                 "Xcamb": Xcamb,
                 "Zcamb": Zcamb,
                 "alpha": alpha,
                 "beta": beta,
                 "phi": phi,
                 "wing_name": f"wing{i}"
             })
    flap_data = []
    for i in range(len(file["flap"])):
        XB, ZB, XC, ZC, Xchord, Zchord, Xcamb, Zcamb, offset, alpha, beta, phi, S = FOIL(file["flap"][i]["panels"], file["flap"][i]["chord"], file["flap"][i]["camber"], file["flap"][i]["camberpos"], file["flap"][i]["thick"], file["flap"][i]["angle"], file["flap"][i]["Xgap"], file["flap"][i]["Zgap"], offset, 0, 0)
        flap_data.append({
                 "XB": XB,
                 "ZB": ZB,
                 "XC": XC,
                 "ZC": ZC,
                 "Xchord": Xchord,
                 "Zchord": Zchord,
                 "Xcamb": Xcamb,
                 "Zcamb": Zcamb,
                 "alpha": alpha,
                 "beta": beta,
                 "phi": phi,
                 "flap_name": f"flap{i}"
             })    
    return {"slat": slat_data, "wing": wing_data, "flap": flap_data}
# %% PLOT GEOMETRY
#
#

def PLTGEOMETRY(geometry):
    
    plt.figure(figsize=(8, 8))
    
    # Plot the slat/s
    for i in range(len(geometry["slat"])):
        XB = geometry["slat"][i]["XB"]
        ZB = geometry["slat"][i]["ZB"]
        XC = geometry["slat"][i]["XC"]
        ZC = geometry["slat"][i]["ZC"]
        Xchord = geometry["slat"][i]["Xchord"]
        Zchord = geometry["slat"][i]["Zchord"]
        Xcamb = geometry["slat"][i]["Xcamb"]
        Zcamb = geometry["slat"][i]["Zcamb"]
        alpha = geometry["slat"][i]["alpha"]
        
        # Localised pivot
        Xpivot = XB[int((len(XC)/2)+1)]
        Zpivot = ZB[int((len(ZC)/2)+1)]
        
        XB = ((XB-Xpivot)*np.cos(-alpha) - (ZB-Zpivot)*np.sin(-alpha)) + Xpivot
        ZB = ((XB-Xpivot)*np.sin(-alpha) + (ZB-Zpivot)*np.cos(-alpha)) + Zpivot
        
        XC = ((XC-Xpivot)*np.cos(-alpha) - (ZC-Zpivot)*np.sin(-alpha)) + Xpivot
        ZC = ((XC-Xpivot)*np.sin(-alpha) + (ZC-Zpivot)*np.cos(-alpha)) + Zpivot
        
        Xchord = ((Xchord-Xpivot)*np.cos(-alpha) - (Zchord-Zpivot)*np.sin(-alpha)) + Xpivot
        Zchord = ((Xchord-Xpivot)*np.sin(-alpha) + (Zchord-Zpivot)*np.cos(-alpha)) + Zpivot
        
        Xcamb = ((Xcamb-Xpivot)*np.cos(-alpha) - (Zcamb-Zpivot)*np.sin(-alpha)) + Xpivot
        Zcamb = ((Xcamb-Xpivot)*np.sin(-alpha) + (Zcamb-Zpivot)*np.cos(-alpha)) + Zpivot
        
        # Globalised pivot
        XBW = geometry["wing"][0]["XB"]
        ZBW = geometry["wing"][0]["ZB"]
        XCW = geometry["wing"][0]["XC"]
        ZCW = geometry["wing"][0]["ZC"]
        alphaW = geometry["wing"][0]["alpha"]
        Xpivot = XBW[int((len(XCW)/2)+1)]
        Zpivot = ZBW[int((len(ZCW)/2)+1)]
      
        XB = ((XB-Xpivot)*np.cos(-alphaW) - (ZB-Zpivot)*np.sin(-alphaW)) + Xpivot
        ZB = ((XB-Xpivot)*np.sin(-alphaW) + (ZB-Zpivot)*np.cos(-alphaW)) + Zpivot
        
        XC = ((XC-Xpivot)*np.cos(-alphaW) - (ZC-Zpivot)*np.sin(-alphaW)) + Xpivot
        ZC = ((XC-Xpivot)*np.sin(-alphaW) + (ZC-Zpivot)*np.cos(-alphaW)) + Zpivot
        
        Xchord = ((Xchord-Xpivot)*np.cos(-alphaW) - (Zchord-Zpivot)*np.sin(-alphaW)) + Xpivot
        Zchord = ((Xchord-Xpivot)*np.sin(-alphaW) + (Zchord-Zpivot)*np.cos(-alphaW)) + Zpivot
        
        Xcamb = ((Xcamb-Xpivot)*np.cos(-alphaW) - (Zcamb-Zpivot)*np.sin(-alphaW)) + Xpivot
        Zcamb = ((Xcamb-Xpivot)*np.sin(-alphaW) + (Zcamb-Zpivot)*np.cos(-alphaW)) + Zpivot
        
        plt.plot(XB, ZB, color='black')
        plt.scatter(XC, ZC, color='black')
        plt.plot(Xchord, Zchord, color='black')
        plt.plot(Xcamb, Zcamb, color='black', linestyle='dashed')
        plt.axis("equal")
    
    # Plot the wing
    XB = geometry["wing"][0]["XB"]
    ZB = geometry["wing"][0]["ZB"]
    XC = geometry["wing"][0]["XC"]
    ZC = geometry["wing"][0]["ZC"]
    Xchord = geometry["wing"][0]["Xchord"]
    Zchord = geometry["wing"][0]["Zchord"]
    Xcamb = geometry["wing"][0]["Xcamb"]
    Zcamb = geometry["wing"][0]["Zcamb"]
    alpha = geometry["wing"][0]["alpha"]
    
    Xpivot = XB[int((len(XC)/2)+1)]
    Zpivot = ZB[int((len(ZC)/2)+1)]
    
    XB = ((XB-Xpivot)*np.cos(-alpha) - (ZB-Zpivot)*np.sin(-alpha)) + Xpivot
    ZB = ((XB-Xpivot)*np.sin(-alpha) + (ZB-Zpivot)*np.cos(-alpha)) + Zpivot
    
    XC = ((XC-Xpivot)*np.cos(-alpha) - (ZC-Zpivot)*np.sin(-alpha)) + Xpivot
    ZC = ((XC-Xpivot)*np.sin(-alpha) + (ZC-Zpivot)*np.cos(-alpha)) + Zpivot
    
    Xchord = ((Xchord-Xpivot)*np.cos(-alpha) - (Zchord-Zpivot)*np.sin(-alpha)) + Xpivot
    Zchord = ((Xchord-Xpivot)*np.sin(-alpha) + (Zchord-Zpivot)*np.cos(-alpha)) + Zpivot
    
    Xcamb = ((Xcamb-Xpivot)*np.cos(-alpha) - (Zcamb-Zpivot)*np.sin(-alpha)) + Xpivot
    Zcamb = ((Xcamb-Xpivot)*np.sin(-alpha) + (Zcamb-Zpivot)*np.cos(-alpha)) + Zpivot

    plt.plot(XB, ZB, color='black', label='Panel')
    plt.scatter(XC, ZC, color='black', label='Collocation')
    plt.plot(Xchord, Zchord, color='black', label='Chord')
    plt.plot(Xcamb, Zcamb, color='black', linestyle='dashed', label='Camber')
    plt.axis("equal")
    plt.title('SOURCE PANEL VORTEX PANEL GEOMETRY [SPVP]')
    plt.legend(loc='upper right')
        
    # Plot the flap/s
    for i in range(len(geometry["flap"])):
        XB = geometry["flap"][i]["XB"]
        ZB = geometry["flap"][i]["ZB"]
        XC = geometry["flap"][i]["XC"]
        ZC = geometry["flap"][i]["ZC"]
        Xchord = geometry["flap"][i]["Xchord"]
        Zchord = geometry["flap"][i]["Zchord"]
        Xcamb = geometry["flap"][i]["Xcamb"]
        Zcamb = geometry["flap"][i]["Zcamb"]
        alpha = geometry["flap"][i]["alpha"]
        
        # Localised pivot
        Xpivot = XB[int((len(XC)/2)+1)]
        Zpivot = ZB[int((len(ZC)/2)+1)]
        
        XB = ((XB-Xpivot)*np.cos(-alpha) - (ZB-Zpivot)*np.sin(-alpha)) + Xpivot
        ZB = ((XB-Xpivot)*np.sin(-alpha) + (ZB-Zpivot)*np.cos(-alpha)) + Zpivot
        
        XC = ((XC-Xpivot)*np.cos(-alpha) - (ZC-Zpivot)*np.sin(-alpha)) + Xpivot
        ZC = ((XC-Xpivot)*np.sin(-alpha) + (ZC-Zpivot)*np.cos(-alpha)) + Zpivot
        
        Xchord = ((Xchord-Xpivot)*np.cos(-alpha) - (Zchord-Zpivot)*np.sin(-alpha)) + Xpivot
        Zchord = ((Xchord-Xpivot)*np.sin(-alpha) + (Zchord-Zpivot)*np.cos(-alpha)) + Zpivot
        
        Xcamb = ((Xcamb-Xpivot)*np.cos(-alpha) - (Zcamb-Zpivot)*np.sin(-alpha)) + Xpivot
        Zcamb = ((Xcamb-Xpivot)*np.sin(-alpha) + (Zcamb-Zpivot)*np.cos(-alpha)) + Zpivot
        
        # Globalised pivot
        XBW = geometry["wing"][0]["XB"]
        ZBW = geometry["wing"][0]["ZB"]
        XCW = geometry["wing"][0]["XC"]
        ZCW = geometry["wing"][0]["ZC"]
        alphaW = geometry["wing"][0]["alpha"]
        Xpivot = XBW[int((len(XCW)/2)+1)]
        Zpivot = ZBW[int((len(ZCW)/2)+1)]
      
        XB = ((XB-Xpivot)*np.cos(-alphaW) - (ZB-Zpivot)*np.sin(-alphaW)) + Xpivot
        ZB = ((XB-Xpivot)*np.sin(-alphaW) + (ZB-Zpivot)*np.cos(-alphaW)) + Zpivot
        
        XC = ((XC-Xpivot)*np.cos(-alphaW) - (ZC-Zpivot)*np.sin(-alphaW)) + Xpivot
        ZC = ((XC-Xpivot)*np.sin(-alphaW) + (ZC-Zpivot)*np.cos(-alphaW)) + Zpivot
        
        Xchord = ((Xchord-Xpivot)*np.cos(-alphaW) - (Zchord-Zpivot)*np.sin(-alphaW)) + Xpivot
        Zchord = ((Xchord-Xpivot)*np.sin(-alphaW) + (Zchord-Zpivot)*np.cos(-alphaW)) + Zpivot
        
        Xcamb = ((Xcamb-Xpivot)*np.cos(-alphaW) - (Zcamb-Zpivot)*np.sin(-alphaW)) + Xpivot
        Zcamb = ((Xcamb-Xpivot)*np.sin(-alphaW) + (Zcamb-Zpivot)*np.cos(-alphaW)) + Zpivot
    
        plt.plot(XB, ZB, color='black')
        plt.scatter(XC, ZC, color='black')
        plt.plot(Xchord, Zchord, color='black')
        plt.plot(Xcamb, Zcamb, color='black', linestyle='dashed')
        plt.axis("equal")
        
    plt.show()
# %% OUTPUT CONSOLE MESSAGE
#
#

def PRECON(file, geometry):
    print("")
    print("")
    print("")
    print("================= PREPROCESSOR COMPLETE  =================")    
    
    for i in range(len(file["slat"])):
        slat_name = file["slat"][i]["slat_name"]
        foil = file["slat"][i]["foil"]
        nx = file["slat"][i]["panels"]
        c = file["slat"][i]["chord"]
        t = file["slat"][i]["thick"]
        m = file["slat"][i]["camber"]
        p = file["slat"][i]["camberpos"]
        alpha = file["slat"][i]["chord"]
        XC = geometry["slat"][i]["XC"]
        XB = geometry["slat"][i]["XB"]
        
        print("")
        print("[" + slat_name + "]                                       ")
        print("__________________________________________________________")
        print(" FOIL                                         : " + foil)
        print(" CHORD PANELS [nx]                            : %2.3f" % nx)
        print(" TOTAL PANELS [nx]                            : %2.3f" % (2*nx))
        print(" CHORD [c]                                    : %2.3f" % c)
        print(" THICKNESS [t]                                : %2.3f" % t)
        print(" CAMBER [m]                                   : %2.3f" % m)
        print(" CAMBER POSITION [p]                          : %2.3f" % p)
        print(" ANGLE OF ATTACK [alpha]                      : %2.3f" % alpha)
        print(" CONTROL POINTS                               : %2.3f" % (len(XC)))
        print(" BOUNDARY POINTS                              : %2.3f" % (len(XB)))
        
    
          
    for i in range(len(file["wing"])):
        wing_name = file["wing"][i]["wing_name"]
        foil = file["wing"][i]["foil"]
        nx = file["wing"][i]["panels"]
        c = file["wing"][i]["chord"]
        t = file["wing"][i]["thick"]
        m = file["wing"][i]["camber"]
        p = file["wing"][i]["camberpos"]
        alpha = file["wing"][i]["chord"]
        XC = geometry["wing"][i]["XC"]
        XB = geometry["wing"][i]["XB"]
        
        print("")
        print("[" + wing_name + "]                                       ")
        print("__________________________________________________________")
        print(" FOIL                                         : " + foil)
        print(" CHORD PANELS [nx]                            : %2.3f" % nx)
        print(" TOTAL PANELS [nx]                            : %2.3f" % (2*nx))
        print(" CHORD [c]                                    : %2.3f" % c)
        print(" THICKNESS [t]                                : %2.3f" % t)
        print(" CAMBER [m]                                   : %2.3f" % m)
        print(" CAMBER POSITION [p]                          : %2.3f" % p)
        print(" ANGLE OF ATTACK [alpha]                      : %2.3f" % alpha)
        print(" CONTROL POINTS                               : %2.3f" % (len(XC)))
        print(" BOUNDARY POINTS                              : %2.3f" % (len(XB)))
          
          
    for i in range(len(file["flap"])):
        flap_name = file["flap"][i]["flap_name"]
        foil = file["flap"][i]["foil"]
        nx = file["flap"][i]["panels"]
        c = file["flap"][i]["chord"]
        t = file["flap"][i]["thick"]
        m = file["flap"][i]["camber"]
        p = file["flap"][i]["camberpos"]
        alpha = file["flap"][i]["chord"]
        XC = geometry["flap"][i]["XC"]
        XB = geometry["flap"][i]["XB"]
        
        print("")
        print("[" + flap_name + "]                                       ")
        print("__________________________________________________________")
        print(" FOIL                                         : " + foil)
        print(" CHORD PANELS [nx]                            : %2.3f" % nx)
        print(" TOTAL PANELS [nx]                            : %2.3f" % (2*nx))
        print(" CHORD [c]                                    : %2.3f" % c)
        print(" THICKNESS [t]                                : %2.3f" % t)
        print(" CAMBER [m]                                   : %2.3f" % m)
        print(" CAMBER POSITION [p]                          : %2.3f" % p)
        print(" ANGLE OF ATTACK [alpha]                      : %2.3f" % alpha)
        print(" CONTROL POINTS                               : %2.3f" % (len(XC)))
        print(" BOUNDARY POINTS                              : %2.3f" % (len(XB)))    
    
    
    for i in range(len(file["flow"])):
        Vinf = file["flow"][i]["velocity"]
        rho = file["flow"][i]["density"]
        mu = file["flow"][i]["viscosity"]
        c = file["wing"][i]["chord"]
        
        print("")
        print("[flow]")    
        print("__________________________________________________________")
        print(" VELOCITY [Vinf]                              : %2.3f" % Vinf)
        print(" DENSITY [rho]                                : %2.3f" % rho)
        print(" VISCOSITY [mu]                               : %2.8f" % mu)
        print(" REYNOLDS [RE]                                : %2.3f" % ((rho*Vinf*c)/mu))

    print("")
    print("")
    print("==========================================================")
