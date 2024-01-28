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
import sys
import keyboard
import matplotlib.pyplot as plt
import numpy as np

from FOIL import BOUNDPTS
from FOIL import COLLOPTS

# %% READ XML DATA FILE
# - Reads the pre-configured single/multi-element geometry from an xml file
# - Allows for a simplfied means of defining multi-elements
# - Reads flow data 
# - Generates a 'file' from the xml

def XML(file):
    base_path = r"XML"
    file_path = os.path.join(base_path, file)
    
    with open(file_path) as file:
        content = file.read()

    order = ["<slat>", "<wing>", "<flap>"]
    
    geometry_data = {key: [] for key in order if key in content}

    for key in geometry_data:
        section = content.split(key)[1:]
        for i, info in enumerate(section, start=1):
            geometry_data[key].append({
                "foil": info.split("<foil>")[1].split("</foil>")[0],
                "panels": int(info.split("<panels>")[1].split("</panels>")[0]),
                "chord": float(info.split("<chord>")[1].split("</chord>")[0]),
                "camber": float(info.split("<camber>")[1].split("</camber>")[0]),
                "camberpos": float(info.split("<camberpos>")[1].split("</camberpos>")[0]),
                "thick": float(info.split("<thick>")[1].split("</thick>")[0]),
                "angle": np.radians(float(info.split("<angle>")[1].split("</angle>")[0])),
                "Xgap": float(info.split("<Xgap>")[1].split("</Xgap>")[0]),
                "Zgap": float(info.split("<Zgap>")[1].split("</Zgap>")[0]),
                "name": f"{key}{i}"
            })
            
    flow_data = ({
        "Vinf": float(info.split("<velocity>")[1].split("</velocity>")[0]),
        "rho": float(info.split("<density>")[1].split("</density>")[0]),
        "mu": float(info.split("<viscosity>")[1].split("</viscosity>")[0]),
        "bi": float(info.split("<boundary>")[1].split("</boundary>")[0]),
        "name": "<flow>"
    })
    
    return geometry_data, flow_data


# %% COMPUTE GEOMETRY
# - Using the 'file' dictionary the discretised geometry of the multi-element is generated
# - Individual foil elements are stored in 'geometry' dictionary
# - Complete foil elements are stored in 'system' dictionary (includes concatenated arrays with false panels)

def GEOMETRY(geometry_data): 
    
    # Ensure that the file is ordered correctly
    order = ["<slat>", "<wing>", "<flap>"]
    keys = [key for key in order if key in geometry_data]
    individual = {key: [] for key in order if key in geometry_data}
    
    Xpivot = 0
    Zpivot = 0
    
    # Generate individual element dictionary "individual"
    for j in range(len(keys)):
        for i in range(len(geometry_data[keys[j]])):
            if j == i == 0:
                offset = 0
                cumulative_offset = 0
            XB, ZB, Xchord, Zchord, Xcamb, Zcamb, offset = BOUNDPTS(0 if keys[j] == "<wing>" else geometry_data[keys[j]][i]["angle"],geometry_data[keys[j]][i]["panels"]+1, geometry_data[keys[j]][i]["chord"], geometry_data[keys[j]][i]["camber"], geometry_data[keys[j]][i]["camberpos"], geometry_data[keys[j]][i]["thick"], geometry_data[keys[j]][i]["Xgap"], geometry_data[keys[j]][i]["Zgap"], offset)
            
            start_ind = cumulative_offset
            end_ind = cumulative_offset + len(XB) - 2
            cumulative_offset = end_ind + 2 
            
            if i+j == (len(keys)-1)+(len(geometry_data[keys[j]])-1):
                length = len(XB)-1
                end_ind = end_ind + 1
            else:
                length = len(XB)
                
            if keys[j] == '<wing>':
                Xpivot = XB[int(len(XB)/2)]
                Zpivot = ZB[int(len(ZB)/2)]
            
            individual[keys[j]].append({
                "XB": XB,
                "ZB": ZB,
                "Xchord": Xchord,
                "Zchord": Zchord,
                "Xcamb": Xcamb,
                "Zcamb": Zcamb,
                "alpha": geometry_data[keys[j]][i]["angle"]*np.ones(length),
                "start_ind": start_ind,
                "end_ind": end_ind,
                "offset": offset,
                "foil": geometry_data[keys[j]][i]["foil"],
                "nx": geometry_data[keys[j]][i]["panels"],
                "c": geometry_data[keys[j]][i]["chord"],
                "t": geometry_data[keys[j]][i]["thick"],
                "m": geometry_data[keys[j]][i]["camber"],
                "p": geometry_data[keys[j]][i]["camberpos"],
                "name": f"{keys[j]}{i}"
            })
        
    # Generate concatenated multi-element dictionary "multi"
    multi = {}
    for key in individual[keys[0]][0].keys():
        arrays = [np.atleast_1d(category_dict[key]) for category_name, category_list in individual.items() for category_dict in category_list]
        concatenated_array = np.concatenate(arrays)
        multi[key] = concatenated_array
       
    XC, ZC, S, beta, phi = COLLOPTS(multi)
        
    multi.update({
        "XC": XC,
        "ZC": ZC,
        "S": S,
        "beta": beta,
        "phi": phi,
        "Xpivot": Xpivot,
        "Zpivot": Zpivot
        })

    return individual, multi

# %% PLOT GEOMETRY
#
#

def PLTGEOMETRY(individual,multi):
    
    plt.figure(figsize=(9, 9))
    
    # Plot the collocation points
    XB = multi["XB"]
    ZB = multi["ZB"]
    XC = multi["XC"]
    ZC = multi["ZC"]
    Xchord = multi["Xchord"]
    Zchord = multi["Zchord"]
    Xcamb = multi["Xcamb"]
    Zcamb = multi["Zcamb"]
    Xpivot = multi["Xpivot"]
    Zpivot = multi["Zpivot"]
    alpha = individual["<wing>"][0]["alpha"][0]
    start_ind = multi["start_ind"]
    end_ind = multi["end_ind"]

    XC = ((XC-Xpivot)*np.cos(-alpha) - (ZC-Zpivot)*np.sin(-alpha)) + Xpivot
    ZC = ((XC-Xpivot)*np.sin(-alpha) + (ZC-Zpivot)*np.cos(-alpha)) + Zpivot
    
    XB = ((XB-Xpivot)*np.cos(-alpha) - (ZB-Zpivot)*np.sin(-alpha)) + Xpivot
    ZB = ((XB-Xpivot)*np.sin(-alpha) + (ZB-Zpivot)*np.cos(-alpha)) + Zpivot
    
    Xchord = ((Xchord-Xpivot)*np.cos(-alpha) - (Zchord-Zpivot)*np.sin(-alpha)) + Xpivot
    Zchord = ((Xchord-Xpivot)*np.sin(-alpha) + (Zchord-Zpivot)*np.cos(-alpha)) + Zpivot
    
    Xcamb = ((Xcamb-Xpivot)*np.cos(-alpha) - (Zcamb-Zpivot)*np.sin(-alpha)) + Xpivot
    Zcamb = ((Xcamb-Xpivot)*np.sin(-alpha) + (Zcamb-Zpivot)*np.cos(-alpha)) + Zpivot

    for i in range(len(XB)-2):
        label = 'Real Panel' if any(start <= i <= end for start, end in zip(start_ind, end_ind)) else 'False Panel'
        c = 'blue' if any(start <= i <= end for start, end in zip(start_ind, end_ind)) else 'red'
        plt.plot([XB[i], XB[i+2]], [ZB[i], ZB[i+2]], color=c, label=label)

    plt.scatter(XC, ZC, color='black', linestyle='dashed', label='Collocation')
    plt.plot(Xchord, Zchord, color='black', linestyle='dotted', label='Chord')
    plt.plot(Xcamb, Zcamb, color='black', linestyle='dashed', label='Camber')
    plt.axis("equal")
    
    # Legend outside the loop to avoid multiple entries
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    plt.title('SOURCE PANEL VORTEX PANEL GEOMETRY [SPVP]')
    plt.show()
    
# %% OUTPUT TO CONSOLE
#
#

def PRECONSOLE(flow_data, multi, duration):
    print("")
    print("")
    print("")
    print("================= PREPROCESSOR COMPLETE  =================")    
    print("")
    print("COMPLETED IN : %2.3f" % duration + " sec")

    for i in range(len(multi["name"])):
        name = multi["name"][i]
        foil = multi["foil"][i]
        nx = multi["nx"][i]
        c = multi["c"][i]
        t = multi["t"][i]
        m = multi["m"][i]
        p = multi["p"][i]
        alpha = multi["alpha"][multi["start_ind"][i]]
        XC = multi["XC"][multi["start_ind"][i]:multi["end_ind"][i]]
        XB = multi["XB"][multi["start_ind"][i]:multi["end_ind"][i]]

        PRTFOIL(name, foil, nx, c, t, m, p, alpha, XC, XB)
        
    Vinf = flow_data["Vinf"]
    rho = flow_data["rho"]
    mu = flow_data["mu"]
    c = multi["c"][i]
        
    PRTFLOW(Vinf, rho, mu, c)
        
    print("")
    print("        [Press (C) to continue or (X) to cancel]          ")
    print("==========================================================")
        
    key = keyboard.read_event(suppress=True).name
    if key == "c":
        status = False
    if key == "x":
        sys.exit()
    return status

# %% FOIL CONSOLE MESSAGE
#
#   
 
def PRTFOIL(name, foil, nx, c, t, m, p, alpha, XC, XB):
    print("")
    print("[" + name + "]                                       ")
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
    
# %% FLOW CONSOLE MESSAGE
#
# 

def PRTFLOW(Vinf, rho, mu, c):
    print("")
    print("[flow]")    
    print("__________________________________________________________")
    print(" VELOCITY [Vinf]                              : %2.3f" % Vinf)
    print(" DENSITY [rho]                                : %2.3f" % rho)
    print(" VISCOSITY [mu]                               : %2.8f" % mu)
    print(" REYNOLDS [RE]                                : %2.3f" % ((rho*Vinf*c)/mu))
