import os, sys, glob
import matplotlib.pyplot as Mplot

"""
variables
"""
raw_path = os.path.realpath(sys.argv[0])
raw_direc = os.path.dirname(raw_path)
work_direc = "MERCURY"
classfile = ".aei"

aei_classfile = os.path.join(raw_direc + f"/{work_direc}", "**", f"*{classfile}")
aei_glob      = glob.glob(aei_classfile, recursive = True)
aei_file_list = sorted([i if "000" in i else '' for i in aei_glob])

#dirpath = "/home/xi-feng/NCU/Meeting/MERCURY/mercury6_smirik/"
#dirpath = "/home/xi-feng/NCU/Meeting/Orbital_Result/mercury_data_800kyr/"

jd = 2453311.5

def orbital_element(textfile):
    rawfile = open(textfile, "rt", encoding = "utf-8")
    rawlist = rawfile.readlines()
    para    = rawlist[3]
    
    DATA    = {}
    DATA["t"] = []
    DATA["a"] = []
    DATA["e"] = []
    DATA["i"] = []
    DATA["Ω"] = []
    DATA["ω"] = []
    
    for data in range(4, len(rawlist)):
        raw_value = rawlist[data].split(" ")
        for k in range(len(raw_value)):
            if '' in raw_value:
                raw_value.remove('')
            else:
                continue    
        try:
            DATA["a"].append(float(raw_value[3]))
            DATA["e"].append(float(raw_value[4]))
            DATA["i"].append(float(raw_value[5]))
            DATA["Ω"].append(float(raw_value[6]))
            DATA["ω"].append(float(raw_value[7]))
            DATA["t"].append(float(raw_value[0]) / 1e6)
        except ValueError:
            pass   
    
    rawfile.close()
    return DATA

def rec_multi_aeifile():
    MULTIDATA = {}
    for aeif in aei_file_list:
        if aeif == '':
            continue
        else:
            bodyid = aeif.split("/")[-1].replace(".aei", "")
            MULTIDATA[bodyid] = orbital_element(aeif)
            print("finish loading", bodyid + classfile)
    return MULTIDATA        

def plot_orbital_element():
    fig = Mplot.figure(figsize = (10, 10))
    ax_a, ax_e, ax_Om, ax_om, ax_i = fig.subplots(5, 1)
    ax_a.set_title(f"10 clones of Eros's orbital evolution")
    ax_i.set_xlabel("Evolution time (Myr)")
    ax_a.set_ylabel("a (AU)")
    #ax_a.set_yscale("log")
    ax_e.set_ylabel("e")
    ax_i.set_ylabel("i (deg)")
    ax_Om.set_ylabel("q (AU)")
    ax_om.set_ylabel("Q (AU)")
    
    rawdata = rec_multi_aeifile()
    for rd in rawdata.keys():
        print(f"case{rd}: {(min(rawdata[rd]['a']), max(rawdata[rd]['a']))}")
    # saved case
    # [5, 19, 22, 26, 27, 34, 41, 61, 62, 98]
    
    caselist = [list(rawdata.keys())[ID] for ID in [5, 19, 22, 26, 27, 34, 41, 61, 62, 98]]
    
    for data in caselist:
        print("loading", data + classfile, "...")
        orbit_t = rawdata[data]["t"]
        print("finish loading timescale", f"{len(orbit_t)} values")
        orbit_a = rawdata[data]["a"]
        print("finish loading semimajor axis", f"{len(orbit_a)} values")
        orbit_e = rawdata[data]["e"]
        print("finish loading eccentricity", f"{len(orbit_e)} values")
        orbit_i = rawdata[data]["i"]
        print("finish loading inclination", f"{len(orbit_i)} values")
        orbit_Ω = rawdata[data]["Ω"]
        print("finish loading longtitude of accending inode", f"{len(orbit_Ω)} values")
        orbit_ω = rawdata[data]["ω"]
        print("finish loading argument of perihelion", f"{len(orbit_ω)} values") 
        
        orbit_q = [orbit_a[i] * (1 - orbit_e[i]) for i in range(len(orbit_a))]
        orbit_Q = [orbit_a[i] * (1 + orbit_e[i]) for i in range(len(orbit_a))]
        
        if data == "000000":
            continue
        else:
            ax_a.plot(orbit_t, orbit_a, linewidth = 1, alpha = 0.75)
            print("finish scattering 1/5")
            ax_e.plot(orbit_t, orbit_e, linewidth = 1, alpha = 0.75)
            print("finish scattering 2/5")
            ax_i.plot(orbit_t, orbit_i, linewidth = 1, alpha = 0.75)
            print("finish scattering 3/5")
            ax_Om.plot(orbit_t, orbit_q, linewidth = 1, alpha = 0.75)
            print("finish scattering 4/5")
            ax_om.plot(orbit_t, orbit_Q, linewidth = 1, alpha = 0.75)
            print("finish scattering 5/5")    
    
    init_data = "000000"
    ax_a.plot(rawdata[init_data]["t"], rawdata[init_data]["a"], color = "black", linewidth = 0.75, alpha = 0.7, label = "original data")
    print("finish scattering 1/5")
    ax_e.plot(rawdata[init_data]["t"], rawdata[init_data]["e"], color = "black", linewidth = 0.75, alpha = 0.7)
    print("finish scattering 2/5")
    ax_i.plot(rawdata[init_data]["t"], rawdata[init_data]["i"], color = "black", linewidth = 0.75, alpha = 0.7)
    print("finish scattering 3/5")
    ax_Om.plot(rawdata[init_data]["t"], [rawdata[init_data]["a"][i] * (1 - rawdata[init_data]["e"][i]) for i in range(len(rawdata[init_data]["a"]))], 
               color = "black", linewidth = 0.75, alpha = 0.7)
    print("finish scattering 4/5")
    ax_om.plot(rawdata[init_data]["t"], [rawdata[init_data]["a"][i] * (1 + rawdata[init_data]["e"][i]) for i in range(len(rawdata[init_data]["a"]))],
               color = "black", linewidth = 0.75, alpha = 0.7)
    print("finish scattering 5/5")    
    
    ax_a.legend()
    Mplot.savefig(f"para_data_mercury_jd{jd}.png")
    Mplot.close()
    
    figai = Mplot.figure(figsize = (10, 10))
    axai = figai.add_subplot()
    axai.set_title("a-i relation of Eros's orbital evolution with 10 clones")
    axai.set_xlabel("a (AU)")
    axai.set_ylabel("i (deg)")
    
    quartor = int(len(rawdata[init_data]["t"]) / 4)
    for cp in caselist:
        axai.scatter([rawdata[cp]["a"][i] for i in range(quartor)], [rawdata[cp]["i"][i] for i in range(quartor)], s = 0.5, c = "blue")
        axai.scatter([rawdata[cp]["a"][i] for i in range(quartor, 2*quartor)], [rawdata[cp]["i"][i] for i in range(quartor, 2*quartor)], s = 0.5, c = "red")
        axai.scatter([rawdata[cp]["a"][i] for i in range(2*quartor, 3*quartor)], [rawdata[cp]["i"][i] for i in range(2*quartor, 3*quartor)], 
                     s = 0.5, c = "orange")
        axai.scatter([rawdata[cp]["a"][i] for i in range(3*quartor, 4*quartor)], [rawdata[cp]["i"][i] for i in range(3*quartor, 4*quartor)], 
                     s = 0.5, c = "green")
    axai.legend(["0 ~ τ/4", "τ/4 ~ τ/2", "τ/2 ~ 3τ/4", "3τ/4 ~ τ"])
    Mplot.savefig(f"paraAvsI_jd{jd}_mercury.png")
    Mplot.close()
    
plot_orbital_element()


