import os, sys, glob
import numpy as n
import matplotlib.pyplot as Mplot

"""
variables
"""
raw_path = os.path.realpath(sys.argv[0])
raw_direc = os.path.dirname(raw_path)
work_direc = "Orbital_Result/Eros/mercury_case2"
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
    ax_a.set_title(f"100 clones of Eros's orbital evolution")
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
    
    caselist = [list(rawdata.keys())[ID] for ID in range(len(list(rawdata.keys())))]
    
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
    Mplot.savefig(f"para_data_mercury_jd{jd}_case2.png")
    Mplot.close()
    
    figai = Mplot.figure(figsize = (10, 10))
    axai = figai.add_subplot()
    axai.set_title("a-i relation of Eros's orbital evolution with 100 clones")
    axai.set_xlabel("a (AU)")
    axai.set_ylabel("i (deg)")
    axai.set_xlim(0.7, 2)
    axai.set_ylim(0, 25)
    
    quartor = int(len(rawdata[init_data]["t"]) / 4)
    for cp in caselist:
        axai.scatter([rawdata[cp]["a"][i] for i in range(quartor)], [rawdata[cp]["i"][i] for i in range(quartor)], s = 0.5, c = "blue")
        axai.scatter([rawdata[cp]["a"][i] for i in range(quartor, 2*quartor)], [rawdata[cp]["i"][i] for i in range(quartor, 2*quartor)], s = 0.5, c = "red")
        axai.scatter([rawdata[cp]["a"][i] for i in range(2*quartor, 3*quartor)], [rawdata[cp]["i"][i] for i in range(2*quartor, 3*quartor)], 
                     s = 0.5, c = "orange")
        axai.scatter([rawdata[cp]["a"][i] for i in range(3*quartor, 4*quartor)], [rawdata[cp]["i"][i] for i in range(3*quartor, 4*quartor)], 
                     s = 0.5, c = "green")
        axai.scatter(rawdata[cp]["a"][0], rawdata[cp]["i"][0], s = 100, c = "black", marker = "*")
    axai.legend(["0 ~ τ/4", "τ/4 ~ τ/2", "τ/2 ~ 3τ/4", "3τ/4 ~ τ", "initial value"])
    Mplot.savefig(f"paraAvsI_jd{jd}_mercury_case2.png")
    Mplot.close()
    
    figae = Mplot.figure(figsize = (10, 10))
    axae = figae.add_subplot()
    axae.set_title("a-e relation of Eros's orbital evolution with 100 clones")
    axae.set_xlabel("a (AU)")
    axae.set_ylabel("e")
    for cp in caselist:
        axae.scatter([rawdata[cp]["a"][i] for i in range(quartor)], [rawdata[cp]["e"][i] for i in range(quartor)], s = 0.5, c = "blue")
        axae.scatter([rawdata[cp]["a"][i] for i in range(quartor, 2*quartor)], [rawdata[cp]["e"][i] for i in range(quartor, 2*quartor)], s = 0.5, c = "red")
        axae.scatter([rawdata[cp]["a"][i] for i in range(2*quartor, 3*quartor)], [rawdata[cp]["e"][i] for i in range(2*quartor, 3*quartor)], 
                     s = 0.5, c = "orange")
        axae.scatter([rawdata[cp]["a"][i] for i in range(3*quartor, 4*quartor)], [rawdata[cp]["e"][i] for i in range(3*quartor, 4*quartor)], 
                     s = 0.5, c = "green")
        axae.scatter(rawdata[cp]["a"][0], rawdata[cp]["e"][0], s = 100, c = "black", marker = "*")
    axae.plot(n.linspace(1.3, 2.1, 1000), [1 - 1.3 / i for i in n.linspace(1.3, 2.1, 1000)], color = "black", linestyle = "--", linewidth = 3)
    axae.plot(n.linspace(1.017, 2.1, 1000), [1 - 1.017 / i for i in n.linspace(1.017, 2.1, 1000)], color = "black", linestyle = "--", linewidth = 3)
    axae.plot(n.linspace(1.58, 2.1, 1000), [1 - 1.58 / i for i in n.linspace(1.58, 2.1, 1000)], color = "black", linestyle = "--", linewidth = 3)
    axae.plot(n.linspace(1.67, 2.1, 1000), [1 - 1.67 / i for i in n.linspace(1.67, 2.1, 1000)], color = "black", linestyle = "--", linewidth = 3)
    axae.plot(n.linspace(0.6, 0.983, 1000), [0.983 / i - 1 for i in n.linspace(0.6, 0.983, 1000)], color = "black", linestyle = "-.", linewidth = 3)
    axae.plot([1.0] * 1000, n.linspace(0.0, 0.64, 1000), color = "black", linestyle = "--", linewidth = 1.5)    
    axae.legend(["0 ~ τ/4", "τ/4 ~ τ/2", "τ/2 ~ 3τ/4", "3τ/4 ~ τ", "initial value"])
    axae.text(2.05, 1 - 1.67 / 2.1 + 0.01, "q$_{1.67}$", {"fontsize": 14})
    axae.text(2.05, 1 - 1.3 / 2.1 + 0.01, "q$_{1.3}$", {"fontsize": 14})
    axae.text(2.05, 1 - 1.58 / 2.1 + 0.01, "q$_{1.58}$", {"fontsize": 14})
    axae.text(2.05, 1 - 1.017 / 2.1 + 0.01, "q$_{1.017}$", {"fontsize": 14})
    axae.text(0.6, 0.983 / 0.6 - 1 + 0.01, "Q$_{0.983}$", {"fontsize": 14})
    Mplot.savefig(f"paraAvsE_jd{jd}_mercury_case2.png")
    Mplot.close()

def classify_element(dirname = ""):
    clone_count = 1

    rawfile_path = os.path.realpath(sys.argv[0])
    rawfile_dir  = os.path.dirname(rawfile_path)
    data_class   = sorted(glob.glob(os.path.join(rawfile_dir + f"/{dirname}", "**", f"*.aei"), recursive = True)) 
    
    # collect clone data
    data_t, data_a, data_e, data_q, data_Q = [], [], [], [], []
    for rawfile in data_class:
        if "000" in rawfile:  
            print(f"open {rawfile}...")
            data_oe = orbital_element(rawfile)
            data_t = data_oe["t"]
            data_a.append(data_oe["a"])
            data_e.append(data_oe["e"])
        else:
            continue
            
    data_q += [[data_a[a][i] * (1 - data_e[a][i]) for i in range(len(data_a[a]))] for a in range(len(data_a))]
    data_Q += [[data_a[a][i] * (1 + data_e[a][i]) for i in range(len(data_a[a]))] for a in range(len(data_a))]
    
    clone_count = len(data_a) - 4
    print("total data:", clone_count)
    
    # classification counting
    print("start classifying the data")
    q_minus, q_balance, q_plus = 0, -4, 0
    
    for cp in range(len(data_a)):
        #asteroid_class_ = []
        print(f"classify case {cp + 1}/{len(data_a)}")
       
        #counting
        if data_q[cp][0] > data_q[cp][-1] and data_q[cp][-1] < 1.017:
            q_minus += 1
        elif data_q[cp][0] < data_q[cp][-1] and data_q[cp][-1] > 1.3:
            q_plus += 1
        else:
            q_balance += 1
            continue
    
    # print the result
    print(f"(q-, q_0, q+) = {(q_minus, q_balance, q_plus)}")
    print(f"probability of ecountering earth's orbit =", round(q_minus / clone_count, 4))
    print(f"probability of the balance of orbit      =", round(q_balance / clone_count, 4))
    print(f"probability of escaping earth's orbit    =", round(q_plus / clone_count, 4))

plot_orbital_element()
#classify_element("Orbital_Result/Eros/")



