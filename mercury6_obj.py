import os
import matplotlib.pyplot as Mplot

"""
variables
"""

dirpath = "/home/xi-feng/NCU/Meeting/MERCURY/mercury6_smirik/"
obj = "1999_US3"
classfile = ".aei"
jd = 2460800.5

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
        print(raw_value)
        DATA["t"].append(float(raw_value[0]))
        DATA["a"].append(float(raw_value[3]))
        DATA["e"].append(float(raw_value[4]))
        DATA["i"].append(float(raw_value[5]))
        DATA["Ω"].append(float(raw_value[6]))
        DATA["ω"].append(float(raw_value[7]))    
    
    rawfile.close()
    return DATA

def plot_orbital_element(rawdata):
    os.chdir("/home/xi-feng/NCU/Meeting/Orbital_Result")
    orbit_t = rawdata["t"]
    print("finish loading timescale")
    orbit_a = rawdata["a"]
    print("finish loading semimajor axis")
    orbit_e = rawdata["e"]
    print("finish loading eccentricity")
    orbit_i = rawdata["i"]
    print("finish loading inclination")
    orbit_Ω = rawdata["Ω"]
    print("finish loading longtitude of accending node")
    orbit_ω = rawdata["ω"]
    print("finish loading argument of perihelion")
    
    fig = Mplot.figure(figsize = (40, 10))
    ax_a, ax_e, ax_i, ax_Om, ax_om = fig.subplots(5, 1)
    ax_om.set_xlabel("Evolution time (yr)")
    ax_a.set_ylabel("a (AU)")
    ax_a.set_ylim(2.6, 2.8)
    ax_e.set_ylabel("e")
    ax_i.set_ylabel("i (deg)")
    ax_Om.set_ylabel("Ω (deg)")
    ax_om.set_ylabel("ω (deg)")
    ax_a.scatter(orbit_t, orbit_a, s = 0.25, c = "black")
    print("finish scattering 1/5")
    ax_e.scatter(orbit_t, orbit_e, s = 0.25, c = "black")
    print("finish scattering 2/5")
    ax_i.scatter(orbit_t, orbit_i, s = 0.25, c = "black")
    print("finish scattering 3/5")
    ax_Om.scatter(orbit_t, orbit_Ω, s = 0.25, c = "black")
    print("finish scattering 4/5")
    ax_om.scatter(orbit_t, orbit_ω, s = 0.25, c = "black")
    print("finish scattering 5/5")
    
    Mplot.savefig(f"para_data_mercury_jd{jd}.png")
    
plot_orbital_element(orbital_element(dirpath + obj + classfile))


