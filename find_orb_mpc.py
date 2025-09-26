import os, sys, glob
import numpy as n
import matplotlib.pyplot as Mplot
import time
import urllib.request as ur
from bs4 import BeautifulSoup as bs
import json

# mpc raw data
data_name = "MPCORB"
data_exte = ".DAT"

# open the raw data
pyfile_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
find_dat   = os.path.join(f"{pyfile_dir}/**/{data_name}*{data_exte}")
glob_dat   = sorted(glob.glob(find_dat, recursive = True))
# print(glob_dat) # checkpoint

# condition on/off
condition = True

# data conductor
mpc_obj = []

# boundary condition
q_min, q_max = 0.8, 2.0
i_min, i_max = 0.0, 45.0
mpc_H_min = 22.65
mpc_H_max = 22.86
lim_condi_code = 3

# mercury6 setup
coor_mode_list = ("Cartesian", "Asteroidal", "Cometary")
selected_mode = coor_mode_list[1]

def obtain_mpcdata(mpcdat):
    # alert for condition filter:" "on"
    if condition == True:
        print("condition filter: 'on'")
    
    # read raw mpc data
    mpcdata = open(mpcdat, "r")
    mpcr    = mpcdata.readlines()
    for elem in range(43, len(mpcr)):
        # remove empty index
        raw_value = mpcr[elem].split(" ")
        for k in range(len(raw_value)):
            if '' in raw_value:
                raw_value.remove('')
            else:
                continue
        
        # obtain data
        try:
            mpc_id = raw_value[0]
            mpc_name = raw_value[-2]
            mpc_a = float(raw_value[10])
            mpc_e = float(raw_value[8])
            mpc_i = float(raw_value[7]) / (180/n.pi)
            mpc_H = float(raw_value[1])
            mpc_M = float(raw_value[4])
            mpc_o = float(raw_value[6])
            mpc_O = float(raw_value[5])
        except ValueError:
            continue
        except IndexError:
            continue
            
        mpc_q = mpc_a * (1 - mpc_e)        
            
        print(f"loading process... {elem - 41}/{len(mpcr) - 42} ({round(((elem - 41)/ (len(mpcr) - 42)) * 100 , 3)}%)") # checkpoint
        # add condition filter
        if condition == True:
            if q_min < mpc_q < q_max:
                if mpc_H_min < mpc_H < mpc_H_max:
                    mpc_obj.append([mpc_id, mpc_a, mpc_e, mpc_i, mpc_q, mpc_H, mpc_name, mpc_O, mpc_o, mpc_M])
                else:
                    continue
            else:
                continue
        else:
            mpc_obj.append([mpc_id, mpc_a, mpc_e, mpc_i, mpc_q, mpc_H, mpc_name, mpc_O, mpc_o, mpc_M])

    mpcdata.close()
    print(f"filted {len(mpc_obj)} data!")
    return mpc_obj

def plot_qi_distribution():
    # obtain orbital element
    array_mpc_obj = n.array(mpc_obj)
    mpc_id, mpc_a, mpc_e, mpc_i, mpc_q, mpc_H, mpc_name, mpc_O, mpc_o, mpc_M = array_mpc_obj.T
    mpc_a = [float(v) for v in mpc_a]
    mpc_name = [v for v in mpc_name]
    len_name = [len(v) for v in mpc_name]
    mpc_e = [float(v) for v in mpc_e]
    mpc_i = [180 - (float(v) * (180/n.pi)) if float(v) > 90 else float(v) * (180/n.pi) for v in mpc_i]
    #mpc_i = [float(v) for v in mpc_i]
    mpc_q = [float(v) for v in mpc_q]
    mpc_H = [float(v) for v in mpc_H]
    mpc_D = [1329 * 10 ** (-float(v)/5) / n.sqrt(0.14) for v in mpc_H] # set up albedo p_v = 0.14 ref. to Fernandez, 2023
    
    #print(max(mpc_i), "at", list(mpc_id)[list(mpc_i).index(max(mpc_i))])
    
    """
    print("=== LIST OF DATA POINT ===")
    print("id", "name" + "\t" * (int(max(len_name)/8)), "a (au)", "e", "i (deg)", "H", "D (km)", sep = "\t")
    for mpc in range(len(mpc_id)):
       if  len(mpc_name[mpc]) >= 8:
           print(mpc_id[mpc], mpc_name[mpc], round(mpc_a[mpc], 4), round(mpc_e[mpc], 4), round(mpc_i[mpc], 4), mpc_H[mpc], round(mpc_D[mpc], 4), sep = "\t")
       else:
           print(mpc_id[mpc], mpc_name[mpc] + "\t" * (int(max(len_name)/8)), 
                 round(mpc_a[mpc], 4), round(mpc_e[mpc], 4), round(mpc_i[mpc], 4), mpc_H[mpc], round(mpc_D[mpc], 4), sep = "\t")
    print("==========================")
    """
    
    # histogram with i vs number
    figni = Mplot.figure(figsize = (10, 10))
    axni = figni.subplots()
    axni.set_title(f"histogram of minor planets' distribution of inclination between q = {q_min} ~ {q_max} au ({mpc_H_min} < H < {mpc_H_max})")
    axni.set_xlim(0.0, 90.0)
    axni.set_xlabel("i (deg)", fontsize = 12)
    axni.set_ylabel("number", fontsize = 12)
    style = {"facecolor": "white", "edgecolor": "black"}
    axni.hist(mpc_i, 90, **style)
    
    # save figure
    figni.savefig("mca_i-n_distribution.png")
    Mplot.close()
    
    figqi = Mplot.figure(figsize = (10, 10)) 
    Mplot.suptitle(f"q-i distribution of minor planets between q = {q_min} ~ {q_max} au (995-m < D < 1005-m)")
    axqi1 = figqi.add_axes([0.07, 0.07, 0.6, 0.6])
    axqi1.set_xlabel("q (au)")
    axqi1.set_ylabel("i (deg)")
    axqi2 = figqi.add_axes([0.07, 0.72, 0.6, 0.2])
    axqi3 = figqi.add_axes([0.72, 0.07, 0.2, 0.6])
    axqi1.hist2d(mpc_q, mpc_i, bins = 50, cmap="gray_r")
    styleqi = {"facecolor": "white", "edgecolor": "black"}
    axqi1.set_xlim(0.79, 2.01)
    axqi1.set_ylim(i_min, i_max)
    axqi2.set_xlim(0.79, 2.01)
    axqi2.set_ylabel("n")
    axqi3.set_ylim(i_min, i_max)
    axqi3.set_xlabel("n")
    axqi2.hist(mpc_q, 50, **styleqi)
    axqi3.hist(mpc_i, 50, orientation='horizontal', **styleqi)
    axqi1.plot([0.8] * len(mpc_q), n.linspace(0, i_max, len(mpc_q)), linestyle = "--", color = "black", linewidth = 1)
    axqi1.plot([1.3] * len(mpc_q), n.linspace(0, i_max, len(mpc_q)), linestyle = "--", color = "black", linewidth = 1)
    axqi1.plot([1.67] * len(mpc_q), n.linspace(0, i_max, len(mpc_q)), linestyle = "--", color = "black", linewidth = 1)
    axqi1.plot([2.0] * len(mpc_q), n.linspace(0, i_max, len(mpc_q)), linestyle = "--", color = "black", linewidth = 1)
    axqi1.text(0.81, i_max - 1, "NEA", fontsize = 8)
    axqi1.text(1.31, i_max - 1, "MCA and near MCA", fontsize = 8)
    axqi1.text(1.68, i_max - 1, "Main Belt", fontsize = 8)
    Mplot.show()
    
    # save figure
    # Mplot.savefig("mca_q-i_distribution.png")
    Mplot.close()

#record as small.in for mercury6
def rec_smallin(rawfile):
    print("begin writing small.in ...")
    with open(rawfile, "w", encoding = "utf-8") as rf:
        # obtain orbital element
        array_mpc_obj = n.array(mpc_obj)
        mpc_id, mpc_a, mpc_e, mpc_i, mpc_q, mpc_H, mpc_name, mpc_O, mpc_o, mpc_M = array_mpc_obj.T
        mpc_a = [float(v) for v in mpc_a]
        mpc_name = [v for v in mpc_name]
        len_name = [len(v) for v in mpc_name]
        mpc_e = [float(v) for v in mpc_e]
        mpc_i = [180 - (float(v) * (180/n.pi)) if float(v) > 90 else float(v) * (180/n.pi) for v in mpc_i]
        #mpc_i = [float(v) for v in mpc_i]
        mpc_o = [float(v) for v in mpc_o]
        mpc_O = [float(v) for v in mpc_O]
        mpc_M = [float(v) for v in mpc_M]
        
        # add condition code from SBDB
        max_condi_code = lim_condi_code
        condi_code = []
        for datacase in mpc_id:
            try:
                print(f"obtain condition code of mpcid {datacase} ({list(mpc_id).index(datacase)+1})")
                sbdb_api = f"https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={datacase}"
                sbdb_url = ur.urlopen(sbdb_api)
                sbdb_soup = bs(sbdb_url, "html.parser").decode("utf-8")
                condi_obj = float(json.loads(sbdb_soup, strict = False)["orbit"]["condition_code"])
                condi_code.append(condi_obj)
            except KeyError:
                continue

        # record small.in
        rf.write(")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n")
        rf.write(") Lines beginning with ')' are ignored.\n")
        rf.write(")---------------------------------------------------------------------\n")
        rf.write(f" style (Cartesian, Asteroidal, Cometary) = {selected_mode}\n")
        rf.write(")---------------------------------------------------------------------\n")
        
        cc = 0
        for i in range(len(condi_code)):
            if condi_code[i] <= max_condi_code:
                rf.write(str(mpc_id[i]) + " " + f"ep=2460800.5d0\n") # 2025-05-05
                rf.write(str(mpc_a[i]) + "d0 ")
                rf.write(str(mpc_e[i]) + "d0 ")
                rf.write(str(mpc_i[i]) + "d0 ")
                rf.write(str(mpc_O[i]) + "d0 ")
                rf.write(str(mpc_o[i]) + "d0 ")
                rf.write(str(mpc_M[i]) + "d0 ")
                rf.write(("0d0" + " ") * 3)
                rf.write("\n")
                cc += 1
            else:
                continue
    print(f"filtered {cc} data within condition code") # cc < 3
    print(f"finish recording {rawfile}")

def rec_yarkovsky_in(rawfile):
    print("begin writing yarkovsky.in ...")
    with open(rawfile, "w", encoding = "utf-8") as rf:
        # obtain orbital element
        array_mpc_obj = n.array(mpc_obj)
        mpc_id, mpc_a, mpc_e, mpc_i, mpc_q, mpc_H, mpc_name, mpc_O, mpc_o, mpc_M = array_mpc_obj.T
        albedo = 0.14
        mpc_D = [(1329*10**(-0.2 * float(v)))/(albedo**0.5) for v in mpc_H]
        
        # add condition code from SBDB
        max_condi_code = lim_condi_code
        condi_code = []
        for datacase in mpc_id:
            try:
                print(f"obtain condition code of mpcid {datacase} ({list(mpc_id).index(datacase)+1})")
                sbdb_api = f"https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={datacase}"
                sbdb_url = ur.urlopen(sbdb_api)
                sbdb_soup = bs(sbdb_url, "html.parser").decode("utf-8")
                condi_obj = float(json.loads(sbdb_soup, strict = False)["orbit"]["condition_code"])
                condi_code.append(condi_obj)
            except KeyError:
                continue
        
        cc = 0
        #yorp parameter
        mpc_density = 2500
        thermal_conductivity = 0.1
        heat_capacity = 800
        obiliquity = 0
        rotational_period = 4
        
        
        for i in range(len(condi_code)):
            if condi_code[i] <= max_condi_code:
                rf.write(f"{mpc_id[i]} {mpc_density} {thermal_conductivity} {heat_capacity} {mpc_D[i]} {obiliquity} {rotational_period} {1.0} {1.0}\n")
                cc += 1
            else:
                continue
    print(f"filtered {cc} data within condition code") # cc < 3
    print(f"finish recording {rawfile}")
        
"""
main code
"""
obtain_mpcdata(glob_dat[0])
#plot_qi_distribution()
rec_smallin("small.in")
rec_yarkovsky_in("yarkovsky.in")

