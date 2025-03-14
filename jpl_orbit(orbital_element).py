import urllib.request as ur
from bs4 import BeautifulSoup as bs
import json
import numpy as n
import matplotlib.pyplot as Mplot
import datetime, julian
import rebound
import time
import os

# data of mars-crossing asteroid (mca)
def info_mca(mca_id):
    url = f"https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={mca_id}&full-prec=True"
    data_url = ur.urlopen(url)
    data_soap = bs(data_url, "html.parser").decode("utf-8")
    obj = json.loads(data_soap, strict=False)
    #print(obj.keys()) # checkpoint
    return obj 

# the mca data's orbital_element
def mca_rv(mca_id):
    selected_mca = info_mca(mca_id)
    elem = selected_mca["orbit"]["elements"]
    a = float(elem[1]["value"]) # semimajor axis
    e = float(elem[0]["value"]) # eccentricity
    i = float(elem[3]["value"]) # inclination
    Ω = float(elem[4]["value"]) # longitude of the ascending node
    ω = float(elem[5]["value"]) # argument of perihelion
    epoch = float(selected_mca["orbit"]["epoch"])
    return [a, e, i, Ω, ω, epoch]

# object's position and velocity
obj = {"90000091": "Encke", "132": "Aethra", "102528": "1999US3"}
obj_k, obj_v = list(obj.keys())[-1], list(obj.values())[-1]

# result
o_a, o_e, o_i, o_Om, o_om, jd = mca_rv(obj_v)

# record positional data
resu_dir  = "Orbital_Result"

"""
========================================================================================================================================================
"""

A, E, I, Om, om = [], [], [], [], []
ts_ = [] 
planet_l = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn"]
pcolor = ["yellow", "silver", "gold", "blue", "red", "green", "lime", "black"]

# rebound
def rebound_simulation():
    create_file = open(f"sim_xyzjd{jd}.txt", "wt", encoding = "utf-8")
    create_file.close()
    # add particles
    while True:
        try:
            select_existed_particle = input("Do you want to use the existed data for initial condition? [Y/N] Ans: ")
            if select_existed_particle == "Y" or "y":
                print("You choose [Y]")
                os.chdir(resu_dir)
                
                while True: 
                    try:
                               
                        # begin with the simulation
                        sim = rebound.Simulation()
                        sim.units = ("AU", "d", "kg")
                            
                        #initial condition
                        for p in planet_l:
                            sim.add(p, date = f"JD{jd}")
                        sim.add(m = 0, a = o_a, e = o_e, inc = o_i, Omega = o_Om, omega = o_om) # Object
                            
                        sim.t = 0
                        sim.dt = 7.0
                        max_count = 1000000
                        
                        #timescale                             
                        ts = [365.24 *1 *t for t in range(max_count)]
                        
                        for t in ts:
                            # print(f"dt = {sim.dt}")
                            sim.move_to_com()
                            sim.integrate(t)
                            
                            particles = sim.particles[-1]
                            
                            a = round(float(particles.a), 4)
                            e = round(float(particles.e), 4)
                            i = round(float(particles.inc)* 180 / n.pi, 4)
                            Ω = round(float(particles.Omega)* 180 / n.pi, 4) 
                            ω = round(float(particles.omega)* 180 / n.pi, 4)
                                
                            A.append(a)
                            E.append(e)
                            I.append(i)
                            Om.append(Ω)
                            om.append(ω)
                            ts_.append(t)

                            print("init_data processing {}%".format(round(ts.index(t)*100/len(ts), 2)))
                        sim.stop()
                        break
                    except Exception as err:
                        print(err.str())
                        break

            elif existed_particle == "N" or "n":
                print("You choose [N]")
                xyzvkeys = []
            break
        except Exception as err:
            print(err.str())
            break

# Plot orbital element with time domain
def plot_orbital_element():
   fig = Mplot.figure(figsize = (10, 10))
   axa, axe, axi, axOm, axom = fig.subplots(5, 1)
   axom.set_xlabel(f"Evolution time (yr)")
   axa.set_ylabel("a (AU)")
   axe.set_ylabel("e")
   axi.set_ylabel("i (deg)")
   axOm.set_ylabel("Ω (deg)")
   axom.set_ylabel("ω (deg)")
                          
   delta_time = [ts_[t] / 365.24 for t in range(len(ts_))]
   axa.scatter(delta_time, A, c = "black", s = 1)
   axe.scatter(delta_time, E, c = "black", s = 1)
   axi.scatter(delta_time, I, c = "black", s = 1)
   axOm.scatter(delta_time, Om, c = "black", s = 1)
   axom.scatter(delta_time, om, c = "black", s = 1)
   Mplot.savefig(f"para_data(jd{jd}).png")
   Mplot.close()
   
   Mplot.figure(figsize = (10, 10))
   q = [A[t] * (1 - E[t]) for t in range(len(A))]
   Mplot.scatter(om, q, c = "black")  
   Mplot.xlabel("ω (deg)") 
   Mplot.ylabel("q (AU)") 
   Mplot.savefig(f"para_omegaVSq(jd{jd}).png")
   Mplot.close() 
   
   Mplot.figure(figsize = (10, 10))
   Mplot.scatter(I, q, c = "black")  
   Mplot.xlabel("i (deg)") 
   Mplot.ylabel("q (AU)") 
   Mplot.savefig(f"para_incVSq(jd{jd}).png")
   Mplot.close()  

"""
========================================================================================================================================================
"""
# main operator

rebound_simulation()
plot_orbital_element()

