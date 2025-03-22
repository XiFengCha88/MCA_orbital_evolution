import astroquery.jplhorizons as aqjpl
import numpy as n
import matplotlib.pyplot as Mplot
import rebound
import os

# data of mars-crossing asteroid (mca)
def info_mca(mca_id, jd):
    obj = aqjpl.Horizons(id = mca_id, epochs = jd).elements()
    #print(obj.keys()) # checkpoint
    return obj 

# the mca data's orbital element
def mca_rv(mca_id):
    selected_mca = info_mca(mca_id, jd)
    #x, y, z = selected_mca["x"], selected_mca["y"], selected_mca["z"]
    #vx, vy, vz = selected_mca["vx"], selected_mca["vy"], selected_mca["vz"]
    #return [x, y, z, vx, vy, vz]
    a = selected_mca["a"]
    e = selected_mca["e"]
    i = selected_mca["incl"]
    Ω = selected_mca["Omega"]
    ω = selected_mca["w"]
    return [a, e, i/(180/n.pi), Ω/(180/n.pi), ω/(180/n.pi)]

# object's position and velocity
obj = {"90000091": "Encke", "132": "Aethra", "102528": "1999 US3"}
obj_k, obj_v = list(obj.keys())[-1], list(obj.values())[-1]
jd = 2457542.5 # initial epoch time

# object's orbital element
o_a, o_e, o_i, o_Om, o_om = mca_rv(obj_v)

# result
obj_t = info_mca(f"{obj_k}", jd)["datetime_jd"]
jdate = [t for t in obj_t]

# modify the datetime
mjdate = [jdate[i] - jdate[0] for i in range(len(jdate))]

# record position data
resu_dir  = "Orbital_Result"

"""
========================================================================================================================================================
"""

# data collection
A, E, I, Om, om = [], [], [], [], []
ts_ = [] 
planet_l = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
pcolor = ["yellow", "silver", "gold", "blue", "red", "green", "lime", "skyblue", "violet", "black"]

# time step
dt = 7.0
max_count = 5000000
step = 1

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
                
                # begin the simulation
                sim = rebound.Simulation()
                sim.units = ("AU", "d", "kg")
                        
                # use integrator in this simulation
                sim.integrator = "MERCURIUS"
                        
                # add particle
                for p in planet_l:
                    sim.add(p, date = f"JD{jd}")
                sim.add(m = 0, a = o_a, e = o_e, inc = o_i, Omega = o_Om, omega = o_om)                           # Object
                        
                # set up initial time and time step   
                sim.t = 0
                sim.dt = dt
                        
                # objects' position
                position = {}
                for p in planet_l:
                    plantary_position = {f"x_{p}": [], f"y_{p}": [], f"z_{p}": []}
                    position[p] = plantary_position
                position[obj_v] = {f"x_{obj_v}": [], f"y_{obj_v}": [], f"z_{obj_v}": []}
                        
                distance = {}
                for p in planet_l:
                    distance[f"r_{p}"] = []
   
                for t in range(max_count):
                    # print(f"dt = {sim.dt}")
                    sim.move_to_hel()
                    sim.steps(step)
                            
                    particles = sim.particles[-1]
                    #print(particles.orbit())
                    
                    # add distance from the object to each planet        
                    for p in planet_l:
                        pos_planet = position[p]
                        pos_planet[f"x_{p}"].append(round(sim.particles[planet_l.index(p)].x, 4))
                        pos_planet[f"y_{p}"].append(round(sim.particles[planet_l.index(p)].y, 4))
                        pos_planet[f"z_{p}"].append(round(sim.particles[planet_l.index(p)].z, 4))
                                
                        opx = abs((position[p][f"x_{p}"][-1] - particles.x)) ** 2 
                        opy = abs((position[p][f"y_{p}"][-1] - particles.y)) ** 2
                        opz = abs((position[p][f"z_{p}"][-1] - particles.z)) ** 2
                        
                        distop = round((opx + opy + opz) ** 0.5, 4)
                        distance[f"r_{p}"].append(distop)
                            
                        position[obj_v][f"x_{obj_v}"].append(round(particles.x, 4))
                        position[obj_v][f"y_{obj_v}"].append(round(particles.y, 4))
                        position[obj_v][f"z_{obj_v}"].append(round(particles.z, 4))
                            
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


                        print("init_data processing {}%".format(round(t*100/max_count, 2)), sim.dt)
                        
                    with open("/home/xi-feng/NCU/Meeting/"+f"sim_xyzJD{jd}.txt", "wt", encoding = "utf-8") as recp:
                        recp.write(str([position, distance]))
                        
                    with open(f"/home/xi-feng/NCU/Meeting/sim_orbital_element_jd{jd}.txt", "wt", encoding = "utf-8") as rec_oe:
                        orbital_element = {}
                        orbital_element["name"] = obj_v
                        orbital_element["element"] = {"t": [],"a": [], "e": [], "i": [], "Ω": [], "ω": []}
                        orbital_element["element"]["t"] = ts_
                        orbital_element["element"]["a"] = A
                        orbital_element["element"]["e"] = E
                        orbital_element["element"]["i"] = I
                        orbital_element["element"]["Ω"] = Om
                        orbital_element["element"]["ω"] = om
                        rec_oe.write(str(orbital_element))
                        
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
def plot_orbital_element(textfile):
   orbital_file = open(textfile, "rt")
   readfile = eval(orbital_file.readline())

   ts_ = readfile["element"]["t"]
   A = readfile["element"]["a"]
   E = readfile["element"]["e"]
   I = readfile["element"]["i"]
   Om = readfile["element"]["Ω"]
   om = readfile["element"]["ω"]

   fig = Mplot.figure(figsize = (40, 10))
   axa, axe, axi, axOm, axom = fig.subplots(5, 1)
   axom.set_xlabel(f"Evolution time (yr)")
   axa.set_ylabel("a (AU)")
   axe.set_ylabel("e")
   axi.set_ylabel("i (deg)")
   axOm.set_ylabel("Ω (deg)")
   axom.set_ylabel("ω (deg)")
                          
   delta_time = [t * step * dt / 365.24 for t in range(max_count)]
   axa.scatter(delta_time, A, c = "black", s = 0.25)
   axe.scatter(delta_time, E, c = "black", s = 0.25)
   axi.scatter(delta_time, I, c = "black", s = 0.25)
   axOm.scatter(delta_time, Om, c = "black", s = 0.25)
   axom.scatter(delta_time, om, c = "black", s = 0.25)
   Mplot.savefig(f"para_data(jd{jd}).png")
   Mplot.close()
   
   Mplot.figure(figsize = (10, 10))
   q = [A[t] * (1 - E[t]) for t in range(len(A))]
   Mplot.scatter(om, q, c = "black", s = 0.25)  
   Mplot.xlabel("ω (deg)") 
   Mplot.ylabel("q (AU)") 
   Mplot.savefig(f"para_omegaVSq(jd{jd}).png")
   Mplot.close() 
   
   Mplot.figure(figsize = (10, 10))
   Mplot.scatter(I, q, c = "black", s = 0.25)  
   Mplot.xlabel("i (deg)") 
   Mplot.ylabel("q (AU)") 
   Mplot.savefig(f"para_incVSq(jd{jd}).png")
   Mplot.close()  
   
   orbital_file.close()

def plot_sim_dist(textfile):
   if "data_frame" not in os.listdir():
       os.mkdir("data_frame")
   os.chdir("data_frame")
   
   simt = [t / 365.24 for t in ts_]
   dp = []
   
   pos = {}
   with open(textfile, "rt") as recp:
       pos = eval(recp.readlines()[0])[0]
   
   for i in range(0, max_count, 10):
       fig = Mplot.figure(figsize = (10, 10))
       ax = fig.add_subplot()
       ax.set_xlim(-10, 10)
       ax.set_ylim(-10, 10)
       ax.set_xlabel("x(AU)")
       ax.set_ylabel("y(AU)")
       #ax.set_zlabel("z(AU)")
       for p in range(len(planet_l)):
           ax.scatter(pos[planet_l[p]][f"x_{planet_l[p]}"][i], pos[planet_l[p]][f"y_{planet_l[p]}"][i], s = 75, c = pcolor[p], label = planet_l[p])
       ax.scatter(pos[obj_v][f"x_{obj_v}"][i], pos[obj_v][f"y_{obj_v}"][i], s = 75, c = pcolor[-1], label = obj_v)
       ax.legend()
       Mplot.savefig(f"para_distri_JD{jd}_{str(i).zfill(7)}.png")
       print("finish loading index", i)
       Mplot.close()

   """
   fig2 = Mplot.figure(figsize = (80, 10))
   ax = fig2.add_subplot()
   
   with open(textfile, "rt") as recp:
      dop = eval(recp.readlines()[0])[1]
      for d in dop.keys():
          dp = dop[d]
          ax.plot(simt, dp, pcolor[list(dop.keys()).index(d)], linewidth = 1, label = d)
   ax.set_xlabel("evolution time (yr)")
   ax.set_ylabel("distance (AU)")
   ax.legend()
   Mplot.savefig("planet_distance.png") 
   Mplot.close()           
   """
"""
========================================================================================================================================================
"""
# main operator

rebound_simulation()
plot_orbital_element(f"/home/xi-feng/NCU/Meeting/sim_orbital_element_jd{jd}.txt")
#plot_sim_dist("/home/xi-feng/NCU/Meeting/"+f"sim_xyzJD{jd}.txt")

    
    
    
    
    
 
