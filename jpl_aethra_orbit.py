import astroquery.jplhorizons as aqjpl
import astropy as ap
import numpy as n
import matplotlib.pyplot as Mplot
import datetime, julian
import rebound
import time
import os

# data of mars-crossing asteroid (mca)
def info_mca(mca_id, jd):
    obj = aqjpl.Horizons(id = mca_id, epochs = jd).vectors()
    #print(obj.keys()) # checkpoint
    return obj 

# the mca data's RA and DEC
def mca_rv(mca_id):
    selected_mca = info_mca(mca_id, jd)
    x, y, z = selected_mca["x"], selected_mca["y"], selected_mca["z"]
    vx, vy, vz = selected_mca["vx"], selected_mca["vy"], selected_mca["vz"]
    return [x, y, z, vx, vy, vz]

# object's position and velocity
obj = {"90000091": "Encke", "132": "Aethra", "102528": "1999 US3"}
obj_k, obj_v = list(obj.keys())[-1], list(obj.values())[-1]
jd = 2400000
# 2460676 # 2025-01-01 0:0:0
file_ = open(f"sim_xyzJD{jd}.txt", "w")
file_.close()

print(f"select positional data of {obj_v}")
obj_rv = mca_rv(obj_k)
x_a, y_a, z_a, vx_a, vy_a, vz_a = obj_rv

# mercury's position and velocity
print("select positional data of Mercury")
mercury_rv = mca_rv("199")
x_me, y_me, z_me, vx_me, vy_me, vz_me = mercury_rv

# venus' position and velocity
print("select positional data of Venus")
venus_rv = mca_rv("299")
x_ve, y_ve, z_ve, vx_ve, vy_ve, vz_ve = venus_rv

# earth' position and velocity
print("select positional data of Earth")
earth_rv = mca_rv("3")
x_ea, y_ea, z_ea, vx_ea, vy_ea, vz_ea = earth_rv

# mars' position and velocity
print("select positional data of Mars")
mars_rv   = mca_rv("4")
x_m, y_m, z_m, vx_m, vy_m, vz_m = mars_rv

# jupiter's position and velocity
print("select positional data of Jupiter")
jupiter_rv = mca_rv("5")
x_j, y_j, z_j, vx_j, vy_j, vz_j = jupiter_rv

# saturn's position and velocity
print("select positional data of Saturn")
saturn_rv = mca_rv("6")
x_s, y_s, z_s, vx_s, vy_s, vz_s = saturn_rv

# result
obj_t = info_mca(f"{obj_k}", jd)["datetime_jd"]
jdate = [t for t in obj_t]

# modify the datetime
mjdate = [jdate[i] - jdate[0] for i in range(len(jdate))]

# record position data
resu_dir  = "Orbital_Result"

# plot the distribution
def plot3Dorbit():
    fig = Mplot.figure(figsize = (10, 10))
    ax = fig.add_subplot(projection = "3d")
    ax.scatter(x_a, y_a, z_a, c = "black", s = 7, alpha = 1, label = f"{obj_v}")
    ax.scatter(x_m, y_m, z_m, c = "red", s = 7, alpha = 1, label = "mars")
    ax.scatter(x_j, y_j, z_j, c = "green", s = 7, alpha = 1, label = "jupiter")
    ax.scatter(0, 0, 0, c = "orange", s = 40, label = "sun")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_zlim(-6, 6)
    ax.set_title(f"{obj_v}'s orbit (1994/01/01 ~ 2024/01/01)")
    ax.set_xlabel("$x$ (AU)", fontsize = 12)
    ax.set_ylabel("$y$ (AU)", fontsize = 12)
    ax.set_zlabel("$z$ (AU)", fontsize = 12)
    ax.legend()
    # Mplot.show()
    Mplot.savefig(f"{obj_v}_obrit(3d).png")
    
def plot2Dorbit():
    fig = Mplot.figure(figsize = (10, 10))
    ax = fig.add_subplot()
    ax.scatter(x_a, y_a, c = "black", s = 7, alpha = 1, label = "aethra")
    ax.scatter(x_m, y_m, c = "red", s = 7, alpha = 1, label = "mars")
    ax.scatter(x_j, y_j, c = "green", s = 7, alpha = 1, label = "jupiter")
    ax.scatter(0, 0, c = "orange", s = 40, label = "sun")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_title(f"{obj_v}'s orbit (1994/01/01 ~ 2024/01/01)")
    ax.set_xlabel("$x$ (AU)", fontsize = 12)
    ax.set_ylabel("$y$ (AU)", fontsize = 12)
    ax.legend()
    Mplot.savefig(f"{obj_v}_obrit(2d).png")

def plotzcross():
    fig = Mplot.figure(figsize = (10, 10))
    ax_1, ax_2, ax_3 = fig.subplots(3, 1)
    ax_1.plot(mjdate, [int(0)]*len(aethra_t), "r--", color = "r")
    ax_2.plot(mjdate, [int(0)]*len(aethra_t), "r--", color = "r")
    ax_3.plot(mjdate, [int(0)]*len(aethra_t), "r--", color = "r")
    ax_1.scatter(mjdate, x_a - x_m, c = "black")
    ax_2.scatter(mjdate, y_a - y_m, c = "black")
    ax_3.scatter(mjdate, z_a - z_m, c = "black")
    ax_3.set_xlabel(f"Julian Date - {jdate[0]}", fontsize = 12)
    ax_1.set_ylabel("$x_a - x_m$ (AU)", fontsize = 12)
    ax_2.set_ylabel("$y_a - y_m$ (AU)", fontsize = 12)
    ax_3.set_ylabel("$z_a - z_m$ (AU)", fontsize = 12)
    Mplot.savefig(f"z_{obj_v[0]}m_diff.png")

def plotvdist():
    delta_time = mjdate
    fig = Mplot.figure(figsize = (10, 10))
    ax_1, ax_2, ax_3 = fig.subplots(3, 1)
    ax_1.plot(delta_time, [int(0)]*(len(aethra_t)), "r--", color = "r")
    ax_2.plot(delta_time, [int(0)]*(len(aethra_t)), "r--", color = "r")
    ax_3.plot(delta_time, [int(0)]*(len(aethra_t)), "r--", color = "r")
    ax_1.scatter(delta_time, vx_a, c = "black")
    ax_2.scatter(delta_time, vy_a, c = "black")
    ax_3.scatter(delta_time, vz_a, c = "black")
    ax_3.set_xlabel(f"Julian Date - {jdate[0]}", fontsize = 12)
    ax_1.set_ylabel("$v_{a, x}$ (AU/day)", fontsize = 12)
    ax_2.set_ylabel("$v_{a, y}$ (AU/day)", fontsize = 12)
    ax_3.set_ylabel("$v_{a, z}$ (AU/day)", fontsize = 12)
    Mplot.savefig(f"v_{obj_v[0]}_diff.png")

"""
========================================================================================================================================================
"""

A, E, I, Om, om = [], [], [], [], []
ts_ = [] 
planet_l = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn"]
pcolor = ["yellow", "silver", "gold", "blue", "red", "green", "lime", "black"]

# rebound
def rebound_simulation():
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
                        """
                        sim.add(m = 1.9884754159566474e+30, x = 0, y = 0, z = 0, vx = 0, vy = 0, vz = 0) # Sun
                        sim.add(m = 3.301109449002709e+23, x = x_me, y = y_me, z = z_me, vx = vx_me, vy = vy_me, vz = vz_me) # Mercury
                        sim.add(m = 4.867466257521635e+24, x = x_ve, y = y_ve, z = z_ve, vx = vx_ve, vy = vy_ve, vz = vz_ve) # Venus
                        sim.add(m = 6.045825576341311e+24, x = x_ea, y = y_ea, z = z_ea, vx = vx_ea, vy = vy_ea, vz = vz_ea) # Earth
                        sim.add(m = 6.417120534329241e+23, x = x_m, y = y_m, z = z_m, vx = vx_m, vy = vy_m, vz = vz_m)       # Mars
                        sim.add(m = 1.8985802402728163e+27, x = x_j, y = y_j, z = z_j, vx = vx_j, vy = vy_j, vz = vz_j)      # Jupiter
                        sim.add(m = 5.6847662661820045e+26, x = x_s, y = y_s, z = z_s, vx = vx_s, vy = vy_s, vz = vz_s)      # Saturn 
                        """
                        for p in planet_l:
                            sim.add(p, date = f"JD{jd}")
                        sim.add(m = 0, x = x_a, y = y_a, z = z_a, vx = vx_a, vy = vy_a, vz = vz_a)                           # Object
                            
                        sim.t = 0
                        sim.dt = 7.0
                        #sim.exit_min_distance = 0.007 # Hill radius Jupiter(AU), Elizabeth Lovegrove, 2012: min(rH) = 3.54283e-4 AU
                        max_count = 2000
                                                     
                        ts = [365.24 *500 *t for t in range(max_count)]
                        
                        # objects' position
                        position = {}
                        for p in planet_l:
                            plantary_position = {f"x_{p}": [], f"y_{p}": [], f"z_{p}": []}
                            position[p] = plantary_position
                        position[obj_v] = {f"x_{obj_v}": [], f"y_{obj_v}": [], f"z_{obj_v}": []}
                        
                        distance = {}
                        for p in planet_l:
                            distance[f"r_{p}"] = []
                        
                        t_count = 0    
                        for t in ts:
                            # print(f"dt = {sim.dt}")
                            sim.move_to_com()
                            sim.integrate(t)
                            
                            particles = sim.particles[-1]
                            
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
                            t_count += 1

                            print("init_data processing {}%".format(round(ts.index(t)*100/len(ts), 2)))
                        
                        with open("/home/xi-feng/NCU/Meeting/"+f"sim_xyzJD{jd}.txt", "wt", encoding = "utf-8") as recp:
                            recp.write(str([position, distance]))
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
   axa.scatter(delta_time, A, c = "black", s = 4)
   axe.scatter(delta_time, E, c = "black", s = 4)
   axi.scatter(delta_time, I, c = "black", s = 4)
   axOm.scatter(delta_time, Om, c = "black", s = 4)
   axom.scatter(delta_time, om, c = "black", s = 4)
   Mplot.savefig("para_data.png")
   Mplot.close()
   
   Mplot.figure(figsize = (10, 10))
   q = [A[t] * (1 - E[t]) for t in range(len(A))]
   Mplot.scatter(om, q, c = "black")  
   Mplot.xlabel("ω (deg)") 
   Mplot.ylabel("q (AU)") 
   Mplot.savefig("para_omegaVSq.png")
   Mplot.close() 
   
   Mplot.figure(figsize = (10, 10))
   Mplot.scatter(I, q, c = "black")  
   Mplot.xlabel("i (deg)") 
   Mplot.ylabel("q (AU)") 
   Mplot.savefig("para_incVSq.png")
   Mplot.close()  

def plot_sim_dist(textfile):
   simt = [t / 365.24 for t in ts_]
   dp = []
   fig = Mplot.figure(figsize = (10, 10))
   ax = fig.add_subplot()
   
   with open(textfile, "rt") as recp:
       pos = eval(recp.readlines()[0])[0]
       for p in pos.keys():
           p_x, p_y, p_z = [], [], []
           p_x, p_y, p_z = pos[p][f"x_{p}"], pos[p][f"y_{p}"], pos[p][f"z_{p}"]
           ax.scatter(p_x, p_y, c = pcolor[list(pos.keys()).index(p)], s = 10, alpha = 0.7, label = p)
   ax.set_xlabel("x(AU)")
   ax.set_ylabel("y(AU)")
   #ax.set_zlabel("z(AU)")
   ax.legend()
   Mplot.savefig("planet_distri.png")
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
plot_orbital_element()
plot_sim_dist("/home/xi-feng/NCU/Meeting/"+f"sim_xyzJD{jd}.txt")

    
    
    
    
    
 
