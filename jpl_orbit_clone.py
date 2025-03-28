import urllib.request as ur
from bs4 import BeautifulSoup as bs
import json
import numpy as n
import matplotlib.pyplot as Mplot
import rebound
import os

# object's data collection
def obj_collection(obj_id):
    url       = f"https://ssd-api.jpl.nasa.gov/sbdb.api?sstr={obj_id}&full-prec=True&cov=mat"
    data_url  = ur.urlopen(url)
    data_soap = bs(data_url, "html.parser").decode("utf-8")
    obj       = json.loads(data_soap, strict=False)
    #print(obj.keys()) # checkpoint
    return obj 

# object's data 
def data_covariant_matrix(mca_id):
    selected_obj = obj_collection(mca_id)
    cov_array    = selected_obj["orbit"]["covariance"]["data"]
    return n.array([i for i in cov_array]).astype(float)

# clone data in number range
def clone_data(mca_id, num = 100):
    # obtain the covariance matrix
    selected_cov_matrix = data_covariant_matrix(mca_id)
    
    # obtain the eigenvectors from every eigenvalue
    cm_eigenvectors = n.matrix(n.linalg.eig(selected_cov_matrix)[1])
    
    # digonalize covarience matrix: (D^-1)(CM)(D)
    diagonalized_cm = cm_eigenvectors.I * selected_cov_matrix * cm_eigenvectors
    
    # give zero if the element in diagonalized matrix is negative 
    for i in range(len(diagonalized_cm)):
        if diagonalized_cm[i, i] < 0:
            diagonalized_cm[i, i] = 0
        else:
            continue
            
    # take the diagonal elements
    diagonal_element = n.diag(diagonalized_cm)
    
    # obtain JPL's Horizons elements
    selected_obj = obj_collection(mca_id)
    cov_elem     = selected_obj["orbit"]["covariance"]["elements"]
    
    init_e       = float(cov_elem[0]["value"])
    init_a       = float(cov_elem[1]["value"])
    init_i       = float(cov_elem[5]["value"])
    init_Ω       = float(cov_elem[3]["value"])
    init_ω       = float(cov_elem[4]["value"])
    
    # clone each element
    clone_e, clone_a, clone_i, clone_Ω, clone_ω = [init_e], [init_a/(1-init_e)], [init_i], [init_Ω], [init_ω]
    for i in range(num):
        gaussian_random = n.matrix(n.random.normal(n.zeros(6), diagonal_element ** 0.5)).reshape(6, 1)
        gaussian_evs    = cm_eigenvectors * gaussian_random
        clone_e.append(float(init_e + gaussian_evs[0]))
        clone_a.append(float((init_a + gaussian_evs[1]) / (1 - clone_e[-1])))
        clone_i.append(float(init_i + gaussian_evs[5]))
        clone_Ω.append(float(init_Ω + gaussian_evs[3]))
        clone_ω.append(float(init_ω + gaussian_evs[4]))
    
    # return the result
    return [clone_a, clone_e, clone_i, clone_Ω, clone_ω]
    
# object's position and velocity
obj = {"90000091": "Encke", "132": "Aethra", "102528": "1999US3"}
obj_k, obj_v = list(obj.keys())[-1], list(obj.values())[-1]
jd = 2457542.5

# result
# clone_data(obj_k)

# record positional data
resu_dir  = "Orbital_Result"

"""
========================================================================================================================================================
""" 
# timestep
dt = 16.8017298426161
step = 1000

# data count
max_points = 10000

planet_l = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
pcolor = ["yellow", "silver", "gold", "blue", "red", "green", "lime", "skyblue", "violet", "black"]


# rebound
def rebound_simulation(integrate_dt, step, max_number, clone_num = 100):
    A, E, I, Om, om = [], [], [], [], []
    ts_ = []
    clone_a, clone_e, clone_i, clone_Om, clone_om = clone_data(obj_k)
    
    # add particles
    while True:
        try:
            select_existed_particle = input("Do you want to use the existed data for initial condition? [Y/N] Ans: ")
            if select_existed_particle == "Y" or "y":
                print("You choose [Y]")
                os.chdir(resu_dir)
     
                # begin with the simulation
                sim = rebound.Simulation()
                sim.units = ("AU", "d", "kg")
                            
                #initial condition
                for p in planet_l:
                    sim.add(p, date = f"JD{jd}")
                for i in range(0, clone_num):
                    sim.add(m = 0, a = clone_a[i], e = clone_e[i], inc = clone_i[i], Omega = clone_Om[i], omega = clone_om[i]) # Object
                            
                sim.t = 0
                sim.dt = integrate_dt
                max_count = max_number
                        
                for t in range(max_number):
                    cp_a, cp_e, cp_i, cp_Om, cp_om = [], [], [], [], []
                    sim.move_to_com()
                    sim.steps(step)
                       
                    for cn in range(9, len(sim.particles)):     
                        particles = sim.particles[cn]
                            
                        a = round(float(particles.a), 4)
                        e = round(float(particles.e), 4)
                        i = round(float(particles.inc)* 180 / n.pi, 4)
                        Ω = round(float(particles.Omega)* 180 / n.pi, 4) 
                        ω = round(float(particles.omega)* 180 / n.pi, 4)
                        
                        cp_a.append(a)
                        cp_e.append(e)
                        cp_i.append(i)
                        cp_Om.append(Ω)
                        cp_om.append(ω)
                                
                    A.append(cp_a)
                    E.append(cp_e)
                    I.append(cp_i)
                    Om.append(cp_Om)
                    om.append(cp_om)
                    ts_.append(sim.t)

                    print("init_data processing {}/{}".format(t, max_number)) #round(ts.index(t)*100/len(ts), 2)
                            
                sim.stop()
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
                    # with open(f"/home/xi-feng/NCU/Meeting/sim_object_postition_jd{jd}.txt", "wt", encoding = "utf-8") as rec_oe:
                    
            elif existed_particle == "N" or "n":
                print("You choose [N]")
                xyzvkeys = []
            break
        except Exception as err:
            print(err.str())
            break

# Plot orbital element with time domain
def plot_orbital_element(textfile):
    rawfile = {}
    with open(textfile, "rt") as recfile:
        rawdata = eval(recfile.readlines()[0])
        rawfile = rawdata["element"]
    
    evolution_time = rawfile["t"]
    semimajor_axis = n.array(rawfile["a"]).T.tolist()
    #eccentricity   = rawfile["e"]
    #inclination    = rawfile["i"]
    #longti_node    = rawfile["Ω"]
    #arg_perihelion = rawfile["ω"]
    
    fig = Mplot.figure(figsize = (20, 10))
    #axa, axe, axi, axOm, axom = fig.subplots(5, 1)
    axa = fig.subplots() 
    
    axa.set_xlabel("Evolution time (yr)")
    axa.set_ylabel("Δa (AU)")
    """
    axe.set_ylabel("e")
    axi.set_ylabel("i (deg)")
    axOm.set_ylabel("Ω (deg)")
    axom.set_ylabel("ω (deg)")
    """
    for cp in range(1, 100):
        print(len(evolution_time), len(semimajor_axis[cp]))
        axa.scatter(evolution_time, [semimajor_axis[cp][i] - semimajor_axis[0][i] for i in range(len(semimajor_axis[cp]))], s = 1, c = "black")
    axa.scatter(evolution_time, [semimajor_axis[0][i] - semimajor_axis[0][i] for i in range(len(semimajor_axis[0]))], s = 1, c = "yellow")
    
    Mplot.savefig(f"clone_para_a_jd{jd}.png")
# main operator
#rebound_simulation(dt, step, max_points)
plot_orbital_element(f"/home/xi-feng/NCU/Meeting/sim_orbital_element_jd{jd}.txt")




