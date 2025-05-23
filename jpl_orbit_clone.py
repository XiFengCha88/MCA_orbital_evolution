import urllib.request as ur
from bs4 import BeautifulSoup as bs
import astroquery.jplhorizons as aqjpl
import json
import numpy as n
import matplotlib.pyplot as Mplot
import rebound
import os, sys, glob
import csv

# data of mars-crossing asteroid (mca)
def info_mca(mca_id, jd):
    obj = aqjpl.Horizons(id = mca_id, epochs = jd).elements()
    #print(obj.keys()) # checkpoint
    return obj 

# the mca data's orbital element
def mca_rv(mca_id):
    selected_mca = info_mca(mca_id, jd)
    a = selected_mca["a"]
    e = selected_mca["e"]
    i = selected_mca["incl"]
    Ω = selected_mca["Omega"]
    ω = selected_mca["w"]
    return [a, e, i/(180/n.pi), Ω/(180/n.pi), ω/(180/n.pi)]

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
    """
    init_e       = float(cov_elem[0]["value"])
    init_a       = float(cov_elem[1]["value"])
    init_i       = float(cov_elem[5]["value"])
    init_Ω       = float(cov_elem[3]["value"])
    init_ω       = float(cov_elem[4]["value"])
    """
    init_a, init_e, init_i, init_Ω, init_ω = mca_rv(mca_id) # aqjpl   
    
    # clone each element
    clone_e, clone_a, clone_i, clone_Ω, clone_ω = [init_e], [init_a], [init_i], [init_Ω], [init_ω]
    for i in range(num):
        gaussian_random = n.matrix(n.random.normal(n.zeros(6), diagonal_element ** 0.5)).reshape(6, 1)
        gaussian_evs    = cm_eigenvectors * gaussian_random
        
        clone_e.append(float(init_e + gaussian_evs[0]))
        clone_a.append(float((init_a + gaussian_evs[1])))
        clone_i.append(float((init_i + gaussian_evs[5]) / (180/n.pi)))
        clone_Ω.append(float((init_Ω + gaussian_evs[3]) / (180/n.pi)))
        clone_ω.append(float((init_ω + gaussian_evs[4]) / (180/n.pi)))
    
    # return the results
    return [clone_a, clone_e, clone_i, clone_Ω, clone_ω]
    
# object's position and velocity
obj = {"90000091": "Encke", "132": "Aethra", "102528": "1999US3", "433": "Eros"}
obj_k, obj_v = list(obj.keys())[-1], list(obj.values())[-1]
pjd = 2457724.5
jd = 2453311.5
""" 
  === EPOCH TIME IN HORIZONS === 
        2P/Encke | jd2459778.5
      132 Aethra | jd2457567.5
  102528 1999US3 | jd2457542.5
        422 Eros | jd2453311.5
"""
# result
# clone_data(obj_k)

# record positional data
resu_dir  = "Orbital_Result"
caseid = "case" + str(5)

"""
========================================================================================================================================================
""" 
# timestep
orbit_period = 643.1403141999031
dt = 0.02 * orbit_period
step = 4000
"""
  === ORBITAL PERIOD (days) ===
  102528 1999US3 | 1680.17298426161
        433 Eros | 643.1403141999031
"""

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
                for ci in range(clone_num):
                    cp_a, cp_e, cp_i, cp_Om, cp_om = [], [], [], [], []
                    sim = rebound.Simulation()
                    sim.units = ("AU", "d", "kg")
                            
                    #initial condition
                    for p in planet_l:
                        sim.add(p, date = f"JD{jd}")
                    sim.add(m = 0, a = clone_a[ci], e = clone_e[ci], inc = clone_i[ci], Omega = clone_Om[ci], omega = clone_om[ci]) # Object
                    
                    for i in range(len(sim.particles)-1):
                        print(planet_l[i], f"x = {sim.particles[i].x}", f"y = {sim.particles[i].y}", f"z = {sim.particles[i].z}", f"at jd  = {jd}", sep = "\t")
                    print("Eros", f"x = {sim.particles[-1].x}", f"y = {sim.particles[-1].y}", f"z = {sim.particles[-1].z}", f"at jd = {jd}", sep = "\t")
 
                    sim.t = 0
                    sim.dt = integrate_dt
                    sim.integrator = "MERCURIUS"
                 
                    max_count = max_number
                        
                    for t in range(max_count):
                        sim.move_to_hel()
                        sim.steps(step)
                            
                        particles = sim.particles[-1]
                            
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
                        
                        print("init_data processing {}/{} ({}/{})".format(t+1, max_number, ci+1, clone_num))
                                
                    A.append(cp_a)
                    E.append(cp_e)
                    I.append(cp_i)
                    Om.append(cp_Om)
                    om.append(cp_om)
                    ts_.append(float(sim.t) / 365.242199)
                            
                    sim = None
                    
                with open(f"/home/xi-feng/NCU/Meeting/clone{clone_num}_sim_orbital_element_jd{jd}_{caseid}.txt", "wt", encoding = "utf-8") as rec_oe:
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
                
                with open(f"/home/xi-feng/NCU/Meeting/clone{clone_num}_sim_position_jd{jd}_{caseid}.csv", "w", encoding = "utf-8") as rec_pos:
                    writer = csv.writer()
                 
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
    
    evolution_time = [dt * step * t / 365.242199 for t in range(max_points)]
    semimajor_axis = rawfile["a"]
    eccentricity   = rawfile["e"]
    inclination    = rawfile["i"]
    longti_node    = rawfile["Ω"]
    arg_perihelion = rawfile["ω"]
    
    fig = Mplot.figure(figsize = (10, 10))
    axa, axe, axOm, axom, axi = fig.subplots(5, 1)
    #axa, axe = fig.subplots(2, 1)
    axa.set_title(f"100 clones of {obj_v}'s orbital evolution")
    
    #axa.set_ylim(0.1, 200)
    axi.set_xlabel("Evolution time (yr)")
    axa.set_ylabel("a (AU)")
    axe.set_ylabel("e")
    axi.set_ylabel("i (deg)")
    axOm.set_ylabel("q (AU)")
    axom.set_ylabel("Q (AU)")
    #axa.set_yscale("log")
    
    perihelion = [[semimajor_axis[cp][i] * (1 - eccentricity[cp][i]) for i in range(len(semimajor_axis[cp]))] for cp in range(0, 100)]
    aphelion   = [[semimajor_axis[cp][i] * (1 + eccentricity[cp][i]) for i in range(len(semimajor_axis[cp]))] for cp in range(0, 100)]
    
    # saved case
    """ 
    === 433 Eros ===
    case1: [17, 22, 25, 37, 57, 65, 68, 72, 80, 93]
    case2: [12, 20, 23, 25, 35, 43, 44, 57, 96, 99]
    case4: [9, 17, 30, 41, 45, 57, 64, 69, 89, 94]
    case5: [1, 2, 9, 29, 34, 37, 61, 74, 78, 88]
    case6: [3, 8, 18, 34, 39, 43, 64, 69, 70, 94]
    """
    datacase = [i for i in range(100)]
    
    for cp in datacase:
        print("case", cp, [min(semimajor_axis[cp]), max(semimajor_axis[cp])])
        axa.plot(evolution_time, semimajor_axis[cp], linewidth = 1, alpha = 0.75)
        axe.plot(evolution_time, eccentricity[cp], linewidth = 1, alpha = 0.75)
        axOm.plot(evolution_time, perihelion[cp], linewidth = 1, alpha = 0.75)
        axom.plot(evolution_time, aphelion[cp], linewidth = 1, alpha = 0.75)
        axi.plot(evolution_time, inclination[cp], linewidth = 1, alpha = 0.75)
        
    axa.plot(evolution_time, semimajor_axis[0], linestyle = "-", color = "black", linewidth = 0.75, alpha = 0.7, label = "original data")
    axe.plot(evolution_time, eccentricity[0], linestyle = "-", color = "black", linewidth = 0.75, alpha = 0.7)
    axi.plot(evolution_time, inclination[0], linestyle = "-", color = "black", linewidth = 0.75, alpha = 0.7)
    axOm.plot(evolution_time, perihelion[0], linestyle = "-", color = "black", linewidth = 0.75, alpha = 0.7)
    axom.plot(evolution_time, aphelion[0], linestyle = "-", color = "black",  linewidth = 0.75, alpha = 0.7)
    axa.legend()
    Mplot.savefig(f"clone_para_a_jd{jd}_{caseid}.png")
    Mplot.close()
    
    #plot a-i relation
    fig_ai = Mplot.figure(figsize = (10, 10))
    ax_ai = fig_ai.subplots()
    ax_ai.set_title(f"a-i relation of {obj_v}'s orbital evolution with 100 clones") 
    ax_ai.set_xlabel("a (AU)")
    ax_ai.set_ylabel("i (deg)")
    ax_ai.set_xlim(0.7, 2)
    ax_ai.set_ylim(0, 25)
    
    for cp in datacase:
        quartor = int(len(semimajor_axis[cp]) / 4)
        ax_ai.scatter([semimajor_axis[cp][i] for i in range(quartor)], [inclination[cp][i] for i in range(quartor)],
                      s = 0.5, c = "blue")
        ax_ai.scatter([semimajor_axis[cp][i] for i in range(quartor, 2 * quartor)], [inclination[cp][i] for i in range(quartor, 2 * quartor)],
                      s = 0.5, c = "red")
        ax_ai.scatter([semimajor_axis[cp][i] for i in range(2 * quartor, 3 * quartor)], [inclination[cp][i] for i in range(2 * quartor, 3 * quartor)],
                      s = 0.5, c = "orange") 
        ax_ai.scatter([semimajor_axis[cp][i] for i in range(3 * quartor, 4 * quartor)], [inclination[cp][i] for i in range(3 * quartor, 4 * quartor)],
                      s = 0.5, c = "green")  
        ax_ai.scatter(semimajor_axis[cp][0], inclination[cp][0], s = 100, c = "black", marker = "*", edgecolor = "black")     
   
    ax_ai.legend(["0 ~ τ/4", "τ/4 ~ τ/2", "τ/2 ~ 3τ/4", "3τ/4 ~ τ", "initial value"]) # τ: simulation time 
    Mplot.savefig(f"paraAvsI_jd{jd}_{caseid}.png")
    Mplot.close()
    
    fig_ae = Mplot.figure(figsize = (10, 10))
    ax_ae = fig_ae.subplots()
    ax_ae.set_title(f"a-e relation of {obj_v}'s orbital evolution with 100 clones") 
    ax_ae.set_xlabel("a (AU)")
    ax_ae.set_ylabel("e")
   
    for cp in datacase:
        quartor = int(len(semimajor_axis[cp]) / 4)
        ax_ae.scatter([semimajor_axis[cp][i] for i in range(quartor)], [eccentricity[cp][i] for i in range(quartor)],
                      s = 0.5, c = "blue")
        ax_ae.scatter([semimajor_axis[cp][i] for i in range(quartor, 2 * quartor)], [eccentricity[cp][i] for i in range(quartor, 2 * quartor)],
                      s = 0.5, c = "red")
        ax_ae.scatter([semimajor_axis[cp][i] for i in range(2 * quartor, 3 * quartor)], [eccentricity[cp][i] for i in range(2 * quartor, 3 * quartor)],
                      s = 0.5, c = "orange") 
        ax_ae.scatter([semimajor_axis[cp][i] for i in range(3 * quartor, 4 * quartor)], [eccentricity[cp][i] for i in range(3 * quartor, 4 * quartor)],
                      s = 0.5, c = "green")     
        ax_ae.scatter(semimajor_axis[cp][0], eccentricity[cp][0], s = 100, c = "black", marker = "*", edgecolor = "black")      
    ax_ae.plot(n.linspace(1.3, 2.1, 1000), [1 - 1.3 / i for i in n.linspace(1.3, 2.1, 1000)], color = "black", linestyle = "--", linewidth = 3)
    ax_ae.plot(n.linspace(1.017, 2.1, 1000), [1 - 1.017 / i for i in n.linspace(1.017, 2.1, 1000)], color = "black", linestyle = "--", linewidth = 3)
    ax_ae.plot(n.linspace(1.58, 2.1, 1000), [1 - 1.58 / i for i in n.linspace(1.58, 2.1, 1000)], color = "black", linestyle = "--", linewidth = 3)
    ax_ae.plot(n.linspace(1.67, 2.1, 1000), [1 - 1.67 / i for i in n.linspace(1.67, 2.1, 1000)], color = "black", linestyle = "--", linewidth = 3)
    ax_ae.plot(n.linspace(0.6, 0.983, 1000), [0.983 / i - 1 for i in n.linspace(0.6, 0.983, 1000)], color = "black", linestyle = "-.", linewidth = 3)
    ax_ae.plot([1.0] * 1000, n.linspace(0.0, 0.64, 1000), color = "black", linestyle = "--", linewidth = 1.5)
    ax_ae.legend(["0 ~ τ/4", "τ/4 ~ τ/2", "τ/2 ~ 3τ/4", "3τ/4 ~ τ", "initial value"]) # τ: simulation time 
    ax_ae.text(2.05, 1 - 1.67 / 2.1 + 0.01, "q$_{1.67}$", {"fontsize": 14})
    ax_ae.text(2.05, 1 - 1.3 / 2.1 + 0.01, "q$_{1.3}$", {"fontsize": 14})
    ax_ae.text(2.05, 1 - 1.58 / 2.1 + 0.01, "q$_{1.58}$", {"fontsize": 14})
    ax_ae.text(2.05, 1 - 1.017 / 2.1 + 0.01, "q$_{1.017}$", {"fontsize": 14})
    ax_ae.text(0.6, 0.983 / 0.6 - 1 + 0.01, "Q$_{0.983}$", {"fontsize": 14})
    Mplot.savefig(f"paraAvsE_jd{jd}_{caseid}.png")
    Mplot.close()
    
# classification  
def classify_element(dirname = ""):
    clone_count = 1
    
    rawfile_path = os.path.realpath(sys.argv[0])
    rawfile_dir  = os.path.dirname(rawfile_path)
    data_class   = sorted(glob.glob(os.path.join(rawfile_dir + f"/{dirname}", "**", f"clone100*.txt"), recursive = True)) 
    
    # collect clone data
    data_t, data_a, data_e, data_q, data_Q = [], [], [], [], []
    for rawfile in data_class:
        if str(jd) in rawfile:  
            print(f"open {rawfile}...")
            with open(rawfile, "r") as rf:
                rawelem = eval(rf.readlines()[0])
                data_t = [dt * step * t / 365.242199 for t in range(max_points)]
                data_a += rawelem["element"]["a"]
                data_e += rawelem["element"]["e"]
        else:
            continue
    data_q += [[data_a[a][i] * (1 - data_e[a][i]) for i in range(len(data_a[a]))] for a in range(len(data_a))]
    data_Q += [[data_a[a][i] * (1 + data_e[a][i]) for i in range(len(data_a[a]))] for a in range(len(data_a))]
    
    clone_count = len(data_a) - (len(data_class) - 2)
    print("total data:", clone_count)
    
    # classification counting
    print("start classifying the data")
    q_minus, q_balance, q_plus = 0, -4, 0
    
    for cp in range(len(data_a)):
        asteroid_class_ = []
        print(f"classify case {cp + 1}/{len(data_a)}")
        for t in data_t:
            index_t = data_t.index(t)
            if data_a[cp][index_t] > 1.0:
                if 1.017 < data_q[cp][index_t] < 1.3:
                    #print(f"become Amors at t = {t} yr (data index {index_t})")
                    asteroid_class_.append("Amor")
                elif data_q[cp][index_t] < 1.017:
                    #print(f"become Apollo at t = {t} yr (data index {index_t})")
                    asteroid_class_.append("Apollo")
                elif 1.3 < data_q[cp][index_t] < 1.58:
                    #print(f"become deep mars-crosser at t = {} yr (data index {index_t}")
                    asteroid_class_.append("dmca")
                elif 1.58 < data_q[cp][index_t] < 1.67:
                    #print(f"become shallow mars-crosser at t = {} yr (data index {index_t}")
                    asteroid_class_.append("smca")
                else:
                    continue
            elif data_a[cp][index_t] < 1.0:
                if data_Q[cp][index_t] > 0.983:
                    #print(f"become Atens at t = {t} yr (data index {index_t})")
                    asteroid_class_.append("Atens")
                elif data_Q[cp][index_t] < 0.983:
                    #print(f"become Atiras at t = {t} yr (data index {index_t})")
                    asteroid_class_.append("Atiras")
                else:
                    continue
            else:
                continue
       
        #counting
        if asteroid_class_ == [asteroid_class_[0]] * len(data_t):
            q_balance += 1
        else: 
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
    
# main operator
#rebound_simulation(dt, step, max_points)
plot_orbital_element(f"/home/xi-feng/NCU/Meeting/clone100_sim_orbital_element_jd{jd}_{caseid}.txt")
#classify_element()



