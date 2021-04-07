import math
import time
import matplotlib.pyplot as plt

#simulator mode, set false for specific scenario (vars at bottom)
perform_cotters = False

# particle constants
p_particle_m = [1000, 1500] # kg / m^3 
d_particle_m = [2.5E-6, 10.0E-6] # m
q_particle_m = [2.0E-17, 2.0E-17] # C
c_particles_m = [29.7E-9 , 44.9E-9] # concentration kg/m^3
h_particle_init_m = [0.1, 0.5]#[1.5E-2, 1.5E-2] # initial height above plate
v_particle_init = 0

# fluid constants
p_fluid_m = [1.293,1.293] # kg / m^3
u_fluid_m = [1.81E-5,1.81E-5] # Ns / m^2 p*v
v_fluid_m = [0.1, 5] # m/s

# plate constants
u_plate_m = [200,2000] # Volts (potential applied to plate)

# cotter's settings
x_terms = [p_particle_m, d_particle_m, q_particle_m, c_particles_m, \
           h_particle_init_m, p_fluid_m, u_fluid_m, u_plate_m, v_fluid_m]
x_terms_names = ['p_particle', 'd_particle', 'q_particle', 'c_particles', \
                'h_particle_init', 'p_fluid', 'u_fluid', 'u_plate', 'v_fluid_m']

# general constants
g = 9.81 # m/s^2
e0 = 8.854E-12 # F*m^-1

dt = 0.00002

# max height to look at surrounding charges
def hSurround(h_particle_init):
    return 2 * h_particle_init

def mParticle(p_particle, d_particle):
    return math.pi / 6 * p_particle * math.pow(d_particle, 3)

# Gravitational Force (const) (downward) --------------------------
def fGravity(p_particle, d_particle):
    return - math.pi / 6.0 * p_particle * g * math.pow(d_particle, 3)

# Buoyant Force (cosnt) (upward) --------------------------
def fBuoyancy(p_fluid, d_particle):
    return math.pi / 6.0 * p_fluid * g * math.pow(d_particle, 3)

# Drag Force (depends on v apparent) (upward) --------------------------
def Cd(v_apparent, d_particle):
    #print(p_fluid * d_particle * math.fabs(v_apparent) / u_fluid)
    if(v_apparent != 0):
        return 24 / (1.39E-5 / v_apparent / d_particle)
    else:
        return 0
    #return 24 / (p_fluid * d_particle * math.fabs(v_apparent) / u_fluid)

def fDrag(v_apparent, u_fluid, d_particle):
    if(v_apparent != 0):
        return 0.5 * Cd(math.fabs(v_apparent), d_particle) * 1.293 * 0.25 * math.pi * math.pow(d_particle, 2) * \
            math.pow(v_apparent, 2)
    else:
        return 0
    #return 6.0 * math.pi * u_fluid * d_particle / 2 * v_apparent

# Electric Field Force due to plate (const) (downward) --------------------------
def fPlate(d, q_particle, u_plate):
    return - q_particle * u_plate / d # E = V / d , F = Eq
    #return  - q_particle * u_plate / 4.0 / math.pi / e0

# Electric field force due to surrounding charges (depends on height) --------------------------
def fSurrCharges(d, q_particle, p_particle, d_particle, c_particles, H):
    return math.pow(q_particle, 2) * c_particles / mParticle(p_particle, d_particle) / 2.0 / e0 * (2*d - H)

def totalForce(d, v, p_particle, d_particle, q_particle, c_particles, h_particle_init, 
                     p_fluid, u_fluid, u_plate):
     #print("g: {:.3e} b: {:.3e} d: {:.3e} plate: {:.3e} surr: {:.3e}".format(\
     #    fGravity(p_particle,d_particle), fBuoyancy(p_fluid,d_particle), fDrag(v,u_fluid,d_particle), \
     #        fPlate(d, q_particle, u_plate), fSurrCharges(d, q_particle, p_particle, d_particle, c_particles, hSurround(h_particle_init))))
     #return fGravity(p_particle, d_particle) + fBuoyancy(p_fluid, d_particle) + \
     return fBuoyancy(p_fluid, d_particle) + \
            fDrag(v, u_fluid, d_particle) + fPlate(d, q_particle, u_plate) + \
            fSurrCharges(d, q_particle, p_particle, d_particle, c_particles, hSurround(h_particle_init))

# model will start with initial v and d
# each time, calculate new v and d based on force calculation
def estimateTravelTime(p_particle, d_particle, q_particle, c_particles, 
                       h_particle_init, p_fluid, u_fluid, u_plate):

    current_h = h_particle_init
    current_v = 0
    current_a = 0
    current_t = 0

    # run model until the current h becomes <= 0 (hits the plate)

    counter = 0
    max_count = 100000000
    while(current_h > 0 and counter < max_count):
        current_a = totalForce(current_h, current_v, p_particle, d_particle, 
                               q_particle, c_particles, h_particle_init, 
                               p_fluid, u_fluid, u_plate) / mParticle(p_particle, d_particle)
        # dv = a*dt
        current_v += current_a * dt
        # dh = v*dt
        current_h += current_v * dt
        current_t += dt

        #print("T: {:5.3f} H: {:5.4e} V: {:5.4e} A: {:5.4e}: ".format(current_t, 
        #    current_h, current_v, current_a))

        # history_h.append(current_h)
        # history_v.append(current_v)
        # history_a.append(current_a)
        # history_t.append(current_t)
        counter += 1
    
    if(counter >= max_count):
        print("Counter maxed out!")
    print("Run time: {:.6f} End Velocity: {:.6f} End Accel: {:.6f}".format(current_t, current_v, current_a))
    return current_t

def estimateMinLength(p_particle, d_particle, q_particle, c_particles, 
                       h_particle_init, p_fluid, u_fluid, u_plate, v_fluid):
    t = estimateTravelTime(p_particle, d_particle, q_particle, c_particles, 
                       h_particle_init, p_fluid, u_fluid, u_plate)
    #print("Min Length: {:.2f}".format(v_fluid * t))
    return v_fluid * t

# Cotter's Method Functions
# note, j is 1 based indexing (that is how Cotter's method is defined)
def codd(j, y):
    j = int(j)
    n = int((len(y) - 2) / 2)
    return 0.25 * ( (y[2*n+2-1] - y[j+n+1-1]) + (y[j+1-1] - y[1-1]))
def ceven(j, y):
    j = int(j)
    n = int((len(y) - 2) / 2)
    return 0.25 * ( (y[2*n+2-1] - y[j+n+1-1]) - (y[j+1-1] - y[1-1]))
def m(j, y):
    return math.fabs(codd(j, y)) + math.fabs(ceven(j, y))
def s(j, y):
    m_list = []
    for i in range(int((len(y)-2) / 2)):
        m_list.append(m(i+1, y))
    sum_m = sum(m_list)
    return m(j, y) / sum_m
def isSignificant(j, y):
    return s(j, y) >= 1.0/len(y)

# x terms contains list of each variable in which min is index 0 and max is index 1
# remember cotter's method stuff is 1 based indexing
y_terms = []

# L H L L L H H H 
# L L H L H L H H < for 3 x's
# L L L H H H L H
if(perform_cotters):
    low_terms, high_terms = list(zip(*x_terms))
    low_terms = list(low_terms)
    high_terms = list(high_terms)

    # Low Case
    y_terms.append(estimateMinLength(*low_terms))

    # Iterate through with all lows and one high
    for i in range(len(low_terms)):
        input_list = low_terms.copy()
        # Make one high
        input_list[i] = high_terms[i] 
        y_terms.append(estimateMinLength(*input_list))

    # Iterate through with all highs and one low
    for i in range(len(low_terms)):
        input_list = high_terms.copy()
        # Make one high
        input_list[i] = low_terms[i] 
        y_terms.append(estimateMinLength(*input_list))

    # High Case
    y_terms.append(estimateMinLength(*high_terms))

    #print(y_terms)

    term_significance_level = []
    term_is_significant = []
    # Perform Significance Analysis
    for i in range(len(low_terms)):
        term_significance_level.append(s(i+1, y_terms))
        term_is_significant.append(isSignificant(i+1, y_terms))

    #print(term_is_significant)

    print("Cotter's Method Results")
    print("{:16s} {:12s} {:11s}".format("Term", "Level", "Significant"))
    for i in range(len(low_terms)):
        print("{:<16s} {:<12.3f} {:<11s}".format(x_terms_names[i], term_significance_level[i]*100, str(term_is_significant[i])))

# Individual Scenario Simulation
else:
    print("Performing analysis on Specific Scenario")

    p_particle_s = 1000 #1000 # kg / m^3 
    d_particle_s =  2.5E-6 #[2.5E-6, 10.0E-6] # m
    q_particle_s = 2.0E-17 # C
    c_particles_s = 44.9E-9 #[29.7E-9 , 44.9E-9] # concentration kg/m^3
    h_particle_init_s = 0.35 # initial height above plate

    # fluid constants
    p_fluid_s = 1.293 # kg / m^3
    u_fluid_s = 1.81E-5 # Ns / m^2 p*v

    v_fluid_s = 2.5 # m/s

    # plate constants
    u_plate_s = 2000 #[200,2000] # Volts (potential applied to plate)

    scenario_terms = [p_particle_s, d_particle_s, q_particle_s, c_particles_s, \
            h_particle_init_s, p_fluid_s, u_fluid_s, u_plate_s, v_fluid_s]

    dist = estimateMinLength(*scenario_terms)

    print("\nInput Variables----------")
    for i in range(len(x_terms_names)):
        print("{:<16s}: {:<12.3e}".format(x_terms_names[i], scenario_terms[i]))
    print("\nOutput-------------------")
    print("minimum distance: {:.2f}".format(dist))
    
print("\n")

