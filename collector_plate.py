import math
import time
import matplotlib.pyplot as plt

# 10 and 2.5 micro meters

# constants
p_particle = 1000 # kg / m^3 
d_particle = 5.0E-6 # m
q_particle = 2.0E-17 # C
c_particles = 1.5E-11 # concentration
h_particle_init = 1.5E-3 # initial height above plate
v_particle_init = 0
m_particle = math.pi / 6 * p_particle * math.pow(d_particle, 3)

p_fluid = 1.293 # kg / m^3
u_fluid = 1.81E-5 #15.0E-6 * p_fluid # Ns / m^2 p*v

o_plate = 4.0E-6 # q / area c/m^2
H = 2 * h_particle_init # max height to look at surrounding charges

g = 9.81 # m/s^2
e0 = 8.854E-12 # F*m^-1

dt = 0.00002

# Gravitational Force (const) (downward)
f_gravity = - math.pi / 6.0 * p_particle * g * math.pow(d_particle, 3)

# Buoyant Force (cosnt) (upward)
f_buoyancy = math.pi / 6.0 * p_fluid * g * math.pow(d_particle, 3)

# Drag Force (depends on v apparent) (upward)
def Cd(v_apparent):
    print(p_fluid * d_particle * math.fabs(v_apparent) / u_fluid)
    return 24 / (p_fluid * d_particle * math.fabs(v_apparent) / u_fluid)

def fDrag(v_apparent):
    # return 0.5 * p_fluid * Cd(v_apparent) * math.pi * \
    #        math.pow(d_particle, 2) * math.pow(v_apparent, 2)
    return 6.0 * math.pi * u_fluid * d_particle / 2 * v_apparent

# Electric Field Force due to plate (const) (downward)
f_plate = - q_particle * o_plate / 4.0 / math.pi / e0

# Electric field force due to surrounding charges (depends on height)
def fSurrCharges(d):
    return math.pow(q_particle, 2) * c_particles / 2.0 / e0 * (2*d - H)

def totalForce(v, d):
    print("g: {:.3e} b: {:.3e} d: {:.3e} plate: {:.3e} surr: {:.3e}".format(\
        f_gravity, f_buoyancy, fDrag(v), f_plate, fSurrCharges(d)))
    return f_gravity + f_buoyancy + fDrag(v) + f_plate + fSurrCharges(d)

# model will start with initial v and d
# each time, calculate new v and d based on force calculation

current_h = h_particle_init
current_v = v_particle_init
current_a = totalForce(current_h, current_v) / m_particle
current_t = 0

# run model until the current h becomes <= 0 (hits the plate)
history_h = [current_h]
history_v = [current_v]
history_a = [current_a]
history_t = [current_t]

while(current_h > 0):
    current_a = totalForce(current_h, current_v) / m_particle
    # dv = a*dt
    current_v += current_a * dt
    # dh = v*dt
    current_h += current_v * dt
    current_t += dt

    print("T: {:5.3f} H: {:5.4e} V: {:5.4e} A: {:5.4e}: ".format(current_t, 
           current_h, current_v, current_a))

    history_h.append(current_h)
    history_v.append(current_v)
    history_a.append(current_a)
    history_t.append(current_t)
plt.subplot(1, 3, 1)
plt.plot(history_t, history_h)
plt.subplot(1, 3, 2)
plt.plot(history_t, history_v)
plt.subplot(1, 3, 3)
plt.plot(history_t, history_a)

plt.show()
    






