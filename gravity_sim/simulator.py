# Euler Cromer method for 2 body simulation
import numpy as np
import matplotlib.pyplot as plt

# Gravitational constant (AU ^3/d^2/ M_sun):
G = 0.00029591220828559
# Simulation time (days)
# t0 = 0
# t1 = 365
dt = 0.5

# r1 - r3: Positions (AU), v1 - v3: Velocities (AU/d), m: Mass (Solar masses)


def simulator(grav_objs):
    sun = grav_objs.sprites()[0]
    earth = grav_objs.sprites()[3]
    
    x_s = np.array([sun.params["r1"], sun.params["r2"], sun.params["r3"]])
    x_e = np.array([earth.params["r1"], earth.params["r2"], earth.params["r3"]])
    v_s = np.array([sun.params["v1"], sun.params["v2"], sun.params["v3"]])
    v_e = np.array([earth.params["v1"], earth.params["v2"], earth.params["v3"]])
    m_s = sun.params["m"]
    m_e = earth.params["m"]

    a_e = -G * m_s * (x_e - x_s) / np.linalg.norm(x_e - x_s) ** 3
    a_s = -G * m_e * (x_s - x_e) / np.linalg.norm(x_e - x_s) ** 3

    v_e = v_e + a_e * dt
    x_e = x_e + v_e * dt
    
    v_s = v_s + a_s * dt
    x_s = x_s + v_s * dt
    
    return x_s, v_s, x_e, v_e


if __name__ == "__main__":
    pass
