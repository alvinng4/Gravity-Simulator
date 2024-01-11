import numpy as np
import matplotlib.pyplot as plt
import numba as nb

# Gravitational constant (AU ^3/d^2/ M_sun):
G = 0.00029591220828559
# Simulation time (days)
# t0 = 0
# t1 = 365

# r1 - r3: Positions (AU)
# v1 - v3: Velocities (AU/d)
# m: Mass (Solar masses)
# a_i = - G M_j (ri - rj) / |r_ij|^3


def initialize_problem(grav_sim, x, v, m):
    objects_count = grav_sim.stats.objects_count
    if len(m) == objects_count:
        pass
    else:
        x = np.zeros((objects_count, 3))
        v = np.zeros((objects_count, 3))
        m = np.zeros(objects_count)
    for j in range(objects_count):
        x[j] = np.array(
            [grav_sim.grav_objs.sprites()[j].params[f"r{i + 1}"] for i in range(3)]
        )
        v[j] = np.array(
            [grav_sim.grav_objs.sprites()[j].params[f"v{i + 1}"] for i in range(3)]
        )
        m[j] = grav_sim.grav_objs.sprites()[j].params["m"]

    return x, v, m


@nb.njit
def ode_n_body_first_order(objects_count, x, m):
    # Allocating memory
    a = np.zeros((objects_count, 3))

    # Differential equations:
    for j in range(0, objects_count):
        for k in range(0, objects_count):
            if j != k:
                R = x[j] - x[k]
                a[j] += -G * m[k] * R / np.linalg.norm(R) ** 3

    return a


@nb.njit
def euler(x, v, a, dt=0.001):
    return x + v * dt, v + a * dt


@nb.njit
def euler_cromer(x, v, a, dt=0.001):
    v = v + a * dt
    x = x + v * dt
    return x, v

@nb.njit
def rk2(objects_count, x, v, a, m, dt):
    # POOR PERFORMANCE NEED FIX
    x_half, v_half = euler(x, v, a, 0.5 * dt)

    k2_v = ode_n_body_first_order(objects_count, x_half, m)
    k2_x = v_half

    v = v + dt * k2_v 
    x = x + dt * k2_x

    return x, v

@nb.njit
def rk4(objects_count, x, v, a, m, dt):
    # POOR PERFORMANCE NEED FIX
    k1_v = a
    k1_x = v

    v1 = v + 0.5 * k1_v * dt
    x1 = x + 0.5 * k1_x * dt
    k2_v = ode_n_body_first_order(objects_count, x1, m)
    k2_x = v1

    v2 = v + 0.5 * k2_v * dt
    x2 = x + 0.5 * k2_x * dt
    k3_v = ode_n_body_first_order(objects_count, x2, m)
    k3_x = v2

    v3 = v + k3_v * dt
    x3 = x + k3_x * dt
    k4_v = ode_n_body_first_order(objects_count, x3, m)
    k4_x = v3


    v = v + dt * (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6.0
    x = x + dt * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6.0

    return x, v

@nb.njit 
def total_energy(objects_count, x, v, m):
    E = 0
    for j in range(0, objects_count):
        E += 0.5 * m[j] * np.linalg.norm(v[j])**2
        for k in range(0, objects_count):
            if j != k:  
                R = x[j] - x[k]
                E += - G * m[j] * m[k] / np.linalg.norm(R)
    return E




