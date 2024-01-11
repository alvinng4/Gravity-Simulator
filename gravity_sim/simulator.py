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
def ode_n_body_first_order(objects_count, x, v, m):
    # Allocating memory
    a = np.zeros((objects_count, 3))

    # Differential equations:
    for j in range(0, objects_count):
        for k in range(0, objects_count):
            if j != k:
                R = x[j] - x[k]
                a[j] += -G * m[k] * R / np.linalg.norm(R) ** 3

    return x, v, a, m


@nb.njit
def euler(objects_count, x, v, a, dt=0.001):
    for j in range(0, objects_count):
        x[j] = x[j] + v[j] * dt
        v[j] = v[j] + a[j] * dt
    return x, v


@nb.njit
def euler_cromer(objects_count, x, v, a, dt=0.001):
    for j in range(0, objects_count):
        v[j] = v[j] + a[j] * dt
        x[j] = x[j] + v[j] * dt
    return x, v

    # def RK4(self):
    # for j in range(0, self.objects_count):

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


if __name__ == "__main__":
    pass
