# Euler Cromer method for 2 body simulation
import numpy as np
import matplotlib.pyplot as plt

# Gravitational constant (AU ^3/d^2/ M_sun):
G = 0.00029591220828559
# Simulation time (days)
# t0 = 0
# t1 = 365
dt = 0.5

# r1 - r3: Positions (AU)
# v1 - v3: Velocities (AU/d)
# m: Mass (Solar masses)
# a_i = - G M_j (ri - rj) / |r_ij|^3

class Simulator:

    def __init__(self, grav_sim):
        self.objects_count = grav_sim.stats.objects_count


    def initialize_problem(self, grav_sim):
        self.x = np.zeros((self.objects_count, 3))
        self.v = np.zeros((self.objects_count, 3))
        self.m = np.zeros(self.objects_count)
        for j in range(self.objects_count):
            self.x[j] = np.array([grav_sim.grav_objs.sprites()[j].params[f"r{i + 1}"] for i in range(3)])
            self.v[j] = np.array([grav_sim.grav_objs.sprites()[j].params[f"v{i + 1}"] for i in range(3)])
            self.m[j] = grav_sim.grav_objs.sprites()[j].params["m"]


    def ode_n_body_first_order(self):
        # Allocating memory
        self.a = self.x * 0

        # Differential equations:
        for j in range(0, self.objects_count):
            for k in range(0, self.objects_count):
                if j != k:
                    R = self.x[j] - self.x[k]
                    self.a[j] += - G * self.m[k] * R / np.linalg.norm(R) ** 3

    def Euler(self):
        for j in range(0, self.objects_count):
            self.x[j] = self.x[j] + self.v[j] * dt
            self.v[j] = self.v[j] + self.a[j] * dt
        

    def Euler_Cromer(self):
        for j in range(0, self.objects_count):
            self.v[j] = self.v[j] + self.a[j] * dt
            self.x[j] = self.x[j] + self.v[j] * dt










if __name__ == "__main__":
    pass
