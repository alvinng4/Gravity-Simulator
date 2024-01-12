import numpy as np
import matplotlib.pyplot as plt
import numba as nb

# Gravitational constant (AU ^3/d^2/ M_sun):
G = 0.00029591220828559
# G = 1.0  # For Testing


# dt: Simulation time (days)
# r1 - r3: Positions (AU)
# v1 - v3: Velocities (AU/d)
# m: Mass (Solar masses)
# a_i = - G M_j (ri - rj) / |r_ij|^3


class Simulator:
    def __init__(self, grav_sim):
        self.stats = grav_sim.stats
        self.settings = grav_sim.settings

        self.m = []
        self.x = []
        self.v = []
        self.a = []

        self.is_initialize = True
        self.set_all_integrators_false()
        self.is_leapfrog = True  # Default integrator
        self.current_integrator = "leapfrog"

    def run_simulation(self, grav_sim):
        self.stats.simulation_time += self.stats.dt
        if self.is_initialize == True:
            self.initialize_problem(grav_sim)

        match self.current_integrator:
            case "euler":
                self.is_initialize = False
                self.a = ode_n_body_first_order(
                    self.stats.objects_count, self.x, self.m
                )
                self.x, self.v = euler(
                    self.x,
                    self.v,
                    self.a,
                    self.settings.dt,
                )
            case "euler_cromer":
                self.is_initialize = False
                self.a = ode_n_body_first_order(
                    self.stats.objects_count, self.x, self.m
                )
                self.x, self.v = euler_cromer(
                    self.x,
                    self.v,
                    self.a,
                    self.settings.dt,
                )
            case "rk2":
                self.is_initialize = False
                self.a = ode_n_body_first_order(
                    self.stats.objects_count, self.x, self.m
                )
                self.x, self.v = rk2(
                    self.stats.objects_count,
                    self.x,
                    self.v,
                    self.a,
                    self.m,
                    self.settings.dt,
                )
            case "rk4":
                self.is_initialize = False
                self.a = ode_n_body_first_order(
                    self.stats.objects_count, self.x, self.m
                )
                self.x, self.v = rk4(
                    self.stats.objects_count,
                    self.x,
                    self.v,
                    self.a,
                    self.m,
                    self.settings.dt,
                )
            case "leapfrog":
                if self.is_initialize == True:
                    self.a = ode_n_body_first_order(
                        self.stats.objects_count, self.x, self.m
                    )
                    self.is_initialize = False

                self.x, self.v, self.a = leapfrog(
                    self.stats.objects_count,
                    self.x,
                    self.v,
                    self.a,
                    self.m,
                    self.settings.dt,
                )
        self.stats.total_energy = total_energy(
            self.stats.objects_count, self.x, self.v, self.m
        )

    def unload_value(self, grav_sim):
        for j in range(self.stats.objects_count):
            grav_sim.grav_objs.sprites()[j].params["r1"] = self.x[j][0]
            grav_sim.grav_objs.sprites()[j].params["r2"] = self.x[j][1]
            grav_sim.grav_objs.sprites()[j].params["r3"] = self.x[j][2]
            grav_sim.grav_objs.sprites()[j].params["v1"] = self.v[j][0]
            grav_sim.grav_objs.sprites()[j].params["v2"] = self.v[j][1]
            grav_sim.grav_objs.sprites()[j].params["v3"] = self.v[j][2]

    def set_all_integrators_false(self):
        self.is_euler = False
        self.is_euler_cromer = False
        self.is_rk2 = False
        self.is_rk4 = False
        self.is_leapfrog = True

    def check_current_integrator(self):
        if self.is_euler == True:
            self.current_integrator = "euler"
        elif self.is_euler_cromer == True:
            self.current_integrator = "euler_cromer"
        elif self.is_rk2 == True:
            self.current_integrator = "rk2"
        elif self.is_rk4 == True:
            self.current_integrator = "rk4"
        elif self.is_leapfrog == True:
            self.current_integrator = "leapfrog"

    def initialize_problem(self, grav_sim):
        objects_count = grav_sim.stats.objects_count
        self.x = np.zeros((objects_count, 3))
        self.v = np.zeros((objects_count, 3))
        self.m = np.zeros(objects_count)
        for j in range(objects_count):
            self.x[j] = np.array(
                [grav_sim.grav_objs.sprites()[j].params[f"r{i + 1}"] for i in range(3)]
            )
            self.v[j] = np.array(
                [grav_sim.grav_objs.sprites()[j].params[f"v{i + 1}"] for i in range(3)]
            )
            self.m[j] = grav_sim.grav_objs.sprites()[j].params["m"]


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
    x_half, v_half = euler(x, v, a, 0.5 * dt)

    k2_v = ode_n_body_first_order(objects_count, x_half, m)
    k2_x = v_half

    v = v + dt * k2_v
    x = x + dt * k2_x

    return x, v


@nb.njit
def rk4(objects_count, x, v, a, m, dt):
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
def leapfrog(objects_count, x, v, a, m, dt):
    a_0 = a
    x = x + v * dt + a_0 * 0.5 * dt * dt
    a_1 = ode_n_body_first_order(objects_count, x, m)
    v = v + (a_0 + a_1) * 0.5 * dt

    return x, v, a_1


@nb.njit
def total_energy(objects_count, x, v, m):
    E = 0
    for j in range(0, objects_count):
        E += 0.5 * m[j] * np.linalg.norm(v[j]) ** 2
        for k in range(0, objects_count):
            if j != k:
                R = x[j] - x[k]
                E += -G * m[j] * m[k] / np.linalg.norm(R)
    return E


def test():
    # Initialize
    R1 = np.array([1.0, 0.0, 0.0])
    R2 = np.array([-1.0, 0.0, 0.0])
    V1 = np.array([0.0, 0.5, 0.0])
    V2 = np.array([0.0, -0.5, 0.0])
    x = np.zeros((2, 3))
    v = np.zeros((2, 3))
    x[0] = R1
    x[1] = R2
    v[0] = V1
    v[1] = V2
    m = [1.0, 1.0]
    t0 = 0.0
    tf = 10000.0
    dt = 0.001

    # Simulation
    npts = int(np.floor((tf - t0) / dt)) + 1
    sol_time = np.linspace(t0, t0 + dt * (npts - 1), npts)
    energy = np.zeros(npts)
    for count, t in enumerate(sol_time):
        a = Simulator.ode_n_body_first_order(2, x, m)
        x, v = Simulator.leapfrog(2, x, v, a, m, dt)
        energy[count] = Simulator.total_energy(2, x, v, m)

    # Plotting
    plt.figure()
    plt.semilogy(sol_time, np.abs((energy - energy[0]) / energy[0]))
    plt.xlabel("Time")
    plt.ylabel("|(E(t)-E0)/E0|")
    plt.show()

    plt.figure()
    plt.semilogy(sol_time, np.abs(energy))
    plt.xlabel("Time")
    plt.ylabel("E(t)")
    plt.show()


if __name__ == "__main__":
    test()
