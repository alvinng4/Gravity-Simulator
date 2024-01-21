import numpy as np
import numba as nb  # Note: nb.njit cannot works on functions inside a class

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
        self.stats.simulation_time += self.settings.dt
        if self.is_initialize == True:
            self.initialize_problem(grav_sim)

        match self.current_integrator:
            case "euler":
                self.is_initialize = False
                self.a = acceleration(self.stats.objects_count, self.x, self.m)
                self.x, self.v = euler(
                    self.x,
                    self.v,
                    self.a,
                    self.settings.dt,
                )
            case "euler_cromer":
                self.is_initialize = False
                self.a = acceleration(self.stats.objects_count, self.x, self.m)
                self.x, self.v = euler_cromer(
                    self.x,
                    self.v,
                    self.a,
                    self.settings.dt,
                )
            case "rk2":
                self.is_initialize = False
                self.a = acceleration(self.stats.objects_count, self.x, self.m)
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
                self.a = acceleration(self.stats.objects_count, self.x, self.m)
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
                    self.a = acceleration(self.stats.objects_count, self.x, self.m)
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
        self.is_leapfrog = False

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


# Note: jit cannot works on functions inside a class
@nb.njit
def acceleration(objects_count, x, m):
    """Calculate the acceleration"""
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

    k2_v = acceleration(objects_count, x_half, m)
    k2_x = v_half

    v = v + dt * k2_v
    x = x + dt * k2_x

    return x, v


@nb.njit
def rk4(objects_count, x, v, a, m, dt):
    k1_v = a
    k1_x = v

    k2_v = acceleration(objects_count, x + 0.5 * k1_x * dt, m)
    k2_x = v + 0.5 * k1_v * dt

    k3_v = acceleration(objects_count, x + 0.5 * k2_x * dt, m)
    k3_x = v + 0.5 * k2_v * dt

    k4_v = acceleration(objects_count, x + k3_x * dt, m)
    k4_x = v + k3_v * dt

    v = v + dt * (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6.0
    x = x + dt * (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6.0

    return x, v


@nb.njit
def leapfrog(objects_count, x, v, a, m, dt):
    a_0 = a
    x = x + v * dt + a_0 * 0.5 * dt * dt
    a_1 = acceleration(objects_count, x, m)
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


# Reference: see Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems
# Chapter 6, Page 100 - 101
def butcher_tableaus_rk(order):
    # Select integrator
    # 45) Runge-Kutta-Fehleberg 4(5)
    # 54) Dormand-Prince 5(4)
    # 78) Runge-Kutta-Fehlberg 7(8)
    # 65) Verner's method 6(5), DVERK

    # RUNGE-KUTTA-FEHLBERG 4(5)
    if order == 45:
        # Order
        p = 4
        pt = 5
        a = np.array(
            [
                [1.0 / 4.0, 0.0, 0.0, 0.0, 0.0],
                [3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0],
                [1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0],
                [439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0],
                [-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0],
            ]
        )
        c = np.array([1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5])
        b = np.array([25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -0.2, 0])
        bt = np.array(
            [
                16.0 / 135.0,
                0.0,
                6656.0 / 12825.0,
                28561.0 / 56430.0,
                -9.0 / 50.0,
                2.0 / 55.0,
            ]
        )

    # DORMAND-PRINCE 5(4)
    elif order == 54:
        # order
        p = 5
        pt = 4
        a = np.array(
            [
                [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0],
                [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0],
                [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0],
                [
                    19372.0 / 6561.0,
                    -25360.0 / 2187.0,
                    64448.0 / 6561.0,
                    -212.0 / 729.0,
                    0.0,
                    0,
                ],
                [
                    9017.0 / 3168.0,
                    -355.0 / 33.0,
                    46732.0 / 5247.0,
                    49.0 / 176.0,
                    -5103.0 / 18656.0,
                    0,
                ],
                [
                    35.0 / 384.0,
                    0.0,
                    500.0 / 1113.0,
                    125.0 / 192.0,
                    -2187.0 / 6784.0,
                    11.0 / 84.0,
                ],
            ]
        )
        c = np.array([1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0])
        b = np.array(
            [
                35.0 / 384.0,
                0.0,
                500.0 / 1113.0,
                125.0 / 192.0,
                -2187.0 / 6784.0,
                11.0 / 84.0,
                0,
            ]
        )
        bt = np.array(
            [
                5179.0 / 57600.0,
                0.0,
                7571.0 / 16695.0,
                393.0 / 640.0,
                -92097.0 / 339200.0,
                187.0 / 2100.0,
                1.0 / 40.0,
            ]
        )

    # RUNGE-KUTTA-FEHLBERG 7(8)
    elif order == 78:
        # Order
        p = 7
        pt = 8
        a = np.array(
            [
                [2.0 / 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [
                    1.0 / 36.0,
                    1.0 / 12.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    1.0 / 24.0,
                    0.0,
                    1.0 / 8.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    5.0 / 12.0,
                    0.0,
                    -25.0 / 16.0,
                    25.0 / 16.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    1.0 / 20.0,
                    0.0,
                    0.0,
                    1.0 / 4.0,
                    1.0 / 5.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    -25.0 / 108.0,
                    0.0,
                    0.0,
                    125.0 / 108.0,
                    -65.0 / 27.0,
                    125.0 / 54.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    31.0 / 300.0,
                    0.0,
                    0.0,
                    0.0,
                    61.0 / 225.0,
                    -2.0 / 9.0,
                    13.0 / 900.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    2.0,
                    0.0,
                    0.0,
                    -53.0 / 6.0,
                    704.0 / 45.0,
                    -107.0 / 9.0,
                    67.0 / 90.0,
                    3.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    -91.0 / 108.0,
                    0.0,
                    0.0,
                    23.0 / 108.0,
                    -976.0 / 135.0,
                    311.0 / 54.0,
                    -19.0 / 60.0,
                    17.0 / 6,
                    -1.0 / 12.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    2383.0 / 4100.0,
                    0.0,
                    0.0,
                    -341.0 / 164.0,
                    4496.0 / 1025.0,
                    -301.0 / 82.0,
                    2133.0 / 4100.0,
                    45.0 / 82.0,
                    45.0 / 164,
                    18.0 / 41.0,
                    0.0,
                    0.0,
                ],
                [
                    3.0 / 205.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    -6.0 / 41.0,
                    -3.0 / 205.0,
                    -3.0 / 41.0,
                    3.0 / 41.0,
                    6.0 / 41.0,
                    0.0,
                    0.0,
                ],
                [
                    -1777.0 / 4100.0,
                    0.0,
                    0.0,
                    -341.0 / 164.0,
                    4496.0 / 1025.0,
                    -289.0 / 82.0,
                    2193.0 / 4100.0,
                    51.0 / 82,
                    33.0 / 164,
                    19.0 / 41.0,
                    0.0,
                    1.0,
                ],
            ]
        )
        c = np.array(
            [
                2.0 / 27.0,
                1.0 / 9.0,
                1.0 / 6.0,
                5.0 / 12.0,
                1.0 / 2.0,
                5.0 / 6.0,
                1.0 / 6.0,
                2.0 / 3.0,
                1.0 / 3.0,
                1.0,
                0.0,
                1.0,
            ]
        )
        b = np.array(
            [
                41.0 / 840.0,
                0.0,
                0.0,
                0.0,
                0.0,
                34.0 / 105.0,
                9.0 / 35.0,
                9.0 / 35.0,
                9.0 / 280.0,
                9.0 / 280.0,
                41.0 / 840.0,
                0.0,
                0.0,
            ]
        )
        bt = np.array(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                34.0 / 105.0,
                9.0 / 35.0,
                9.0 / 35.0,
                9.0 / 280.0,
                9.0 / 280.0,
                0.0,
                41.0 / 840.0,
                41.0 / 840.0,
            ]
        )

    # VERNER 6(5) DVERK
    elif order == 65:
        # Order
        p = 6
        pt = 7
        a = np.array(
            [
                [1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [4.0 / 75.0, 16.0 / 75.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [5.0 / 6.0, -8.0 / 3.0, 5.0 / 2.0, 0.0, 0.0, 0.0, 0.0],
                [-165.0 / 64.0, 55.0 / 6.0, -425.0 / 64.0, 85.0 / 96.0, 0.0, 0.0, 0.0],
                [
                    12.0 / 5.0,
                    -8.0,
                    4015.0 / 612.0,
                    -11.0 / 36.0,
                    88.0 / 255.0,
                    0.0,
                    0.0,
                ],
                [
                    -8263.0 / 15000.0,
                    124.0 / 75.0,
                    -643.0 / 680.0,
                    -81.0 / 250.0,
                    2484.0 / 10625.0,
                    0.0,
                    0.0,
                ],
                [
                    3501.0 / 1720.0,
                    -300.0 / 43.0,
                    297275.0 / 52632.0,
                    -319.0 / 2322.0,
                    24068.0 / 84065.0,
                    0.0,
                    3850.0 / 26703.0,
                ],
            ]
        )
        c = np.array(
            [1.0 / 6.0, 4.0 / 15.0, 2.0 / 3.0, 5.0 / 6.0, 1.0, 1.0 / 15.0, 1.0]
        )
        b = np.array(
            [
                3.0 / 40.0,
                0.0,
                875.0 / 2244.0,
                23.0 / 72.0,
                264.0 / 1955.0,
                0.0,
                125.0 / 11592.0,
                43.0 / 616.0,
            ]
        )
        bt = np.array(
            [
                13.0 / 160.0,
                0.0,
                2375.0 / 5984.0,
                5.0 / 16.0,
                12.0 / 85.0,
                3.0 / 44.0,
                0.0,
                0.0,
            ]
        )

    return a, b, bt, c, p, pt
