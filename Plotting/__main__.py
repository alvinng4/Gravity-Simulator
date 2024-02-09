from pathlib import Path
import sys

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
import timeit

import numba as nb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from gravity_sim import simulator

# Gravitational constant (AU ^3/d^2/ M_sun):
G = 0.00029591220828559


class Plotter:
    def __init__(self):
        self.available_integrators = [
            "euler",
            "euler_cromer",
            "rk4",
            "leapfrog",
            "rkf45",
            "dopri",
            "rkf78",
            "dverk",
        ]
        self.available_systems = [
            "circular_binary_orbit",
            "sun_earth",
            "figure-8",
            "solar_system",
        ]
        while True:
            print("Available integrators: ")
            print(self.available_integrators)
            self.integrator = input("Enter integrator: ").strip().lower()
            if self.integrator not in self.available_integrators:
                print("Invalid integrator. Please try again.")
                continue
            break

        while True:
            try:
                self.tf = float(input("Enter tf (days): "))
                if self.tf <= 0:
                    raise ValueError
                break
            except ValueError:
                print("Invalid value. Please try again.")

        while True:
            print("Available systems:")
            print(self.available_systems)
            self.system = input("Enter a system: ").strip().lower()
            if self.system not in self.available_systems:
                print("Invalid system. Please try again.")
                continue
            break

        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            while True:
                try:
                    self.dt = float(input("Enter dt (days): "))
                    if self.dt <= 0:
                        raise ValueError
                    break
                except ValueError:
                    print("Invalid value. Please try again.")

        elif self.integrator in ["rkf45", "dopri", "rkf78", "dverk"]:
            while True:
                try:
                    self.tolerance = float(input("Enter tolerance: "))
                    if self.tolerance <= 0:
                        raise ValueError
                    break
                except ValueError:
                    print("Invalid value. Please try again.")

    def run_prog(self):
        self._initialize_system()
        self._simulation()
        self._plot_trajectory()

    def _initialize_system(self):
        self.t0 = 0.0
        if self.integrator in ["rkf45", "dopri", "rkf78", "dverk"]:
            match self.integrator:
                case "rkf45":
                    self.order = 45
                case "dopri":
                    self.order = 54
                case "rkf78":
                    self.order = 78
                case "dverk":
                    self.order = 65

        match self.system:
            case "circular_binary_orbit":
                R1 = np.array([1.0, 0.0, 0.0])
                R2 = np.array([-1.0, 0.0, 0.0])
                V1 = np.array([0.0, 0.5, 0.0])
                V2 = np.array([0.0, -0.5, 0.0])
                self.x = np.array([R1, R2])
                self.v = np.array([V1, V2])
                self.m = [1.0 / G, 1.0 / G]
                self.objects_count = 2

            case "sun_earth":
                R1 = np.array(
                    [4.980992013803802e-07, -2.910880021966631e-06, 1.649263431256283e-10]
                )
                R2 = np.array(
                    [-0.165850747934970, 0.969227871585065, -0.000054915079739]
                )
                V1 = np.array(
                    [5.176134557401363e-08, 8.891729869525767e-09, -2.054090638389034e-12]
                )
                V2 = np.array(
                    [-0.017234835658800, -0.002960655317675, 0.000000683945021]
                )
                self.x = np.array([R1, R2])
                self.v = np.array([V1, V2])
                self.m = [
                    1.0,
                    3.00329789031573e-06,
                ]
                self.objects_count = 2

            case "figure-8":
                R1 = np.array([0.970043, -0.24308753, 0.0])
                R2 = np.array([-0.970043, 0.24308753, 0.0])
                R3 = np.array([0.0, 0.0, 0.0])
                V1 = np.array([0.466203685, 0.43236573, 0.0])
                V2 = np.array([0.466203685, 0.43236573, 0.0])
                V3 = np.array([-0.93240737, -0.86473146, 0.0])
                self.x = np.array([R1, R2, R3])
                self.v = np.array([V1, V2, V3])
                self.m = [1.0 / G, 1.0 / G, 1.0 / G]
                self.objects_count = 3

            case "solar_system":
                R1 = np.array(
                    [-0.007967955691534, -0.002906227441573, 0.000210305430155]
                )
                R2 = np.array(
                    [-0.282598326953863, 0.197455979595808, 0.041774335580637]
                )
                R3 = np.array(
                    [-0.723210370166638, -0.079483020263124, 0.040428714281743]
                )
                R4 = np.array(
                    [-0.173819201725705, 0.966324555023514, 0.000155390185490]
                )
                R5 = np.array(
                    [-0.301326239258265, -1.454029331393295, -0.023005314339914]
                )
                R6 = np.array(
                    [3.485202469657675, 3.552136904413157, -0.092710354427984]
                )
                R7 = np.array(
                    [8.988104223143450, -3.719064854634689, -0.293193777732359]
                )
                R8 = np.array(
                    [12.263024178975050, 15.297387924805450, -0.102054902688356]
                )
                R9 = np.array(
                    [29.835014609847410, -1.793812957956852, -0.650640113225459]
                )

                V1 = np.array(
                    [0.000004875094764, -0.000007057133214, -0.000000045734537]
                )
                V2 = np.array(
                    [-0.022321659001897, -0.021572071031763, 0.000285519341050]
                )
                V3 = np.array(
                    [0.002034068201002, -0.020208286265930, -0.000394563984386]
                )
                V4 = np.array(
                    [-0.017230012325382, -0.002967721342619, 0.000000638212538]
                )
                V5 = np.array(
                    [0.014248322593453, -0.001579236181581, -0.000382372279616]
                )
                V6 = np.array(
                    [-0.005470970658852, 0.005642487338479, 0.000098961906021]
                )
                V7 = np.array(
                    [0.001822013845554, 0.005143470425888, -0.000161723590489]
                )
                V8 = np.array(
                    [-0.003097615358317, 0.002276781932346, 0.000048604332222]
                )
                V9 = np.array(
                    [0.000167653661182, 0.003152098732862, -0.000068775010957]
                )

                self.x = np.array([R1, R2, R3, R4, R5, R6, R7, R8, R9])
                self.v = np.array([V1, V2, V3, V4, V5, V6, V7, V8, V9])
                self.m = [
                    1.0,
                    1.66051140935277e-07,
                    2.44827371182131e-06,
                    3.00329789031573e-06,
                    3.22773848604808e-07,
                    0.000954532562518104,
                    0.00028579654259599,
                    4.3655207025844e-05,
                    5.1499991953912e-05,
                ]
                self.objects_count = 9

        # Prevent the error message from numba package:
        # "Encountered the use of a type that is scheduled for deprecation: type 'reflected list' found for argument 'm' of function '...'."
        self.m = nb.typed.List(self.m)

    def _simulation(self):
        start = timeit.default_timer()

        # Simulation
        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            self.npts = int(np.floor((self.tf - self.t0) / self.dt)) + 1
            self.sol_state = np.zeros((self.npts, self.objects_count * 3 * 2))
            self.sol_time = np.linspace(
                self.t0, self.t0 + self.dt * (self.npts - 1), self.npts
            )
            self.energy = np.zeros(self.npts)

        match self.integrator:
            case "euler":
                for count in range(self.npts):
                    self.a = simulator.acceleration(self.objects_count, self.x, self.m)
                    self.x, self.v = simulator.euler(self.x, self.v, self.a, self.dt)
                    self.sol_state[count] = np.concatenate(
                        (
                            np.reshape(self.x, self.objects_count * 3),
                            np.reshape(self.v, self.objects_count * 3),
                        )
                    )
                    self.energy[count] = simulator.total_energy(
                        self.objects_count, self.x, self.v, self.m
                    )
            case "euler_cromer":
                for count in range(self.npts):
                    self.a = simulator.acceleration(self.objects_count, self.x, self.m)
                    self.x, self.v = simulator.euler_cromer(
                        self.x, self.v, self.a, self.dt
                    )
                    self.sol_state[count] = np.concatenate(
                        (
                            np.reshape(self.x, self.objects_count * 3),
                            np.reshape(self.v, self.objects_count * 3),
                        )
                    )
                    self.energy[count] = simulator.total_energy(
                        self.objects_count, self.x, self.v, self.m
                    )
            case "rk4":
                for count in range(self.npts):
                    self.x, self.v = simulator.rk4(
                        self.objects_count, self.x, self.v, self.m, self.dt
                    )
                    self.sol_state[count] = np.concatenate(
                        (
                            np.reshape(self.x, self.objects_count * 3),
                            np.reshape(self.v, self.objects_count * 3),
                        )
                    )
                    self.energy[count] = simulator.total_energy(
                        self.objects_count, self.x, self.v, self.m
                    )
            case "leapfrog":
                self.a = simulator.acceleration(self.objects_count, self.x, self.m)
                for count in range(self.npts):
                    self.x, self.v, self.a = simulator.leapfrog(
                        self.objects_count, self.x, self.v, self.a, self.m, self.dt
                    )
                    self.sol_state[count] = np.concatenate(
                        (
                            np.reshape(self.x, self.objects_count * 3),
                            np.reshape(self.v, self.objects_count * 3),
                        )
                    )
                    self.energy[count] = simulator.total_energy(
                        self.objects_count, self.x, self.v, self.m
                    )

            case "rkf45" | "dopri" | "rkf78" | "dverk":
                (
                    self.power,
                    self.power_test,
                    self.coeff,
                    self.weights,
                    self.weights_test,
                ) = simulator.butcher_tableaus_rk(self.order)
                self.a = simulator.acceleration(self.objects_count, self.x, self.m)
                self.initial_dt = simulator._initial_time_step_rk_embedded(
                    self.objects_count,
                    self.power,
                    self.x,
                    self.v,
                    self.a,
                    self.m,
                    self.tolerance,
                    self.tolerance,
                )
                self.sol_state, self.sol_time = plotting_rk_embedded(
                    self.objects_count,
                    self.x,
                    self.v,
                    self.m,
                    self.tf,
                    self.initial_dt,
                    self.power,
                    self.power_test,
                    self.coeff,
                    self.weights,
                    self.weights_test,
                    self.tolerance,
                    self.tolerance,
                )

        stop = timeit.default_timer()
        print(f"Run time: {stop - start:.3f} s")
        print("Plotting...")

    def _plot_trajectory(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect="equal")
        # Plot trajectory and initial positions with the same color:
        for ibody in range(self.objects_count):
            traj = ax.plot(
                self.sol_state[:, ibody * 3], self.sol_state[:, 1 + ibody * 3]
            )
            ax.plot(
                self.sol_state[-1, ibody * 3],
                self.sol_state[-1, 1 + ibody * 3],
                "o",
                color=traj[0].get_color(),
            )
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        plt.show()


# def test_two_vectors(integrator: str, tf: float, dt: float = None, tolerance:float = None):
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(1, 1, 1)
# ax1.semilogy(sol_time, np.abs((energy - energy[0]) / energy[0]))
# ax1.set_xlabel("Time")
# ax1.set_ylabel("|(E(t)-E0)/E0|")
# plt.show()

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(1, 1, 1)
# ax2.semilogy(sol_time, np.abs(energy))
# ax2.set_xlabel("Time")
# ax2.set_ylabel("E(t)")
# ax2.yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
# ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot (111, aspect='equal ')
# ax.plot(sol_state [:,0], sol_state [:,1],"b-")
# ax.plot(sol_state [:,0 + 3], sol_state [:,1 + 3],"g-")
# plt.show()


@nb.njit
def plotting_rk_embedded(
    objects_count: int,
    x,
    v,
    m,
    tf: float,
    initial_dt,
    power,
    power_test,
    coeff,
    weights,
    weights_test,
    abs_tolerance: float,
    rel_tolerance: float,
):
    """
    Perform simulation using rk_embedded methods
    Modified for plotting

    :return: sol_state, sol_time
    :rtype: numpy.array
    """
    # Initializing
    t = 0.0
    dt = initial_dt
    stages = len(weights)
    min_power = min([power, power_test])
    error_estimation_delta_weights = weights - weights_test

    # Safety factors for step-size control:
    safety_fac_max = 6.0
    safety_fac_min = 0.33
    safety_fac = 0.38 ** (1.0 / (1.0 + min_power))

    # Initialize vk and xk
    vk = np.zeros((stages, objects_count, 3))
    xk = np.zeros((stages, objects_count, 3))

    # Allocate for dense output:
    npts = 100000
    sol_state = np.zeros((npts, objects_count * 2 * 3))
    sol_time = np.zeros(npts)

    # Initial values
    sol_state[0] = np.concatenate(
        (np.reshape(x, objects_count * 3), np.reshape(v, objects_count * 3))
    )
    sol_time[0] = t

    # Launch integration:
    count = 0
    while True:
        # Calculate xk and vk
        vk[0] = simulator.acceleration(objects_count, x, m)
        xk[0] = np.copy(v)
        for stage in range(1, stages):
            temp_v = np.zeros((objects_count, 3))
            temp_x = np.zeros((objects_count, 3))
            for j in range(stage):
                temp_v += coeff[stage - 1][j] * vk[j]
                temp_x += coeff[stage - 1][j] * xk[j]
            vk[stage] = simulator.acceleration(objects_count, x + dt * temp_x, m)
            xk[stage] = v + dt * temp_v

        # Calculate x_1, v_1 and also delta x, delta v for error estimation
        temp_v = np.zeros((objects_count, 3))
        temp_x = np.zeros((objects_count, 3))
        error_estimation_delta_x = np.zeros((objects_count, 3))
        error_estimation_delta_v = np.zeros((objects_count, 3))
        for stage in range(stages):
            temp_v += weights[stage] * vk[stage]
            temp_x += weights[stage] * xk[stage]
            error_estimation_delta_v += (
                error_estimation_delta_weights[stage] * vk[stage]
            )
            error_estimation_delta_x += (
                error_estimation_delta_weights[stage] * xk[stage]
            )
        v_1 = v + dt * temp_v
        x_1 = x + dt * temp_x
        error_estimation_delta_v *= dt
        error_estimation_delta_x *= dt

        # Error calculation
        tolerance_scale_v = (
            abs_tolerance + np.maximum(np.abs(v), np.abs(v_1)) * rel_tolerance
        )
        tolerance_scale_x = (
            abs_tolerance + np.maximum(np.abs(x), np.abs(x_1)) * rel_tolerance
        )

        # Sum up all the elements of x/tol and v/tol, square and divide by the total number of elements
        sum = np.sum(np.square(error_estimation_delta_x / tolerance_scale_x)) + np.sum(
            np.square(error_estimation_delta_v / tolerance_scale_v)
        )
        error = np.sqrt(sum / (objects_count * 3 * 2))

        if error <= 1 or dt == tf * 1e-12:
            t += dt
            x = x_1
            v = v_1
            count += 1

            # Store step:
            sol_state[count] = np.concatenate(
                (np.reshape(x, objects_count * 3), np.reshape(v, objects_count * 3))
            )
            sol_time[count] = t

            # Check buffer size and extend if needed :
            if (count + 1) == len(sol_state):
                sol_state = np.concatenate((sol_state, np.zeros((npts, objects_count * 2 * 3))))
                sol_time = np.concatenate((sol_time, np.zeros(npts)))

        dt_new = dt * safety_fac / error ** (1.0 / (1.0 + min_power))
        # Prevent dt to be too small or too large relative to the last time step
        if dt_new > safety_fac_max * dt:
            dt = safety_fac_max * dt
        elif dt_new < safety_fac_min * dt:
            dt = safety_fac_min * dt
        elif dt_new / tf < 1e-12:
            dt = tf * 1e-12
        else:
            dt = dt_new

        # Correct overshooting:
        if t + dt > tf:
            dt = tf - t

        if t >= tf:
            return sol_state[0:count], sol_time[0:count]


if __name__ == "__main__":
    plotter = Plotter()
    plotter.run_prog()
