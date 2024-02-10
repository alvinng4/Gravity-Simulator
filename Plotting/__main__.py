from os import environ

# Remove the "Hello from the pygame community." message when starting the program.
environ["PYGAME_HIDE_SUPPORT_PROMPT"] = "1"
import math
from pathlib import Path
import sys

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
sys.path.insert(0, path + "/gravity_sim")
import timeit
import re

import matplotlib.pyplot as plt
import numba as nb
import numpy as np

from gravity_sim import simulator
from gravity_sim.grav_obj import Grav_obj


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
            "3d_helix",
            "sun_earth_moon",
            "figure-8",
            "pyth-3-body",
            "solar_system",
            "solar_system_plus",
        ]
        self.solar_like_systems = [
            "sun_earth_moon",
            "solar_system",
            "solar_system_plus",
        ]
        self.solar_like_systems_colors = {
            "Sun": "orange",
            "Mercury": "slategrey",
            "Venus": "wheat",
            "Earth": "skyblue",
            "Mars": "red",
            "Jupiter": None,
            "Saturn": None,
            "Uranus": "paleturquoise",
            "Neptune": "blue",
            "Moon": "grey",
            "Pluto": None,
            "Ceres": None,
            "Vesta": None,
        }

    def run_prog(self):
        try:
            while True:
                while True:
                    print("\nN-Body gravity simulator")
                    print("Exit the program anytime by hitting Ctrl + C")
                    self._read_user_input()
                    self._print_user_input()
                    if self._ask_user_permission("Proceed?"):
                        break

                self._initialize_system()
                print("")
                print("Simulating the system...")
                self._simulation()
                print("")
                print("Computing energy...")
                self._compute_energy()
                print("")
                print("Plotting...")
                if self.unit == "years":
                    self.sol_time /= 365.24
                self._plot_trajectory()
                self._plot_rel_energy()
                self._plot_tot_energy()

                if not self._ask_user_permission(
                    "All plotting is done. Restart simulation?"
                ):
                    print("Exiting the program...")
                    break

        except KeyboardInterrupt:
            print("")
            print("Keyboard Interrupt detected (Cltr + C). Exiting the program...")

    def _read_user_input(self):
        while True:
            print("")
            print("Available systems:")
            for i, system in enumerate(self.available_systems):
                print(f"{i + 1}. {system}")
            self.system = input("Enter system (Number or name): ")
            if matches := re.search(
                r"^\s*([1-9]|circular_binary_orbit|3d_helix|sun_earth_moon|figure-8|pyth-3-body|solar_system_plus|solar_system)\s*$",
                self.system,
                re.IGNORECASE,
            ):
                if matches.group(1):
                    try:
                        int(matches.group(1))
                    except ValueError:
                        if matches.group(1).lower() in self.available_systems:
                            self.system = matches.group(1).lower()
                            break
                    else:
                        if (int(matches.group(1)) - 1) in range(
                            len(self.available_systems)
                        ):
                            self.system = self.available_systems[
                                int(matches.group(1)) - 1
                            ]
                            break

            print("Invalid input. Please try again.")

        while True:
            print("")
            print("Available integrators: ")
            for i, integrator in enumerate(self.available_integrators):
                print(f"{i + 1}. {integrator}")
            self.integrator = input("Enter integrator (Number or name): ")
            if matches := re.search(
                r"^\s*([1-9]|euler_cromer|euler|rk4|leapfrog|rkf45|dopri|rkf78|dverk)\s*$",
                self.integrator,
                re.IGNORECASE,
            ):
                if matches.group(1):
                    try:
                        int(matches.group(1))
                    except ValueError:
                        if matches.group(1).lower() in self.available_integrators:
                            self.integrator = matches.group(1).lower()
                            break
                    else:
                        if (int(matches.group(1)) - 1) in range(
                            len(self.available_integrators)
                        ):
                            self.integrator = self.available_integrators[
                                int(matches.group(1)) - 1
                            ]
                            break

            print("Invalid input. Please try again.")

        while True:
            print("")
            self.tf = input("Enter tf (d/yr): ")
            if matches := re.search(
                r"([0-9]*\.?[0-9]*)(?:\.|\W*)*(day|year|d|y)?", self.tf, re.IGNORECASE
            ):
                if matches.group(1):
                    self.tf = float(matches.group(1))
                else:
                    print("Invalid input. Please try again.")
                    continue

                if matches.group(2) not in ["year", "y"]:
                    self.unit = "days"
                    break
                else:
                    self.unit = "years"
                    self.tf *= 365.24
                    break

        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            while True:
                print("")
                self.dt = input("Enter dt (d/yr): ")
                if matches := re.search(
                    r"([0-9]*\.?[0-9]*)(?:\.|\W*)*(day|year|d|y)?",
                    self.dt,
                    re.IGNORECASE,
                ):
                    if matches.group(1):
                        self.dt = float(matches.group(1))
                    else:
                        print("Invalid input. Please try again.")
                        continue

                    if matches.group(2) not in ["year", "y"]:
                        self.dt_unit = "days"
                        break
                    else:
                        self.dt_unit = "years"
                        self.dt *= 365.24
                        break

        elif self.integrator in ["rkf45", "dopri", "rkf78", "dverk"]:
            while True:
                try:
                    print("")
                    self.tolerance = float(input("Enter tolerance: "))
                    if self.tolerance <= 0:
                        raise ValueError
                    break
                except ValueError:
                    print("Invalid value. Please try again.")

    def _print_user_input(self):
        print("")
        print(f"System: {self.system}")
        print(f"Integrator: {self.integrator}")
        if self.unit == "years":
            print(f"tf: {self.tf / 365.24} years")
        else:
            print(f"tf: {self.tf} days")
        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            if self.dt_unit == "years":
                print(f"dt: {self.dt / 365.24} years")
            else:
                print(f"dt: {self.dt} days")
        elif self.integrator in ["rkf45", "dopri", "rkf78", "dverk"]:
            print(f"tolerance: {self.tolerance}")

    @staticmethod
    def _ask_user_permission(msg):
        while True:
            if matches := re.search(
                r"^\s*(yes|no|y|n)$", input(f"{msg} (Y/N): "), re.IGNORECASE
            ):
                if matches.group(1).lower() in ["y", "yes"]:
                    return True
                elif matches.group(1).lower() in ["n", "no"]:
                    return False

            print("Invalid input. Please try again.")

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
                self.m = [1.0 / Grav_obj.G, 1.0 / Grav_obj.G]
                self.objects_count = 2

            case "3d_helix":
                R1 = np.array([0.0, 0.0, -1.0])
                R2 = np.array([-math.sqrt(3.0) / 2.0, 0.0, 0.5])
                R3 = np.array([math.sqrt(3.0) / 2.0, 0.0, 0.5])
                v0 = math.sqrt(1.0 / math.sqrt(3))
                V1 = np.array([-v0, 0.5, 0.0])
                V2 = np.array([0.5 * v0, 0.5, (math.sqrt(3.0) / 2.0) * v0])
                V3 = np.array([0.5 * v0, 0.5, -(math.sqrt(3.0) / 2.0) * v0])
                self.x = np.array([R1, R2, R3])
                self.v = np.array([V1, V2, V3])
                self.m = [1.0 / Grav_obj.G, 1.0 / Grav_obj.G, 1.0 / Grav_obj.G]
                self.objects_count = 3

            case "sun_earth_moon":
                self.m = [
                    Grav_obj.SOLAR_SYSTEM_MASSES["Sun"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Earth"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Moon"],
                ]
                R_CM = (
                    1
                    / np.sum(self.m)
                    * (
                        self.m[0] * np.array(Grav_obj.SOLAR_SYSTEM_POS["Sun"])
                        + self.m[1] * np.array(Grav_obj.SOLAR_SYSTEM_POS["Earth"])
                        + self.m[2] * np.array(Grav_obj.SOLAR_SYSTEM_POS["Moon"])
                    )
                )
                V_CM = (
                    1
                    / np.sum(self.m)
                    * (
                        self.m[0] * np.array(Grav_obj.SOLAR_SYSTEM_VEC["Sun"])
                        + self.m[1] * np.array(Grav_obj.SOLAR_SYSTEM_VEC["Earth"])
                        + self.m[2] * np.array(Grav_obj.SOLAR_SYSTEM_VEC["Moon"])
                    )
                )
                R1 = np.array(Grav_obj.SOLAR_SYSTEM_POS["Sun"] - R_CM)
                R2 = np.array(Grav_obj.SOLAR_SYSTEM_POS["Earth"] - R_CM)
                R3 = np.array(Grav_obj.SOLAR_SYSTEM_POS["Moon"] - R_CM)
                V1 = np.array(Grav_obj.SOLAR_SYSTEM_VEC["Sun"] - V_CM)
                V2 = np.array(Grav_obj.SOLAR_SYSTEM_VEC["Earth"] - V_CM)
                V3 = np.array(Grav_obj.SOLAR_SYSTEM_VEC["Moon"] - V_CM)
                self.x = np.array([R1, R2, R3])
                self.v = np.array([V1, V2, V3])

                self.objects_count = 3
                self.objs_name = ["Sun", "Earth", "Moon"]

            case "figure-8":
                R1 = np.array([0.970043, -0.24308753, 0.0])
                R2 = np.array([-0.970043, 0.24308753, 0.0])
                R3 = np.array([0.0, 0.0, 0.0])
                V1 = np.array([0.466203685, 0.43236573, 0.0])
                V2 = np.array([0.466203685, 0.43236573, 0.0])
                V3 = np.array([-0.93240737, -0.86473146, 0.0])
                self.x = np.array([R1, R2, R3])
                self.v = np.array([V1, V2, V3])
                self.m = [1.0 / Grav_obj.G, 1.0 / Grav_obj.G, 1.0 / Grav_obj.G]
                self.objects_count = 3

            case "pyth-3-body":
                R1 = np.array([1.0, 3.0, 0.0])
                R2 = np.array([-2.0, -1.0, 0.0])
                R3 = np.array([1.0, -1.0, 0.0])
                V1 = np.array([0.0, 0.0, 0.0])
                V2 = np.array([0.0, 0.0, 0.0])
                V3 = np.array([0.0, 0.0, 0.0])
                self.x = np.array([R1, R2, R3])
                self.v = np.array([V1, V2, V3])
                self.m = [3.0 / Grav_obj.G, 4.0 / Grav_obj.G, 5.0 / Grav_obj.G]
                self.objects_count = 3

            case "solar_system":
                R1 = Grav_obj.SOLAR_SYSTEM_POS["Sun"]
                R2 = Grav_obj.SOLAR_SYSTEM_POS["Mercury"]
                R3 = Grav_obj.SOLAR_SYSTEM_POS["Venus"]
                R4 = Grav_obj.SOLAR_SYSTEM_POS["Earth"]
                R5 = Grav_obj.SOLAR_SYSTEM_POS["Mars"]
                R6 = Grav_obj.SOLAR_SYSTEM_POS["Jupiter"]
                R7 = Grav_obj.SOLAR_SYSTEM_POS["Saturn"]
                R8 = Grav_obj.SOLAR_SYSTEM_POS["Uranus"]
                R9 = Grav_obj.SOLAR_SYSTEM_POS["Neptune"]

                V1 = Grav_obj.SOLAR_SYSTEM_VEC["Sun"]
                V2 = Grav_obj.SOLAR_SYSTEM_VEC["Mercury"]
                V3 = Grav_obj.SOLAR_SYSTEM_VEC["Venus"]
                V4 = Grav_obj.SOLAR_SYSTEM_VEC["Earth"]
                V5 = Grav_obj.SOLAR_SYSTEM_VEC["Mars"]
                V6 = Grav_obj.SOLAR_SYSTEM_VEC["Jupiter"]
                V7 = Grav_obj.SOLAR_SYSTEM_VEC["Saturn"]
                V8 = Grav_obj.SOLAR_SYSTEM_VEC["Uranus"]
                V9 = Grav_obj.SOLAR_SYSTEM_VEC["Neptune"]

                self.x = np.array([R1, R2, R3, R4, R5, R6, R7, R8, R9])
                self.v = np.array([V1, V2, V3, V4, V5, V6, V7, V8, V9])
                self.m = [
                    Grav_obj.SOLAR_SYSTEM_MASSES["Sun"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Mercury"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Venus"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Earth"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Mars"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Jupiter"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Saturn"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Uranus"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Neptune"],
                ]
                self.objects_count = 9
                self.objs_name = [
                    "Sun",
                    "Mercury",
                    "Venus",
                    "Earth",
                    "Mars",
                    "Jupiter",
                    "Saturn",
                    "Uranus",
                    "Neptune",
                ]

            case "solar_system_plus":
                R1 = Grav_obj.SOLAR_SYSTEM_POS["Sun"]
                R2 = Grav_obj.SOLAR_SYSTEM_POS["Mercury"]
                R3 = Grav_obj.SOLAR_SYSTEM_POS["Venus"]
                R4 = Grav_obj.SOLAR_SYSTEM_POS["Earth"]
                R5 = Grav_obj.SOLAR_SYSTEM_POS["Mars"]
                R6 = Grav_obj.SOLAR_SYSTEM_POS["Jupiter"]
                R7 = Grav_obj.SOLAR_SYSTEM_POS["Saturn"]
                R8 = Grav_obj.SOLAR_SYSTEM_POS["Uranus"]
                R9 = Grav_obj.SOLAR_SYSTEM_POS["Neptune"]
                R10 = Grav_obj.SOLAR_SYSTEM_POS["Pluto"]
                R11 = Grav_obj.SOLAR_SYSTEM_POS["Ceres"]
                R12 = Grav_obj.SOLAR_SYSTEM_POS["Vesta"]

                V1 = Grav_obj.SOLAR_SYSTEM_VEC["Sun"]
                V2 = Grav_obj.SOLAR_SYSTEM_VEC["Mercury"]
                V3 = Grav_obj.SOLAR_SYSTEM_VEC["Venus"]
                V4 = Grav_obj.SOLAR_SYSTEM_VEC["Earth"]
                V5 = Grav_obj.SOLAR_SYSTEM_VEC["Mars"]
                V6 = Grav_obj.SOLAR_SYSTEM_VEC["Jupiter"]
                V7 = Grav_obj.SOLAR_SYSTEM_VEC["Saturn"]
                V8 = Grav_obj.SOLAR_SYSTEM_VEC["Uranus"]
                V9 = Grav_obj.SOLAR_SYSTEM_VEC["Neptune"]
                V10 = Grav_obj.SOLAR_SYSTEM_VEC["Pluto"]
                V11 = Grav_obj.SOLAR_SYSTEM_VEC["Ceres"]
                V12 = Grav_obj.SOLAR_SYSTEM_VEC["Vesta"]

                self.x = np.array([R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12])
                self.v = np.array([V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12])
                self.m = [
                    Grav_obj.SOLAR_SYSTEM_MASSES["Sun"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Mercury"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Venus"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Earth"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Mars"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Jupiter"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Saturn"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Uranus"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Neptune"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Pluto"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Ceres"],
                    Grav_obj.SOLAR_SYSTEM_MASSES["Vesta"],
                ]
                self.objects_count = 12
                self.objs_name = [
                    "Sun",
                    "Mercury",
                    "Venus",
                    "Earth",
                    "Mars",
                    "Jupiter",
                    "Saturn",
                    "Uranus",
                    "Neptune",
                    "Pluto",
                    "Ceres",
                    "Vesta",
                ]

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
            self.progress_percentage = 0

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
                    if ((count + 1) / self.npts) * 100 > self.progress_percentage:
                        self.progress_percentage = int(((count + 1) / self.npts) * 100)
                        self._progress_bar(self.progress_percentage)
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
                    if ((count + 1) / self.npts) * 100 > self.progress_percentage:
                        self.progress_percentage = int(((count + 1) / self.npts) * 100)
                        self._progress_bar(self.progress_percentage)
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
                    if ((count + 1) / self.npts) * 100 > self.progress_percentage:
                        self.progress_percentage = int(((count + 1) / self.npts) * 100)
                        self._progress_bar(self.progress_percentage)
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
                    if ((count + 1) / self.npts) * 100 > self.progress_percentage:
                        self.progress_percentage = int(((count + 1) / self.npts) * 100)
                        self._progress_bar(self.progress_percentage)

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
                self.sol_state, self.sol_time = self._rk_embedded(
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

    def _plot_trajectory(self):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect="equal")
        if self.system in self.solar_like_systems:
            # Plot trajectory and initial positions with the same color:
            for i in range(self.objects_count):
                traj = ax1.plot(
                    self.sol_state[:, i * 3],
                    self.sol_state[:, 1 + i * 3],
                    color=self.solar_like_systems_colors[self.objs_name[i]],
                )
                # Check if the system have names for individual objects
                ax1.plot(
                    self.sol_state[-1, i * 3],
                    self.sol_state[-1, 1 + i * 3],
                    "o",
                    color=traj[0].get_color(),
                    label=self.objs_name[i],
                )
        else:
            # Plot trajectory and initial positions with the same color:
            for i in range(self.objects_count):
                traj = ax1.plot(self.sol_state[:, i * 3], self.sol_state[:, 1 + i * 3])
                # Check if the system have names for individual objects
                ax1.plot(
                    self.sol_state[-1, i * 3],
                    self.sol_state[-1, 1 + i * 3],
                    "o",
                    color=traj[0].get_color(),
                )

        ax1.set_xlabel("X (AU)")
        ax1.set_ylabel("Y (AU)")

        if self.system in self.solar_like_systems:
            fig1.legend(loc=7)
            fig1.tight_layout()
            fig1.subplots_adjust(right=0.8)

        plt.show()

    def _plot_rel_energy(self):
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.semilogy(
            self.sol_time, np.abs((self.energy - self.energy[0]) / self.energy[0])
        )
        ax2.set_xlabel(f"Time ({self.unit})")
        ax2.set_ylabel("|(E(t)-E0)/E0|")
        plt.show()

    def _plot_tot_energy(self):
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.semilogy(self.sol_time, np.abs(self.energy))
        ax3.set_xlabel(f"Time ({self.unit})")
        ax3.set_ylabel("E(t)")
        plt.show()

    def _compute_energy(self):
        self.energy = np.zeros(len(self.sol_state))
        self.progress_percentage = 0
        self.npts = len(self.sol_state)

        start = timeit.default_timer()
        for count in range(self.npts):
            x = self.sol_state[count]
            for i in range(self.objects_count):
                # KE
                self.energy[count] += (
                    0.5
                    * self.m[i]
                    * np.linalg.norm(
                        x[
                            (self.objects_count + i)
                            * 3 : (self.objects_count + 1 + i)
                            * 3
                        ]
                    )
                    ** 2
                )
                # PE
                for j in range(self.objects_count):
                    if i < j:
                        self.energy[count] -= (
                            Grav_obj.G
                            * self.m[i]
                            * self.m[j]
                            / np.linalg.norm(
                                x[i * 3 : (i + 1) * 3] - x[j * 3 : (j + 1) * 3]
                            )
                        )

            if ((count + 1) / self.npts) * 100 > self.progress_percentage:
                self.progress_percentage = int(((count + 1) / self.npts) * 100)
                self._progress_bar(self.progress_percentage)
        stop = timeit.default_timer()
        print(f"Run time:{(stop - start):.3f} s")

    @staticmethod
    def _progress_bar(percentage):
        if percentage != 100:
            fill = "█" * int(percentage / 2)
            bar = fill + "-" * (50 - int(percentage / 2))
            print(f"\r|{bar}| {percentage:3}% Completed ", end="")
        else:
            bar = "█" * 50
            print(f"\r|{bar}| 100% Completed ")

    @staticmethod
    def _rk_embedded(
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

        # Progress percentage for progress bar
        # Type: int
        progress_percentage = 0

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
            sum = np.sum(
                np.square(error_estimation_delta_x / tolerance_scale_x)
            ) + np.sum(np.square(error_estimation_delta_v / tolerance_scale_v))
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
                    sol_state = np.concatenate(
                        (sol_state, np.zeros((npts, objects_count * 2 * 3)))
                    )
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

            if (t / tf * 100) > (progress_percentage + 1):
                progress_percentage = math.floor(t / tf * 100)
                Plotter._progress_bar(progress_percentage)

            if t >= tf:
                Plotter._progress_bar(100)
                return sol_state[0:count], sol_time[0:count]


if __name__ == "__main__":
    plotter = Plotter()
    plotter.run_prog()
