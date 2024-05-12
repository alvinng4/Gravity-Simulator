import argparse
import ctypes
import csv
import datetime
import math
from pathlib import Path
import platform
import re
import sys
import timeit

import matplotlib.pyplot as plt
import numpy as np

from simulator import Simulator
from progress_bar import Progress_bar


class Plotter:
    SIDEREAL_DAYS_PER_YEAR = 365.256363004

    def __init__(self):
        # --------------------Read command line arguments--------------------
        self._read_command_line_arg()

        # Use c library to perform simulation
        self.is_c_lib = self.args.numpy
        if self.is_c_lib:
            try:
                if platform.system() == "Windows":
                    self.c_lib = ctypes.cdll.LoadLibrary(
                        str(Path(__file__).parent / "c_lib.dll")
                    )
                elif platform.system() == "Darwin":
                    self.c_lib = ctypes.cdll.LoadLibrary(
                        str(Path(__file__).parent / "c_lib.dylib")
                    )
                elif platform.system() == "Linux":
                    self.c_lib = ctypes.cdll.LoadLibrary(
                        str(Path(__file__).parent / "c_lib.so")
                    )
            except:
                print("System message: Loading c_lib failed. Running with numpy.")
                self.is_c_lib = False

        self.store_every_n = self.args.store_every_n
        if self.store_every_n < 1:
            raise argparse.ArgumentTypeError(
                "Store every nth points should be larger than 1!"
            )

        # --------------------Initialize attributes--------------------
        self.simulator = Simulator(self)
        self.tolerance = None
        self.dt = None
        self.default_systems = [
            "circular_binary_orbit",
            "eccentric_binary_orbit",
            "3d_helix",
            "sun_earth_moon",
            "figure-8",
            "pyth-3-body",
            "solar_system",
            "solar_system_plus",
            "custom",
        ]
        self.available_integrators = [
            "euler",
            "euler_cromer",
            "rk4",
            "leapfrog",
            "rkf45",
            "dopri",
            "dverk",
            "rkf78",
            "ias15",
        ]
        self.available_integrators_to_printable_names = {
            "euler": "Euler",
            "euler_cromer": "Euler_Cromer",
            "rk4": "RK4",
            "leapfrog": "LeapFrog",
            "rkf45": "RKF45",
            "dopri": "DOPRI",
            "dverk": "DVERK",
            "rkf78": "RKF78",
            "ias15": "IAS15",
        }
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
            "Jupiter": "brown",
            "Saturn": "gold",
            "Uranus": "paleturquoise",
            "Neptune": "blue",
            "Moon": "grey",
            "Pluto": None,
            "Ceres": None,
            "Vesta": None,
        }
        self.recommended_settings = {
            # "template": ["tf", "tf unit", "tolerance"],
            "circular_binary_orbit": [50, "days", 1e-9],
            "eccentric_binary_orbit": [2.6, "years", 1e-9],
            "3d_helix": [20, "days", 1e-9],
            "sun_earth_moon": [1, "years", 1e-9],
            "figure-8": [20, "days", 1e-9],
            "pyth-3-body": [70, "days", 1e-9],
            "solar_system": [200, "years", 1e-9],
            "solar_system_plus": [250, "years", 1e-9],
        }

    def run_prog(self):
        # Catch KeyboardInterrupt
        try:
            # Restart program once the loop is finished.
            while True:
                print("\nGravity simulator")
                print("Exit the program anytime by hitting Ctrl + C\n")
                self._user_interface_before_simulation()
                if self.is_simulate == True:
                    self.computed_energy = False
                    self._launch_simulation()

                else:
                    self.computed_energy = True
                    (
                        self.simulator.objects_count,
                        self.simulator.sol_time,
                        self.simulator.sol_dt,
                        self.simulator.energy,
                        self.simulator.sol_state,
                    ) = self._read_simulation_data()

                self.data_size = len(self.simulator.sol_time)
                if self.data_size > 50000:
                    if self.ask_user_permission(
                        f"There are {self.data_size} lines of data. Do you want to trim the data?"
                    ):
                        self.trim_data()
                    else:
                        print()

                if self.tf_unit == "years":
                    self.sol_time_in_tf_unit = (
                        self.simulator.sol_time / self.SIDEREAL_DAYS_PER_YEAR
                    )
                else:
                    self.sol_time_in_tf_unit = self.simulator.sol_time

                self._user_interface_after_simulation()

        except KeyboardInterrupt:
            print("\nKeyboard Interrupt detected (Cltr + C). Exiting the program...")

    def _user_interface_before_simulation(self):
        while True:
            print("Select an action:")
            print("1. Launch simulation")
            print("2. Read simulation data")

            try:
                action = int(input("Enter action (Number): "))
                if action < 1 or action > 2:
                    raise ValueError
                else:
                    break
            except ValueError:
                print("Invalid input. Please try again.")
                print()
                continue

        print()
        match action:
            case 1:
                self.is_simulate = True
            case 2:
                self.is_simulate = False

    def _user_interface_after_simulation(self):
        while True:
            print("Select an action:")
            print("1. Plot 2D trajectory (xy plane)")
            print("2. Plot 3D trajectory")
            print("3. Plot relative energy error")
            print("4. Plot dt")
            print("5. Read data size")
            print("6. Trim data")
            print("7. Save simulation data")
            print("8. Compare relative energy error")
            print("9. Restart program")
            print("10. Exit")

            try:
                action = int(input("Enter action (Number): "))
                if action < 1 or action > 10:
                    raise ValueError
            except ValueError:
                print("Invalid input. Please try again.")
                print()
                continue

            print()
            match action:
                case 1:
                    self._plot_2d_trajectory()
                case 2:
                    self._plot_3d_trajectory()
                case 3:
                    if not self.computed_energy:
                        self.simulator.compute_energy()
                        self.computed_energy = True
                    self._plot_rel_energy()
                case 4:
                    self._plot_dt()
                case 5:
                    print(f"There are {self.data_size} lines of data.")
                    print()
                case 6:
                    print(f"There are {self.data_size} lines of data.")
                    self.trim_data()
                case 7:
                    if not self.computed_energy:
                        self.simulator.compute_energy()
                        self.computed_energy = True
                    self._save_result()
                case 8:
                    if not self.computed_energy:
                        self.simulator.compute_energy()
                        self.computed_energy = True
                    self._plot_compare_rel_energy()
                case 9:
                    break
                case 10:
                    print("Exiting the program...")
                    sys.exit(0)

    def _launch_simulation(self):
        while True:
            self._read_user_simulation_input()
            self._print_user_simulation_input()
            if self.ask_user_permission("Proceed?"):
                print("")
                break

        self.simulator.initialize_system(self)
        self.simulator.simulation()

    def _plot_2d_trajectory(self):
        print("Plotting 2D trajectory (xy plane)...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect="equal")

        ax.set_xlabel("$x$ (AU)")
        ax.set_ylabel("$y$ (AU)")
        
        # Get specific colors if the system is solar-like
        match self.system:
            case "sun_earth_moon":
                objs_name = ["Sun", "Earth", "Moon"]
            case "solar_system":
                objs_name = [
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
                objs_name = [
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

        if self.system in self.solar_like_systems:
            for i in range(self.simulator.objects_count):
                traj = ax.plot(
                    self.simulator.sol_state[:, i * 3],
                    self.simulator.sol_state[:, 1 + i * 3],
                    color=self.solar_like_systems_colors[objs_name[i]],
                )
                ax.plot(
                    self.simulator.sol_state[-1, i * 3],
                    self.simulator.sol_state[-1, 1 + i * 3],
                    "o",
                    color=traj[0].get_color(),
                    label=objs_name[i],
                )
        else:
            for i in range(self.simulator.objects_count):
                traj = ax.plot(
                    self.simulator.sol_state[:, i * 3],
                    self.simulator.sol_state[:, 1 + i * 3],
                )
                ax.plot(
                    self.simulator.sol_state[-1, i * 3],
                    self.simulator.sol_state[-1, 1 + i * 3],
                    "o",
                    color=traj[0].get_color(),
                )

        if self.system in self.solar_like_systems:
            fig.legend(loc="center right", borderaxespad=0.2)
            fig.tight_layout()

        plt.show()
        print()

    def _plot_3d_trajectory(self):
        print("Plotting 3D trajectory...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_box_aspect([1.0, 1.0, 1.0])  
        ax.set_xlabel("$x$ (AU)")
        ax.set_ylabel("$y$ (AU)")
        ax.set_zlabel("$z$ (AU)")

        # Get specific colors if the system is solar-like
        match self.system:
            case "sun_earth_moon":
                objs_name = ["Sun", "Earth", "Moon"]
            case "solar_system":
                objs_name = [
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
                objs_name = [
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

        if self.system in self.solar_like_systems:
            for i in range(self.simulator.objects_count):
                traj = ax.plot(
                    self.simulator.sol_state[:, i * 3],
                    self.simulator.sol_state[:, i * 3 + 1],
                    self.simulator.sol_state[:, i * 3 + 2],
                    color=self.solar_like_systems_colors[objs_name[i]],
                )
                ax.plot(
                    self.simulator.sol_state[-1, i * 3],
                    self.simulator.sol_state[-1, i * 3 + 1],
                    self.simulator.sol_state[-1, i * 3 + 2],
                    "o",
                    color=traj[0].get_color(),
                    label=objs_name[i],
                )
        else:
            for i in range(self.simulator.objects_count):
                traj = ax.plot(
                    self.simulator.sol_state[:, i * 3],
                    self.simulator.sol_state[:, i * 3 + 1],
                    self.simulator.sol_state[:, i * 3 + 2],
                )
                ax.plot(
                    self.simulator.sol_state[-1, i * 3],
                    self.simulator.sol_state[-1, i * 3 + 1],
                    self.simulator.sol_state[-1, i * 3 + 2],
                    "o",
                    color=traj[0].get_color(),
                )
        
        if self.system in self.solar_like_systems:
            fig.legend(loc="center right", borderaxespad=0.2)
            fig.tight_layout()

        self._set_3d_axes_equal(ax)
        plt.show()
        print()

    @staticmethod
    def _set_3d_axes_equal(ax):
        """
        Make axes of 3D plot have equal scale

        Input
        ax: a matplotlib axis, e.g., as output from plt.gca().

        Reference: karlo, https://stackoverflow.com/questions/13685386/how-to-set-the-equal-aspect-ratio-for-all-axes-x-y-z
        """

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5*max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
        
    def _plot_rel_energy(self):
        if not self.computed_energy:
            if self.ask_user_permission(
                "WARNING: Energy has not been computed. Compute energy?"
            ):
                self.simulator.compute_energy()
                self.computed_energy = True

        print("Plotting relative energy error...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.semilogy(
            self.sol_time_in_tf_unit,
            np.abs(
                (self.simulator.energy - self.simulator.energy[0])
                / self.simulator.energy[0]
            ),
        )
        ax.set_title("Relative energy error against time")
        ax.set_xlabel(f"Time ({self.tf_unit})")
        ax.set_ylabel("$|(E(t)-E_0)/E_0|$")

        plt.show()
        print()

    def _plot_tot_energy(self):
        # WARNING: The unit is in solar masses, AU and day
        print("Plotting total energy...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.semilogy(self.sol_time_in_tf_unit, np.abs(self.simulator.energy))
        ax.set_title("Total energy against time")
        ax.set_xlabel(f"Time ({self.tf_unit})")
        ax.set_ylabel("$E(t)$")

        plt.show()
        print()

    def _plot_dt(self):
        """
        Plot dt(days)
        """
        print("Plotting dt...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(self.sol_time_in_tf_unit, self.simulator.sol_dt, s=0.1)
        ax.set_yscale("log")
        ax.set_title("dt against time")
        ax.set_xlabel(f"Time ({self.tf_unit})")
        ax.set_ylabel("dt(days)")

        plt.show()
        print()

    def _plot_compare_rel_energy(self):
        if not self.computed_energy:
            if self.ask_user_permission(
                "WARNING: Energy has not been computed. Compute energy?"
            ):
                self.simulator.compute_energy()
                self.computed_energy = True

        filepaths_to_labels = {}
        filepaths_to_labels["current"] = input("Enter plot label for current set of data in memory: ")

        read_folder_path = Path(__file__).parent / "results"
        
        while True:
            num_of_data = len(filepaths_to_labels)
            if num_of_data >= 2:
                if not self.ask_user_permission(f"There are {num_of_data} sets of data. Continue adding new files?"):
                    break                    
            while True:
                print()
                read_file_path = input("Enter absolute path of new file to compare, or the complete file name if it is inside gravity_plot/results: ")
                read_file_path = read_folder_path / read_file_path

                if read_file_path.is_file():
                    label = input("Enter plot label for this file: ")
                    filepaths_to_labels[read_file_path] = label
                    break
                else:
                    print("File not found! Please try again.")

        print()
        print("Plotting relative energy error...(Please check the window)")

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(f"Time ({self.tf_unit})")
        ax.set_ylabel("$|(E(t)-E_0)/E_0|$")

        for filepath in filepaths_to_labels:
            if filepath == "current":
                ax.semilogy(
                    self.sol_time_in_tf_unit,
                    np.abs((self.simulator.energy - self.simulator.energy[0]) / self.simulator.energy[0]),
                    label=filepaths_to_labels[filepath],
                )
            else:
                try:
                    with open(filepath, "r") as file:
                        reader = csv.reader(file)

                        # Allocate memory
                        sol_time = np.zeros(50000)
                        energy = np.zeros(50000)

                        i = 0
                        for row in reader:
                            sol_time[i] = row[0]
                            energy[i] = row[2]
                            i += 1

                            # Extending memory buffer
                            if i % 50000 == 0:
                                sol_time = np.concatenate(
                                    (sol_time, np.zeros(50000))
                                )
                                energy = np.concatenate(
                                    (energy, np.zeros(50000))
                                )

                        sol_time = sol_time[:i]
                        energy = energy[:i]

                        if self.tf_unit == "years":
                            sol_time /= self.SIDEREAL_DAYS_PER_YEAR

                        ax.semilogy(
                            sol_time,
                            np.abs((energy - energy[0]) / energy[0]),
                            label=filepaths_to_labels[filepath],
                        )

                except FileNotFoundError:
                    sys.exit("Error: file is not found. Exiting the program")

        fig.legend(loc=7)
        fig.tight_layout()
        fig.subplots_adjust(right=0.8)
        plt.show()

        print()

    def _read_user_simulation_input(self):
        # --------------------Check and Input systems--------------------
        while True:
            self.available_systems = self.default_systems.copy()
            file_path = Path(__file__).parent / "customized_systems.csv"
            with open(file_path, "a"):  # Create file if not exist
                pass
            with open(file_path, "r+") as file:
                reader = csv.reader(file)
                for row in reader:
                    self.available_systems.append(row[0])

            print("Available systems:")
            for i, system in enumerate(self.available_systems):
                print(f"{i + 1}. {system}")
            self.system = input("Enter system (Number or name): ")
            # Temporary string
            temp_str = ""
            for i in range(len(self.available_systems)):
                temp_str += "|" + self.available_systems[i]
            if matches := re.search(
                rf"^\s*([1-9][0-9]*|{temp_str})\s*$",
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

            print("\nInvalid input. Please try again.\n")

        # --------------------Customize system--------------------
        if self.system == "custom":
            print("\nCustomizing system...")
            while True:
                system_name = input("Enter the name of the system: ")
                if matches := re.search(
                    r"^\s*(\S+)\s*$",
                    system_name,
                ):
                    if "," in matches.group(1):
                        print(
                            "Invalid input. Please do not use comma (,) inside the name."
                        )
                    elif matches.group(1) in self.available_systems:
                        print("System name already exist! Please try another one.")
                    else:
                        system_name = matches.group(1)
                        break

            self.system = system_name
            while True:
                try:
                    objects_count = int(input("Number of objects: ").strip())
                    if objects_count <= 0:
                        raise ValueError
                    else:
                        break
                except ValueError:
                    print("Invalid input. Please try again.")

            print("Note: The default unit is M_sun, AU and day, G=0.00029591220828411.")

            masses = []
            for i in range(objects_count):
                while True:
                    try:
                        masses.append(
                            float(
                                input(
                                    f"Please enter the mass for object {i + 1}: "
                                ).strip()
                            )
                        )
                        break
                    except ValueError:
                        print("Invalid input! Please try again.")
            state_vec = []
            for i in range(objects_count):
                for j in range(3):
                    match j:
                        case 0:
                            variable = "x"
                        case 1:
                            variable = "y"
                        case 2:
                            variable = "z"
                    while True:
                        try:
                            state_vec.append(
                                float(
                                    input(
                                        f"Please enter {variable} for object {i + 1}: "
                                    ).strip()
                                )
                            )
                            break
                        except ValueError:
                            print("Invalid input! Please try again.")
            for i in range(objects_count):
                for j in range(3):
                    match j:
                        case 0:
                            variable = "x"
                        case 1:
                            variable = "y"
                        case 2:
                            variable = "z"
                    while True:
                        try:
                            state_vec.append(
                                float(
                                    input(
                                        f"Please enter v{variable} for object {i + 1}: "
                                    ).strip()
                                )
                            )
                            break
                        except ValueError:
                            print("Invalid input! Please try again.")
            file_path = Path(__file__).parent / "customized_systems.csv"
            with open(file_path, "a", newline="") as file:
                writer = csv.DictWriter(
                    file,
                    fieldnames=["system_name", "objects_count", "masses", "state_vec"],
                )
                writer.writerow(
                    {
                        "system_name": system_name,
                        "objects_count": objects_count,
                        "masses": masses,
                        "state_vec": state_vec,
                    }
                )

        # --------------------Recommended settings for systems--------------------
        elif self.system in self.default_systems:
            print("")
            if self.ask_user_permission(
                "Do you want to use the recommended settings for this system?"
            ):
                print("")
                self.integrator = "ias15"
                self.tf, self.tf_unit, self.tolerance = self.recommended_settings[
                    self.system
                ]
                if self.tf_unit == "years":
                    self.tf *= self.SIDEREAL_DAYS_PER_YEAR

                return None

        print("")

        # --------------------Input integrators--------------------
        while True:
            print("Available integrators: ")
            for i, integrator in enumerate(self.available_integrators):
                print(
                    f"{i + 1}. {self.available_integrators_to_printable_names[integrator]}"
                )
            self.integrator = input("Enter integrator (Number or name): ")
            # Temporary string
            temp_str = ""
            for i in range(len(self.available_integrators)):
                temp_str += "|" + self.available_integrators[i]
            if matches := re.search(
                rf"^\s*([1-9]{temp_str})\s*$",
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

            print("\nInvalid input. Please try again.\n")

        # --------------------Input tf--------------------
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
                    self.tf_unit = "days"
                    break
                else:
                    self.tf_unit = "years"
                    self.tf *= self.SIDEREAL_DAYS_PER_YEAR
                    break

        # --------------------Input dt--------------------
        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            while True:
                print("")
                self.dt = input("Enter dt (d/yr): ")
                if matches := re.search(
                    r"([0-9]*\.?[0-9]*)(?:\.|\W*)*(day|year|d|y)?\s*",
                    self.dt,
                    re.IGNORECASE,
                ):
                    if matches.group(1):
                        self.dt = float(matches.group(1))
                    else:
                        print("\nInvalid input. Please try again.")
                        continue

                    if matches.group(2) not in ["year", "y"]:
                        self.dt_unit = "days"
                        break
                    else:
                        self.dt_unit = "years"
                        self.dt *= self.SIDEREAL_DAYS_PER_YEAR
                        break

        elif self.integrator in ["rkf45", "dopri", "dverk", "rkf78", "ias15"]:
            while True:
                try:
                    print("")
                    self.tolerance = float(input("Enter tolerance: "))
                    if self.tolerance <= 0:
                        raise ValueError
                    break
                except ValueError:
                    print("Invalid value. Please try again.")
        print("")

    def _print_user_simulation_input(self):
        print(f"System: {self.system}")
        print(
            f"Integrator: {self.available_integrators_to_printable_names[self.integrator]}"
        )
        if self.tf_unit == "years":
            print(f"tf: {self.tf / self.SIDEREAL_DAYS_PER_YEAR:g} years")
        else:
            print(f"tf: {self.tf} days")
        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            if self.dt_unit == "years":
                print(f"dt: {self.dt / self.SIDEREAL_DAYS_PER_YEAR} years")
            else:
                print(f"dt: {self.dt} days")
        elif self.integrator in ["rkf45", "dopri", "dverk", "rkf78", "ias15"]:
            print(f"Tolerance: {self.tolerance}")
        print(f"Use c_lib: {self.is_c_lib}")
        print(f"Store every nth point: {self.store_every_n}")

        if self.integrator in ["euler", "euler_cromer", "rk4", "leapfrog"]:
            npts = int(np.floor((self.tf / self.dt))) + 1  # + 1 for t0

            store_npts = npts
            if self.store_every_n != 1:
                store_npts = (
                    math.floor((npts - 1) / self.store_every_n) + 1
                )  # + 1 for t0

            print(f"Estimated number of points to be stored: {store_npts}")

        print("")

    def _read_command_line_arg(self):
        parser = argparse.ArgumentParser(description="N-body gravity simulator")
        parser.add_argument(
            "--numpy",
            "-n",
            action="store_false",
            help="disable c_lib and use numpy",
        )
        parser.add_argument(
            "--store_every_n",
            "-s",
            type=int,
            default=1,
            help="Store every nth points",
        )
        self.args = parser.parse_args()

    def trim_data(self):
        # --------------------Get user input--------------------
        while True:
            try:
                desired_trim_size = (
                    input('\nEnter desired data size (Enter "cancel" to cancel): ')
                    .strip()
                    .lower()
                )
                if desired_trim_size == "cancel":
                    print()
                    is_trim_data = False
                    break
                else:
                    desired_trim_size = int(desired_trim_size)
            except ValueError:
                print("Invalid input! Please try again.")
                continue

            if desired_trim_size > self.data_size:
                print("Value too big! Please try again.")
                continue
            if desired_trim_size == self.data_size:
                if self.ask_user_permission(
                    "The input value is equal to the original data size. Continue?"
                ):
                    print()
                    is_trim_data = False
                    break
                else:
                    continue
            elif desired_trim_size < 2:
                print("Value too small! Please try again.")
                continue
            else:
                divide_factor = math.ceil(self.data_size / desired_trim_size)
                trim_size = math.floor(self.data_size / divide_factor) + 1
                if self.ask_user_permission(
                    f"The trimmed data size would be {trim_size}. Continue?"
                ):
                    print()
                    is_trim_data = True
                    break

        # --------------------Trim data--------------------

        # Note: trim size is calculated in the the user input section:
        # divide_factor = math.ceil(self.data_size / desired_trim_size)
        # trim_size = math.floor(self.data_size / divide_factor) + 1

        if is_trim_data == True:
            trimmed_sol_time = np.zeros(trim_size)
            trimmed_sol_dt = np.zeros(trim_size)
            trimmed_sol_state = np.zeros(
                (trim_size, self.simulator.objects_count * 3 * 2)
            )
            if self.computed_energy == True:
                trimmed_energy = np.zeros(trim_size)

            j = 0
            for i in range(self.data_size):
                if i % divide_factor == 0:
                    trimmed_sol_time[j] = self.simulator.sol_time[i]
                    trimmed_sol_dt[j] = self.simulator.sol_dt[i]
                    trimmed_sol_state[j] = self.simulator.sol_state[i]
                    if self.computed_energy == True:
                        trimmed_energy[j] = self.simulator.energy[i]
                    j += 1

            if trimmed_sol_time[-1] != self.simulator.sol_time[-1]:
                trimmed_sol_time[-1] = self.simulator.sol_time[-1]
                trimmed_sol_dt[-1] = self.simulator.sol_dt[-1]
                trimmed_sol_state[-1] = self.simulator.sol_state[-1]
                if self.computed_energy == True:
                    trimmed_energy[-1] = self.simulator.energy[-1]

            self.data_size = len(trimmed_sol_time)
            print(f"Trimmed data size = {self.data_size}")
            print()

            self.simulator.sol_time = trimmed_sol_time
            self.simulator.sol_dt = trimmed_sol_dt
            self.simulator.sol_state = trimmed_sol_state
            if self.computed_energy == True:
                self.simulator.energy = trimmed_energy

            if self.tf_unit == "years":
                self.sol_time_in_tf_unit = (
                    self.simulator.sol_time / self.SIDEREAL_DAYS_PER_YEAR
                )
            else:
                self.sol_time_in_tf_unit = self.simulator.sol_time

    def _save_result(self):
        """
        Save the result in a csv file
        Unit: Solar masses, AU, day
        Format: time, dt, total energy, x1, y1, z1, x2, y2, z2, ... vx1, vy1, vz1, vx2, vy2, vz2, ...
        """
        if not self.computed_energy:
            if self.ask_user_permission(
                "WARNING: Energy has not been computed. Compute energy?"
            ):
                self.simulator.compute_energy()
                self.computed_energy = True

        print("Storing simulation results...")
        file_path = Path(__file__).parent / "results"
        file_path.mkdir(parents=True, exist_ok=True)
        if self.is_simulate:
            file_path = (
                Path(__file__).parent
                / "results"
                / (
                    str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
                    + f"_{self.system}_"
                    + f"{(self.tf / self.SIDEREAL_DAYS_PER_YEAR):g}{self.tf_unit[0]}_"
                    + f"{self.integrator}"
                    + ".csv"
                )
            )
        else:
            file_path = (
                Path(__file__).parent
                / "results"
                / (
                    str(datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
                    + f"_{self.system}"
                    + ".csv"
                )
            )

        if self.computed_energy:
            progress_bar = Progress_bar()
            with progress_bar:
                with open(file_path, "w", newline="") as file:
                    writer = csv.writer(file)
                    for count in progress_bar.track(range(self.data_size)):
                        row = np.insert(
                            self.simulator.sol_state[count],
                            0,
                            self.simulator.energy[count],
                        )
                        row = np.insert(row, 0, self.simulator.sol_dt[count])
                        row = np.insert(row, 0, self.simulator.sol_time[count])
                        writer.writerow(row.tolist())
            print(f"Storing completed. Please check {file_path}")
        else:
            if self.ask_user_permission(
                "WARNING: Energy has not been computed. The energy data will be stored as zeros. Proceed?"
            ):
                progress_bar = Progress_bar()
                with progress_bar:
                    with open(file_path, "w", newline="") as file:
                        writer = csv.writer(file)
                        for count in progress_bar.track(range(self.data_size)):
                            row = np.insert(self.simulator.sol_state[count], 0, 0)
                            row = np.insert(row, 0, self.simulator.sol_dt[count])
                            row = np.insert(row, 0, self.simulator.sol_time[count])
                            writer.writerow(row.tolist())

                print(f"Storing completed. Please check {file_path}")

        print("")

    @staticmethod
    def ask_user_permission(msg):
        while True:
            if matches := re.search(
                r"^\s*(yes|no|y|n)\s*$", input(f"{msg} (y/n): "), re.IGNORECASE
            ):
                if matches.group(1).lower() in ["y", "yes"]:
                    return True
                elif matches.group(1).lower() in ["n", "no"]:
                    return False

            print("Invalid input. Please try again.\n")

    def _read_simulation_data(self):
        self.system = "read_data"
        read_folder_path = Path(__file__).parent / "results"

        while True:
            read_file_path = input(
                "Enter absolute path to the file, or the complete file name if it is inside gravity_plot/results: "
            ).strip()
            read_file_path = read_folder_path / read_file_path

            if read_file_path.is_file():
                break
            else:
                print("File not found! Please try again.")

        # Get object_count
        try:
            with open(read_file_path, "r") as file:
                reader = csv.reader(file)
                objects_count = round((len(next(reader)) - 3) / (3 * 2))
        except FileNotFoundError:
            sys.exit("Error: file is not found. Exiting the program")

        # Read data
        start = timeit.default_timer()
        try:
            with open(read_file_path, "r") as file:
                reader = csv.reader(file)

                # Allocate memory
                sol_time = np.zeros(50000)
                sol_dt = np.zeros(50000)
                energy = np.zeros(50000)
                sol_state = np.zeros((50000, objects_count * 3 * 2))

                print()
                print("Reading data...")

                i = 0
                for row in reader:
                    sol_time[i] = row[0]
                    sol_dt[i] = row[1]
                    energy[i] = row[2]
                    for j in range(objects_count):
                        sol_state[i][j * 3 + 0] = row[3 + j * 3 + 0]
                        sol_state[i][j * 3 + 1] = row[3 + j * 3 + 1]
                        sol_state[i][j * 3 + 2] = row[3 + j * 3 + 2]
                        sol_state[i][objects_count * 3 + j * 3 + 0] = row[
                            3 + objects_count * 3 + j * 3 + 0
                        ]
                        sol_state[i][objects_count * 3 + j * 3 + 1] = row[
                            3 + objects_count * 3 + j * 3 + 1
                        ]
                        sol_state[i][objects_count * 3 + j * 3 + 2] = row[
                            3 + objects_count * 3 + j * 3 + 2
                        ]

                    i += 1

                    # Extending memory buffer
                    if i % 50000 == 0:
                        sol_time = np.concatenate((sol_time, np.zeros(50000)))
                        sol_dt = np.concatenate((sol_dt, np.zeros(50000)))
                        energy = np.concatenate((energy, np.zeros(50000)))
                        sol_state = np.concatenate(
                            (sol_state, np.zeros((50000, objects_count * 3 * 2)))
                        )

                stop = timeit.default_timer()
                print("Reading completed.\n")
                print(f"Run time: {(stop - start):.3f} s")
                print("")

                while True:
                    self.tf_unit = input("Enter tf unit for plotting (d/yr): ")
                    if matches := re.search(
                        r"(day|year|d|y)?", self.tf_unit, re.IGNORECASE
                    ):
                        if matches.group(1) not in ["year", "y"]:
                            self.tf_unit = "days"
                        else:
                            self.tf_unit = "years"

                        if self.ask_user_permission(
                            f"Unit for tf is {self.tf_unit}. Proceed?"
                        ):
                            break

                return (
                    objects_count,
                    sol_time[:i],
                    sol_dt[:i],
                    energy[:i],
                    sol_state[:i],
                )

        except FileNotFoundError:
            sys.exit("Error: file is not found. Exiting the program")

if __name__ == "__main__":
    plotter = Plotter()
    plotter.run_prog()
