import csv
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np

from simulator import Simulator


class Plotter:
    def __init__(self):
        # Use c library to perform simulation
        self.is_ctypes = True

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
            "template": ["tf", "tf unit", "tolerance"],
            "circular_binary_orbit": [50, "days", 1e-6],
            "eccentric_binary_orbit": [2.6, "years", 1e-6],
            "3d_helix": [20, "days", 1e-6],
            "sun_earth_moon": [1, "years", 1e-9],
            "figure-8": [20, "days", 1e-6],
            "pyth-3-body": [70, "days", 1e-6],
            "solar_system": [200, "years", 1e-6],
            "solar_system_plus": [250, "years", 1e-6],
        }

    def run_prog(self):
        # Catch KeyboardInterrupt
        try:
            # Restart once all the progress is finished
            while True:
                # Main program

                # Read user input
                while True:
                    print("\nGravity simulator")
                    print("Exit the program anytime by hitting Ctrl + C\n")
                    self._read_user_input()
                    self._print_user_input()
                    if self.ask_user_permission("Proceed?"):
                        print("")
                        break

                # Launch simulation
                self.simulator = Simulator(self)
                self.simulator.initialize_system(self)
                if self.is_ctypes == True:
                    self.simulator.simulation_ctypes()
                else:
                    self.simulator.simulation_numpy()
                if self.unit == "years":
                    self.simulator.sol_time /= 365.24

                if len(self.simulator.sol_time) > 20000:
                    if self.ask_user_permission(
                        f"There are {len(self.simulator.sol_time)} lines of data. Do you want to trim the data?"
                    ):
                        self.simulator.trim_data()

                # Plot the result
                if self.ask_user_permission("Plot trajectory?"):
                    print("")
                    self._plot_trajectory()

                self.is_compute_energy = False
                if self.ask_user_permission("Compute energy?"):
                    print("")
                    self.is_compute_energy = True
                    self.simulator.compute_energy()

                    if self.ask_user_permission("Plot relative energy error?"):
                        print("")
                        self._plot_rel_energy()
                # self._plot_tot_energy() # Warnings: The unit is in solar masses, AU and day.

                # Store data
                print(f"Lines of data = {len(self.simulator.sol_time)}")
                if self.ask_user_permission("Save simulation data?"):
                    print("")
                    self.simulator.save_result(self.is_compute_energy)

                # Ask permission to restart the whole program
                if not self.ask_user_permission(
                    "All plotting is done. Restart simulation?"
                ):
                    print("Exiting the program...")
                    break

        except KeyboardInterrupt:
            print("\nKeyboard Interrupt detected (Cltr + C). Exiting the program...")

    def _plot_trajectory(self):
        """
        Plot the trajectory
        """
        print("Plotting trajectory...(Please check the window)\n")
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect="equal")
        if self.system in self.solar_like_systems:
            # Plot trajectory and initial positions with the same color:
            for i in range(self.simulator.objects_count):
                traj = ax1.plot(
                    self.simulator.sol_state[:, i * 3],
                    self.simulator.sol_state[:, 1 + i * 3],
                    color=self.solar_like_systems_colors[self.simulator.objs_name[i]],
                )
                # Check if the system have names for individual objects
                ax1.plot(
                    self.simulator.sol_state[-1, i * 3],
                    self.simulator.sol_state[-1, 1 + i * 3],
                    "o",
                    color=traj[0].get_color(),
                    label=self.simulator.objs_name[i],
                )
        else:
            # Plot trajectory and initial positions with the same color:
            for i in range(self.simulator.objects_count):
                traj = ax1.plot(
                    self.simulator.sol_state[:, i * 3],
                    self.simulator.sol_state[:, 1 + i * 3],
                )
                # Check if the system have names for individual objects
                ax1.plot(
                    self.simulator.sol_state[-1, i * 3],
                    self.simulator.sol_state[-1, 1 + i * 3],
                    "o",
                    color=traj[0].get_color(),
                )

        ax1.set_title("Trajectory")
        ax1.set_xlabel("X (AU)")
        ax1.set_ylabel("Y (AU)")

        if self.system in self.solar_like_systems:
            fig1.legend(loc=7)
            fig1.tight_layout()
            fig1.subplots_adjust(right=0.8)

        plt.show()

    def _plot_rel_energy(self):
        """
        Plot the relative energy error
        """
        print("Plotting relative energy error...(Please check the window)")
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.semilogy(
            self.simulator.sol_time,
            np.abs(
                (self.simulator.energy - self.simulator.energy[0])
                / self.simulator.energy[0]
            ),
        )
        ax2.set_title("Relative energy error against time")
        ax2.set_xlabel(f"Time ({self.unit})")
        ax2.set_ylabel("|(E(t)-E0)/E0|")

        plt.show()
        print("")

    def _plot_tot_energy(self):
        """
        Plot the total energy
        Warning: The unit is in solar masses, AU and day
        """
        print("Plotting total energy...(Please check the window)")
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.semilogy(self.simulator.sol_time, np.abs(self.simulator.energy))
        ax3.set_title("Total energy against time")
        ax3.set_xlabel(f"Time ({self.unit})")
        ax3.set_ylabel("E(t)")

        plt.show()

    def _read_user_input(self):
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

        elif self.system in self.default_systems:
            print("")
            if self.ask_user_permission(
                "Do you want to use the recommended settings for this system?"
            ):
                print("")
                self.integrator = "rkf78"
                self.tf, self.unit, self.tolerance = self.recommended_settings[
                    self.system
                ]
                if self.unit == "years":
                    self.tf *= 365.24

                return None

        print("")
        while True:
            print("Available integrators: ")
            for i, integrator in enumerate(self.available_integrators):
                print(f"{i + 1}. {integrator}")
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
                        self.dt *= 365.24
                        break

        elif self.integrator in ["rkf45", "dopri", "dverk", "rkf78"]:
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

    def _print_user_input(self):
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
        elif self.integrator in ["rkf45", "dopri", "dverk", "rkf78"]:
            print(f"tolerance: {self.tolerance}")
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


if __name__ == "__main__":
    plotter = Plotter()
    plotter.run_prog()
