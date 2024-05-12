import csv
import math
from pathlib import Path
import sys 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter

class Plotter:

    @staticmethod
    def plot_2d_trajectory(grav_plot):
        print("Plotting 2D trajectory (xy plane)...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect="equal")

        ax.set_xlabel("$x$ (AU)")
        ax.set_ylabel("$y$ (AU)")
        
        # Get specific colors if the system is solar-like
        match grav_plot.system:
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

        if grav_plot.system in grav_plot.solar_like_systems:
            # Plot with solar system colors
            for i in range(grav_plot.simulator.objects_count):
                traj = ax.plot(
                    grav_plot.simulator.sol_state[:, i * 3],
                    grav_plot.simulator.sol_state[:, 1 + i * 3],
                    color=grav_plot.solar_like_systems_colors[objs_name[i]],
                )
                # Plot the last position as a filled circle
                ax.plot(
                    grav_plot.simulator.sol_state[-1, i * 3],
                    grav_plot.simulator.sol_state[-1, 1 + i * 3],
                    "o",
                    color=traj[0].get_color(),
                    label=objs_name[i],
                )
        else:
            # Plot with random colors
            for i in range(grav_plot.simulator.objects_count):
                traj = ax.plot(
                    grav_plot.simulator.sol_state[:, i * 3],
                    grav_plot.simulator.sol_state[:, 1 + i * 3],
                )
                # Plot the last position as a filled circle
                ax.plot(
                    grav_plot.simulator.sol_state[-1, i * 3],
                    grav_plot.simulator.sol_state[-1, 1 + i * 3],
                    "o",
                    color=traj[0].get_color(),
                )

        if grav_plot.system in grav_plot.solar_like_systems:
            fig.legend(loc="center right", borderaxespad=0.2)
            fig.tight_layout()

        plt.show()
        print()

    @staticmethod
    def plot_3d_trajectory(grav_plot):
        print("Plotting 3D trajectory...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_box_aspect([1.0, 1.0, 1.0])  
        ax.set_xlabel("$x$ (AU)")
        ax.set_ylabel("$y$ (AU)")
        ax.set_zlabel("$z$ (AU)")

        # Get specific colors if the system is solar-like
        match grav_plot.system:
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

        if grav_plot.system in grav_plot.solar_like_systems:
            # Plot with solar system colors
            for i in range(grav_plot.simulator.objects_count):
                traj = ax.plot(
                    grav_plot.simulator.sol_state[:, i * 3],
                    grav_plot.simulator.sol_state[:, i * 3 + 1],
                    grav_plot.simulator.sol_state[:, i * 3 + 2],
                    color=grav_plot.solar_like_systems_colors[objs_name[i]],
                )
                # Plot the last position as a filled circle
                ax.plot(
                    grav_plot.simulator.sol_state[-1, i * 3],
                    grav_plot.simulator.sol_state[-1, i * 3 + 1],
                    grav_plot.simulator.sol_state[-1, i * 3 + 2],
                    "o",
                    color=traj[0].get_color(),
                    label=objs_name[i],
                )
        else:
            # Plot with random colors
            for i in range(grav_plot.simulator.objects_count):
                traj = ax.plot(
                    grav_plot.simulator.sol_state[:, i * 3],
                    grav_plot.simulator.sol_state[:, i * 3 + 1],
                    grav_plot.simulator.sol_state[:, i * 3 + 2],
                )
                # Plot the last position as a filled circle
                ax.plot(
                    grav_plot.simulator.sol_state[-1, i * 3],
                    grav_plot.simulator.sol_state[-1, i * 3 + 1],
                    grav_plot.simulator.sol_state[-1, i * 3 + 2],
                    "o",
                    color=traj[0].get_color(),
                )
        
        if grav_plot.system in grav_plot.solar_like_systems:
            fig.legend(loc="center right", borderaxespad=0.2)
            fig.tight_layout()

        Plotter.set_3d_axes_equal(ax)
        plt.show()
        print()
    
    @staticmethod
    def set_3d_axes_equal(ax):
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
        plot_radius = 0.5 * max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    @staticmethod
    def ask_user_input_animation(grav_plot):
        while True:
            while True:
                try:
                    fps = int(input("Enter FPS (int): "))
                    if fps <= 0:
                        raise ValueError
                    print()
                    break
                except ValueError:
                    print("Invalid input! Please try again.")
                    print()

            print(f"There are {grav_plot.data_size} lines of data.")
            print(f"For FPS = {fps}, the gif would last for about {grav_plot.data_size / fps:.2f} s if plot every nth point = 1.")

            while True:
                try:
                    plot_every_nth_point = int(input("Plot every nth point (int): "))
                    if plot_every_nth_point <= 0 or plot_every_nth_point > grav_plot.data_size:
                        raise ValueError
                    print()
                    break
                except ValueError:
                    print("Invalid input! Please try again.")
                    print()

            total_frame_size = math.floor(grav_plot.data_size / plot_every_nth_point) + 1
            
            while True:
                try:
                    file_name = input("Enter file name without extension (carefully, the program cannot check the validity of the filename): ")
                    print()
                    break
                except ValueError:
                    print("Invalid input! Please try again.")
                    print()
            
            while True:
                try:
                    dpi = int(input("Enter dots per inch (dpi) in int (recommended value is 200): "))
                    if dpi <= 0:
                        raise ValueError
                    print()
                    break
                except ValueError:
                    print("Invalid input! Please try again.")
                    print()

            print(f"FPS = {fps}")
            print(f"Plot every nth point: {plot_every_nth_point}")
            print(f"Estimated total frame: about {total_frame_size}")
            print(f"Estimated time length: {(total_frame_size / fps):.2f} s")
            print(f"File name: {file_name}")
            print(f"dpi: {dpi}")

            if grav_plot.ask_user_permission("Proceed?"):
                print()
                break
            else:
                print()

        return fps, plot_every_nth_point, file_name, dpi
    
    @staticmethod 
    def animation_2d_traj_gif(grav_plot, fps, plot_every_nth_point, file_name, dpi):
        print("Animating 2D trajectory (xy plane) in .gif...")
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect="equal")

        # Get specific colors if the system is solar-like
        match grav_plot.system:
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

        file_path = Path(__file__).parent / "results"
        file_path.mkdir(parents=True, exist_ok=True)
        file_path /= (file_name + ".gif")

        writer = PillowWriter(fps=fps)
        with writer.saving(fig, file_path, dpi):
            # Plot once every nth point
            for i in range(grav_plot.data_size):
                if i % plot_every_nth_point != 0 and i != (grav_plot.data_size - 1):
                    continue

                if grav_plot.system in grav_plot.solar_like_systems:
                    # Plot with solar system colors
                    # Plot the trajectory from the beginning to current position
                    for j in range(grav_plot.simulator.objects_count):
                        traj = ax.plot(
                            grav_plot.simulator.sol_state[:(i + 1), j * 3],
                            grav_plot.simulator.sol_state[:(i + 1), j * 3 + 1],
                            color=grav_plot.solar_like_systems_colors[objs_name[j]],
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            grav_plot.simulator.sol_state[i, j * 3],
                            grav_plot.simulator.sol_state[i, j * 3 + 1],
                            "o",
                            color=traj[0].get_color(),
                            label=objs_name[j],
                        )
                    # Add legend at the beginning
                    if i == 0:
                        fig.legend(loc="center right", borderaxespad=0.2)
                else:
                    # Plot with random colors
                    # Plot the trajectory from the beginning to current position
                    for j in range(grav_plot.simulator.objects_count):
                        traj = ax.plot(
                            grav_plot.simulator.sol_state[:(i + 1), j * 3],
                            grav_plot.simulator.sol_state[:(i + 1), j * 3 + 1],
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            grav_plot.simulator.sol_state[i, j * 3],
                            grav_plot.simulator.sol_state[i, j * 3 + 1],
                            "o",
                            color=traj[0].get_color(),
                        )

                # Set axis labels before capturing the frame
                ax.set_xlabel("$x$ (AU)")
                ax.set_ylabel("$y$ (AU)")

                # Capture the frame
                writer.grab_frame()

                # Clear the plot to prepare for the next frame
                ax.clear()

        print(f"Output completed! Please check {file_path}")
        print()

    @staticmethod
    def animation_3d_traj_gif(grav_plot, fps, plot_every_nth_point, file_name, dpi):
        print("Animating 3D trajectory in .gif...")
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_box_aspect([1.0, 1.0, 1.0])  

        # Get specific colors if the system is solar-like
        match grav_plot.system:
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

        file_path = Path(__file__).parent / "results"
        file_path.mkdir(parents=True, exist_ok=True)
        file_path /= (file_name + ".gif")

        writer = PillowWriter(fps=fps)
        with writer.saving(fig, file_path, dpi):
            # Plot once every nth point
            for i in range(grav_plot.data_size):
                if i % plot_every_nth_point != 0 and i != (grav_plot.data_size - 1):
                    continue

                if grav_plot.system in grav_plot.solar_like_systems:
                    # Plot with solar system colors
                    # Plot the trajectory from the beginning to current position
                    for j in range(grav_plot.simulator.objects_count):
                        traj = ax.plot(
                            grav_plot.simulator.sol_state[:(i + 1), j * 3],
                            grav_plot.simulator.sol_state[:(i + 1), j * 3 + 1],
                            grav_plot.simulator.sol_state[:(i + 1), j * 3 + 2],
                            color=grav_plot.solar_like_systems_colors[objs_name[j]],
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            grav_plot.simulator.sol_state[i, j * 3],
                            grav_plot.simulator.sol_state[i, j * 3 + 1],
                            grav_plot.simulator.sol_state[i, j * 3 + 2],
                            "o",
                            color=traj[0].get_color(),
                            label=objs_name[j],
                        )
                    # Add legend at the beginning
                    if i == 0:
                        fig.legend(loc="center right") 
                else:
                    # Plot with random colors
                    # Plot the trajectory from the beginning to current position
                    for j in range(grav_plot.simulator.objects_count):
                        traj = ax.plot(
                            grav_plot.simulator.sol_state[:(i + 1), j * 3],
                            grav_plot.simulator.sol_state[:(i + 1), j * 3 + 1],
                            grav_plot.simulator.sol_state[:(i + 1), j * 3 + 2],
                        )
                        # Plot the current position as a filled circle
                        ax.plot(
                            grav_plot.simulator.sol_state[i, j * 3],
                            grav_plot.simulator.sol_state[i, j * 3 + 1],
                            grav_plot.simulator.sol_state[i, j * 3 + 2],
                            "o",
                            color=traj[0].get_color(),
                        )

                # Set axis labels and setting 3d axes scale before capturing the frame
                ax.set_xlabel("$x$ (AU)")
                ax.set_ylabel("$y$ (AU)")
                ax.set_zlabel("$z$ (AU)")
                Plotter.set_3d_axes_equal(ax)

                # Capture the frame
                writer.grab_frame()

                # Clear the plot to prepare for the next frame
                ax.clear()

        print(f"Output completed! Please check {file_path}")
        print()
        
    @staticmethod      
    def plot_rel_energy(grav_plot):
        if not grav_plot.computed_energy:
            if grav_plot.ask_user_permission(
                "WARNING: Energy has not been computed. Compute energy?"
            ):
                grav_plot.simulator.compute_energy()
                grav_plot.computed_energy = True

        print("Plotting relative energy error...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.semilogy(
            grav_plot.sol_time_in_tf_unit,
            np.abs(
                (grav_plot.simulator.energy - grav_plot.simulator.energy[0])
                / grav_plot.simulator.energy[0]
            ),
        )
        ax.set_title("Relative energy error against time")
        ax.set_xlabel(f"Time ({grav_plot.tf_unit})")
        ax.set_ylabel("$|(E(t)-E_0)/E_0|$")

        plt.show()
        print()

    @staticmethod
    def plot_tot_energy(grav_plot):
        # WARNING: The unit is in solar masses, AU and day
        print("Plotting total energy...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.semilogy(grav_plot.sol_time_in_tf_unit, np.abs(grav_plot.simulator.energy))
        ax.set_title("Total energy against time")
        ax.set_xlabel(f"Time ({grav_plot.tf_unit})")
        ax.set_ylabel("$E(t)$")

        plt.show()
        print()

    @staticmethod
    def plot_dt(grav_plot):
        """
        Plot dt(days)
        """
        print("Plotting dt...(Please check the window)")
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(grav_plot.sol_time_in_tf_unit, grav_plot.simulator.sol_dt, s=0.1)
        ax.set_yscale("log")
        ax.set_title("dt against time")
        ax.set_xlabel(f"Time ({grav_plot.tf_unit})")
        ax.set_ylabel("dt(days)")

        plt.show()
        print()

    @staticmethod
    def plot_compare_rel_energy(grav_plot):
        if not grav_plot.computed_energy:
            if grav_plot.ask_user_permission(
                "WARNING: Energy has not been computed. Compute energy?"
            ):
                grav_plot.simulator.compute_energy()
                grav_plot.computed_energy = True

        filepaths_to_labels = {}
        filepaths_to_labels["current"] = input("Enter plot label for current set of data in memory: ")

        read_folder_path = Path(__file__).parent / "results"
        
        while True:
            num_of_data = len(filepaths_to_labels)
            if num_of_data >= 2:
                if not grav_plot.ask_user_permission(f"There are {num_of_data} sets of data. Continue adding new files?"):
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
        ax.set_xlabel(f"Time ({grav_plot.tf_unit})")
        ax.set_ylabel("$|(E(t)-E_0)/E_0|$")

        for filepath in filepaths_to_labels:
            if filepath == "current":
                ax.semilogy(
                    grav_plot.sol_time_in_tf_unit,
                    np.abs((grav_plot.simulator.energy - grav_plot.simulator.energy[0]) / grav_plot.simulator.energy[0]),
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

                        if grav_plot.tf_unit == "years":
                            sol_time /= grav_plot.SIDEREAL_DAYS_PER_YEAR

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
