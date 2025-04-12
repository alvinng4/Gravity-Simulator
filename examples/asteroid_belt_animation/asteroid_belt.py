"""
Demonstration on using the gravity simulator to simulate the asteroid belt.

Note:   1. You will need two more packages to run this script. The first one is Pillow,
           which is used to create the gif. The second one is rich, which is used to show
           the progress bar. You can avoid using rich simply by removing the "track" function
           in the for loops.
        2. Technically you can also create nice looking solar system animations by setting N = 0 and
           expanding the axes limits.
        3. We are using "massless" acceleration method because we are considering
           a short time frame, so the gravitational effects from the asteroids are negligible.
           (Even if you were to include the gravitational effect from the asteroids, it would be
           smaller than the round-off error, so it is not worth the effort.)

Warning: Do not run multiple instances of this program at the same time, unless you made copies
         of the whole directory. Otherwise, the final data may overwrite each other.

Author: Ching-Yin Ng
"""

import glob
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import PIL
from grav_sim import GravitySimulatorAPI
from rich.progress import track


# Option 0: asteroid belt
# Option 1: add a star near the asteroid belt
# Option 2: same as option 1 but different position
OPTION = 1

N = 50000
FPS = 30
DPI = 200
N_FRAMES = 500

TF = 5.0 # Number of years
DT = 0.1 # Time step in days

ACC_METHOD = "massless"
INTEGRATOR = "rk4"

OUTPUT_INITIAL = False
SNAPSHOTS_DIR = Path("snapshots/")
SNAPSHOTS_DIR.mkdir(parents=True, exist_ok=True)

FRAMES_DIR = Path("frames/")
FRAMES_DIR.mkdir(parents=True, exist_ok=True)

DELETE_SNAPSHOTS = True # Delete snapshots after creating the gif

def main():
    # ---------- Initialization ---------- #
    print("Initializing the system...", end="")
    
    gs = GravitySimulatorAPI()
    
    system = gs.get_built_in_system("solar_system")
    system.remove([1, 7, 8])  # Remove Mercury, Uranus, and Neptune
    objects_name = ["Sun", "Venus", "Earth", "Mars", "Jupiter", "Saturn"]
    colors = [gs.SOLAR_SYSTEM_COLORS[name] for name in objects_name]
    marker_sizes = [6.0, 1.5, 2.0, 1.5, 4.0, 3.5]

    #################################################
    # Adding a star to the system

    if OPTION == 1:
        system.add_keplerian(
            semi_major_axis=5.5,
            eccentricity=0.7,
            inclination=0.05,
            argument_of_periapsis=0.07,
            longitude_of_ascending_node=0.07,
            true_anomaly=0.35,
            m=1.0,
            primary_particle_id=0,
        )
        objects_name.append("Star 1")
        colors.append("orange")
        marker_sizes.append(6.0)

    elif OPTION == 2:
        system.add_keplerian(
            semi_major_axis=5.0,
            eccentricity=0.7,
            inclination=0.4,
            argument_of_periapsis=4.0,
            longitude_of_ascending_node=4.0,
            true_anomaly=4.0,
            m=1.0,
            primary_particle_id=0,
        )
        objects_name.append("Star 1")
        colors.append("orange")
        marker_sizes.append(6.0)

    #################################################

    massive_objects_count = system.num_particles

    rng = np.random.default_rng()
    a = rng.uniform(2.1, 3.2, size=N)
    ecc = np.abs(rng.normal(loc=0.0, scale=0.12, size=N))
    inc = np.abs(rng.normal(loc=0.0, scale=0.3, size=N))
    argument_of_periapsis = rng.uniform(0, 2 * np.pi, size=N)
    long_asc_node = rng.uniform(0, 2 * np.pi, size=N)
    true_anomaly = rng.uniform(0, 2 * np.pi, size=N)

    # m = 0.0 if we assume asteroids are massless
    for i in range(N):
        system.add_keplerian(
            semi_major_axis=a[i],
            eccentricity=ecc[i],
            inclination=inc[i],
            argument_of_periapsis=argument_of_periapsis[i],
            longitude_of_ascending_node=long_asc_node[i],
            true_anomaly=true_anomaly[i],
            m=0.0,
            primary_particle_id=0,
        )

    system.center_of_mass_correction()
    print("Done!")

    # ---------- Simulation ---------- #
    tf = gs.years_to_days(TF)

    acc_param, integrator_param, output_param, settings = gs.get_new_parameters()
    acc_param.method = ACC_METHOD

    integrator_param.integrator = INTEGRATOR
    integrator_param.dt = 0.1

    output_param.method = "csv"
    output_param.output_interval = tf / N_FRAMES
    output_param.output_initial = OUTPUT_INITIAL
    output_param.output_dir = SNAPSHOTS_DIR

    gs.launch_simulation(system, acc_param, integrator_param, output_param, settings, tf)

    # ---------- Animation ---------- #
    print("Drawing frames...", end="")

    snapshot_files = sorted(glob.glob(str(SNAPSHOTS_DIR / "snapshot_*.csv")))
    actual_n_frames = len(snapshot_files)
    if actual_n_frames == 0:
        raise ValueError(f"No snapshot files found at {SNAPSHOTS_DIR}.")
    else:
        print(f"Found {actual_n_frames} snapshot files.")

    # We can't use the library function to read all
    # data in the memory as it is too large.
    # We have to read them manually one by one.
    for i in track(range(actual_n_frames)):
        snapshot_file = snapshot_files[i]
        data = np.genfromtxt(snapshot_file, delimiter=",", skip_header=5)
        x = data[:, 2:5]

        fig = plt.figure()
        plt.style.use("dark_background")
        ax = fig.add_subplot(111, projection="3d")

        xlim_min = -3
        xlim_max = 3
        ylim_min = -3
        ylim_max = 3
        zlim_min = -3
        zlim_max = 3

        ax.grid(False)
        ax.xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        ax.yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        ax.zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.zaxis.set_visible(False)
        ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

        # Plotting massive objects
        for j in range(massive_objects_count):
            ax.scatter(
                x[j, 0],
                x[j, 1],
                x[j, 2],
                marker="o",
                label=objects_name[j],
                color=colors[j],
                s=marker_sizes[j],
            )

        # Plotting massless objects
        ax.scatter(
            x[massive_objects_count:, 0],
            x[massive_objects_count:, 1],
            x[massive_objects_count:, 2],
            color="white",
            marker=".",
            s=0.1,
            alpha=0.2,
        )

        # Add legend
        legend = ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
        legend.facecolor = "transparent"

        # Adjust figure for the legend
        fig.subplots_adjust(right=0.7)
        fig.tight_layout()

        # Set axis labels and setting 3d axes scale before capturing the frame
        # ax.set_xlabel("$x$ (AU)")
        # ax.set_ylabel("$y$ (AU)")
        # ax.set_zlabel("$z$ (AU)")

        ax.set_xlim3d([xlim_min, xlim_max])
        ax.set_ylim3d([ylim_min, ylim_max])
        ax.set_zlim3d([zlim_min, zlim_max])

        # Set equal aspect ratio to prevent distortion
        ax.set_aspect("equal")

        # Capture the frame
        plt.savefig(FRAMES_DIR / f"frames_{i:04d}.png", dpi=DPI)
        plt.close("all")
    plt.close("all")

    print()
    print("Combining frames to gif...")

    def frames_generator():
        for i in range(actual_n_frames):
            yield PIL.Image.open(FRAMES_DIR / f"frames_{i:04d}.png")

    frames = frames_generator()
    next(frames).save(
        FRAMES_DIR / "asteroid_belt.gif",
        save_all=True,
        append_images=frames,
        loop=0,
        duration=(1000 // FPS),
    )

    for i in range(actual_n_frames):
        (FRAMES_DIR / f"frames_{i:04d}.png").unlink()

    print(f"Output completed! Please check {FRAMES_DIR / 'asteroid_belt.gif'}")
    print()
    
    if DELETE_SNAPSHOTS:
        gs.delete_snapshots(output_dir=SNAPSHOTS_DIR)


if __name__ == "__main__":
    main()
