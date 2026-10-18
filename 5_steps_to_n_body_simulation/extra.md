This is an extra component to 5 steps to n-body simulation. I figured some of you may be interested in the visuals,
so in this section, I will show you how to make 3D plots and produce animations using the `matplotlib` library.

## Plotting in 3D
It would be nice to plot the trajectory in 3D. To do this, we will need two functions,
one to set the 3D axes limits in equal aspect ratio, and one to plot the trajectory in 3D.

??? Note "Code (Click to expand)"
    ```python title="common.py"
    def set_3d_axes_equal(ax: plt.Axes) -> None:
        """
        Make axes of 3D plot have equal scale

        Parameters
        ----------
        ax : matplotlib axis
            The axis to set equal scale

        Reference
        ---------
        karlo, https://stackoverflow.com/questions/13685386/how-to-set-the-equal-aspect-ratio-for-all-axes-x-y-z
        """

        x_limits = ax.get_xlim3d()  # type: ignore
        y_limits = ax.get_ylim3d()  # type: ignore
        z_limits = ax.get_zlim3d()  # type: ignore

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5 * max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])  # type: ignore
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])  # type: ignore
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])  # type: ignore


    def plot_3d_trajectory(
        sol_x: np.ndarray,
        labels: list,
        colors: list,
        legend: bool,
    ) -> None:
        """
        Plot the 3D trajectory.

        Parameters
        ----------
        sol_x : np.ndarray
            Solution position array with shape (N_steps, num_particles, 3).
        labels : list
            List of labels for the particles.
        colors : list
            List of colors for the particles.
        legend : bool
            Whether to show the legend.
        """

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("$x$ (AU)")
        ax.set_ylabel("$y$ (AU)")
        ax.set_zlabel("$z$ (AU)")  # type: ignore

        for i in range(sol_x.shape[1]):
            traj = ax.plot(
                sol_x[:, i, 0],
                sol_x[:, i, 1],
                sol_x[:, i, 2],
                color=colors[i],
            )
            # Plot the last position with marker
            ax.scatter(
                sol_x[-1, i, 0],
                sol_x[-1, i, 1],
                sol_x[-1, i, 2],
                marker="o",
                color=traj[0].get_color(),
                label=labels[i],
            )

        set_3d_axes_equal(ax)

        if legend:
            ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
            fig.subplots_adjust(right=0.7)

        plt.show()
    ```

To make the plot look better, I added a few dwarf planets like Pluto. The initial condition
is already available in the `get_initial_conditions` function (input `solar_system_plus`).
In step 4 or step 5, change `plot_2d_trajectory` to `plot_3d_trajectory`, and you should see a 3D plot:

<img src="../../examples/media/solar_plus_3d.png" alt="3D Trajectory" width="400"/>

## Animations

Animations generally consists of two steps:

1. Generating the frames
2. Combining the frames into a video

If you know how to make a plot, you already know how to generate frames.
First, let us copy `step5.py` to `extra.py` and modify it slightly to include
`solar_system_plus`.

??? Note "Code (Click to expand)"
    ```python title="extra.py"
    OPTION = 2

    # Default units is AU, days, and M_sun

    # Solar system
    if OPTION == 0:
        INITIAL_CONDITION = "solar_system"
        TF = 200.0 * 365.24  # 200 years to days
        TOLERANCE = 1e-8
        OUTPUT_INTERVAL = 0.01 * 365.24  # 0.01 year to days
        INITIAL_DT = 1.0

    # Pyth-3-body
    elif OPTION == 1:
        INITIAL_CONDITION = "pyth-3-body"
        TF = 70.0
        TOLERANCE = 1e-13
        OUTPUT_INTERVAL = 0.001
        INITIAL_DT = 0.01

    elif OPTION == 2:
        INITIAL_CONDITION = "solar_system_plus"
        TF = 250.0 * 365.24  # 200 years to days
        TOLERANCE = 1e-8
        OUTPUT_INTERVAL = 0.5 * 365.24  # 0.5 year to days
        INITIAL_DT = 1.0

    else:
        raise ValueError(
            "Invalid option. Choose 0 for solar system, 1 for Pyth-3-body, or 2 for solar system plus."
        )
    ```

!!! Warning "Output interval"
    Be careful with the output interval. In this section, we are
    converting all solution output into frames. A few hundred frames should be enough.

We need a place to store the frames. Lets store it at
`5_steps_to_n_body_simulation/figures/frames` (or anywhere you like).

```python title="extra.py"
from pathlib import Path

FRAMES_DIR = Path(__file__).parent.parent / "figures" / "frames"
FRAMES_DIR.mkdir(parents=True, exist_ok=True)
```

For simplicity, let us draw the frames in the `main` function.
Note that we need to set the min and max values for the axes.
```python title="extra.py"
    print("Drawing frames...")
    x_min, x_max = np.min(sol_x[:, :, 0]), np.max(sol_x[:, :, 0])
    y_min, y_max = np.min(sol_x[:, :, 1]), np.max(sol_x[:, :, 1])
    z_min, z_max = np.min(sol_x[:, :, 2]), np.max(sol_x[:, :, 2])
    xyz_min = np.min([x_min, y_min, z_min])
    xyz_max = np.max([x_max, y_max, z_max])

    for n in range(output_count):
        print(f"Progress: {n + 1} / {output_count}", end="\r")

        # Draw the trajectory
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("$x$ (AU)")
        ax.set_ylabel("$y$ (AU)")
        ax.set_zlabel("$z$ (AU)")  # type: ignore

        for i in range(sol_x.shape[1]):
            traj = ax.plot(
                sol_x[:n, i, 0],
                sol_x[:n, i, 1],
                sol_x[:n, i, 2],
                color=colors[i],
            )
            # Plot the last position with marker
            ax.scatter(
                sol_x[n, i, 0],
                sol_x[n, i, 1],
                sol_x[n, i, 2],
                marker="o",
                color=traj[0].get_color(),
                label=labels[i],
            )

        ax.set_xlim3d(xyz_min, xyz_max)
        ax.set_ylim3d(xyz_min, xyz_max)
        ax.set_zlim3d(xyz_min, xyz_max)

        # Set equal aspect ratio to prevent distortion
        ax.set_aspect("equal")

        if legend:
            ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
            fig.subplots_adjust(right=0.7)
            fig.tight_layout()

        plt.savefig(FRAMES_DIR / f"frames_{n:05d}.png")
        plt.close("all")
    print("\nDone!")
```

Now, we can combine the frames into a gif file using the `PIL` library.
(If you want a mp4 file instead, you can use the `ffmpeg` library.)
You can install `PIL` using 

```
pip install pillow
```

The code is given below. The logic is simple: we open the first frame and append the rest of 
the frames to it. After the animation is done, we delete all frames by using `Path.unlink()`.

!!! Tip "Generator function"
    We used a generator function `frames_generator`. If you use a list instead,
    you may run out of memory very quickly.

```python title="extra.py"
    print("Combining frames to gif...")

    def frames_generator():
        for i in range(output_count):
            yield PIL.Image.open(FRAMES_DIR / f"frames_{i:05d}.png")

    fps = 30
    frames = frames_generator()
    next(frames).save(
        FRAMES_DIR / "animation.gif",
        save_all=True,
        append_images=frames,
        loop=0,
        duration=(1000 // fps),
    )

    for i in range(output_count):
        (FRAMES_DIR / f"frames_{i:05d}.png").unlink()

    print(f"Output completed! Please check {FRAMES_DIR / 'animation.gif'}")
    print()
```

You should see the animation:

<iframe width="560" height="315" src="https://www.youtube.com/embed/Jqjbyi8ZXvk?si=OolUY32zeqB80Tu3" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

## Conclusion

In this section, we learned how to make 3D plots and animations using the `matplotlib` library.

## Full scripts
The full scripts are available at `5_steps_to_n_body_simulation/python/`,
or https://github.com/alvinng4/grav_sim/blob/main/5_steps_to_n_body_simulation/python/

??? note "extra.py (Click to expand)"
    ```python linenums="1" title="5_steps_to_n_body_simulation/python/extra.py"
    --8<-- "5_steps_to_n_body_simulation/python/extra.py"
    ```

??? note "common.py (Click to expand)"
    ```python linenums="1" title="5_steps_to_n_body_simulation/python/common.py"
    --8<-- "5_steps_to_n_body_simulation/python/common.py"
    ```
