"""
Gravitational system
"""

import ctypes
import csv
from pathlib import Path
from typing import Optional

import numpy as np

from . import plotting
from . import utils


class System:
    """
    Class to represent a gravitational system.
    """

    # Built-in systems
    BUILT_IN_SYSTEMS = [
        "circular_binary_orbit",
        "eccentric_binary_orbit",
        "3d_helix",
        "sun_earth_moon",
        "figure-8",
        "pyth-3-body",
        "solar_system",
        "solar_system_plus",
    ]

    # Gravitational constant (AU^3 day^-2 M_sun^-1)
    G_DEFAULT = 0.000295912208284119496676630

    def __init__(self, c_lib: ctypes.CDLL) -> None:
        """Initialize System object"""
        self.c_lib = c_lib
        self.num_particles = 0
        self.particle_ids = np.zeros((0,), dtype=np.int32)
        self.x = np.zeros((0, 3), dtype=np.float64)
        self.v = np.zeros((0, 3), dtype=np.float64)
        self.m = np.zeros((0,), dtype=np.float64)
        self.G = self.G_DEFAULT

    def add(
        self,
        x: list | np.ndarray,
        v: list | np.ndarray,
        m: float | list | np.ndarray,
    ) -> None:
        """Add one or multiple particle to the system

        Parameters
        ----------
        x : list or np.ndarray
            Position vector(s) of the object(s)
        v : list or np.ndarray
            Velocity vector(s) of the object(s)
        m : float
            Mass(es) of the object(s)
        """
        self.x = np.vstack((self.x, np.array(x, dtype=np.float64)))
        self.v = np.vstack((self.v, np.array(v, dtype=np.float64)))
        self.m = np.hstack((self.m, m))

        if isinstance(m, (list, np.ndarray)):
            num_new_particles = len(m)
        else:
            num_new_particles = 1

        max_id = self.particle_ids.max() if self.num_particles > 0 else -1
        new_ids = np.arange(max_id + 1, max_id + 1 + num_new_particles, dtype=np.int32)
        self.particle_ids = np.hstack((self.particle_ids, new_ids))
        self.num_particles += num_new_particles

    def add_keplerian(
        self,
        semi_major_axis: float,
        eccentricity: float,
        inclination: float,
        argument_of_periapsis: float,
        longitude_of_ascending_node: float,
        true_anomaly: float,
        m: float,
        primary_particle_id: int,
    ):
        """
        Add a celestial body to the system using Keplerian elements

        Warning: This method use the G value from the system. Make sure
                 to set the correct G value before using this method.

        Parameters
        ----------
        semi_major_axis : float
            Semi-major axis
        eccentricity : float
            Eccentricity
        inclination : float
            Inclination
        argument_of_periapsis : float
            Argument of periapsis
        longitude_of_ascending_node : float
            Longitude of ascending node
        true_anomaly : float
            True anomaly
        m : float
            Mass
        primary_object_index : int
            Particle id of the primary object

        Raises
        ------
        ValueError
            If particle_id is not found in the system
        """

        primary_particle_idx = np.where(self.particle_ids == primary_particle_id)[0]
        if len(primary_particle_idx) == 0:
            raise ValueError(
                f"Particle ID {primary_particle_id} not found in the system"
            )
        primary_particle_idx = primary_particle_idx[0]

        x, v = utils.keplerian_to_cartesian(
            self.c_lib,
            semi_major_axis=semi_major_axis,
            eccentricity=eccentricity,
            inclination=inclination,
            argument_of_periapsis=argument_of_periapsis,
            longitude_of_ascending_node=longitude_of_ascending_node,
            true_anomaly=true_anomaly,
            total_mass=float(self.m[primary_particle_idx] + m),
            G=self.G,
        )
        self.add(x + self.x[primary_particle_idx], v + self.v[primary_particle_idx], m)

    def remove(
        self,
        particle_ids: int | list | np.ndarray,
    ) -> None:
        """Remove particle(s) from the system

        Parameters
        ----------
        particle_ids : int, list or np.ndarray
            particle_id(s) to remove

        Raises
        ------
        TypeError
            If indices is not an int, list or np.ndarray
        """
        if isinstance(particle_ids, int):
            particle_ids_arr = np.array([particle_ids])
        elif isinstance(particle_ids, list):
            particle_ids_arr = np.array(particle_ids)
        elif not isinstance(particle_ids, np.ndarray):
            raise TypeError("indices must be an int, list or np.ndarray")
        else:
            particle_ids_arr = particle_ids

        mask = np.isin(self.particle_ids, particle_ids_arr, invert=True)

        self.particle_ids = self.particle_ids[mask]
        self.x = self.x[mask]
        self.v = self.v[mask]
        self.m = self.m[mask]
        self.num_particles = len(self.particle_ids)

    def save(
        self,
        file_path: Optional[str | Path] = None,
    ) -> None:
        """Save system to a CSV file

        Parameters
        ----------
        file_path : str, optional
            File path to save the system, by default None

        Note
        ----
        This function has the following side effects:
        - Prints a message f"System \"{self.name}\" successfully saved to \"{file_path}\""
        """
        if file_path is None:
            file_path = Path(__file__).parent / "customized_systems.csv"

        with open(file_path, "w", newline="") as file:
            writer = csv.writer(file)

            # Write header
            writer.writerow(["particle_id", "m", "x", "y", "z", "vx", "vy", "vz"])

            # Write data
            for i in range(self.num_particles):
                writer.writerow(
                    [
                        self.particle_ids[i],
                        self.m[i],
                        self.x[i, 0],
                        self.x[i, 1],
                        self.x[i, 2],
                        self.v[i, 0],
                        self.v[i, 1],
                        self.v[i, 2],
                    ]
                )

        print(f'System successfully saved to "{file_path}"')

    @staticmethod
    def load_system(
        c_lib: ctypes.CDLL,
        file_path: str | Path,
    ) -> "System":
        """Load system from a CSV file

        Parameters
        ----------
        system_name : str
            Name of the system to load
        file_path : str
            File path to load the system from

        Raises
        ------
        FileNotFoundError
            If the file is not found
        """
        if file_path is not None:
            if isinstance(file_path, str):
                file_path = Path(file_path)

            if not file_path.is_file():
                raise FileNotFoundError(f"File not found: {file_path}")

        with open(file_path, "r") as file:
            reader = csv.reader(file)
            next(reader)  # Skip header

            system = System(c_lib)
            particle_ids = []
            for row in reader:
                system.add(
                    x=[float(row[2]), float(row[3]), float(row[4])],
                    v=[float(row[5]), float(row[6]), float(row[7])],
                    m=float(row[1]),
                )
                particle_ids.append(int(row[0]))
            system.particle_ids = np.array(particle_ids, dtype=np.int32)

        return system

    def center_of_mass_correction(self) -> None:
        """Set center of mass of position and V_CM to zero"""
        r_cm = np.sum(self.m[:, np.newaxis] * self.x, axis=0) / np.sum(self.m)
        v_cm = np.sum(self.m[:, np.newaxis] * self.v, axis=0) / np.sum(self.m)

        self.x -= r_cm
        self.v -= v_cm

    @staticmethod
    def get_built_in_system(c_lib: ctypes.CDLL, system_name: str) -> "System":
        """Get a built-in system from the C library

        Parameters
        ----------
        c_lib : ctypes.CDLL
            C dynamic-link library object
        system_name : str
            Name of the built-in system to be loaded.

        Returns
        -------
        System
            A System object with the loaded built-in system data.
        """
        # if system_name not in System.BUILT_IN_SYSTEMS:
        #     raise ValueError(
        #         f"Unknown system name: {system_name}. Available systems: {System.BUILT_IN_SYSTEMS}"
        #     )

        new_system = System(c_lib)

        num_particles = ctypes.c_int32()
        particle_ids = ctypes.POINTER(ctypes.c_int32)()
        x = ctypes.POINTER(ctypes.c_double)()
        v = ctypes.POINTER(ctypes.c_double)()
        m = ctypes.POINTER(ctypes.c_double)()
        G = ctypes.c_double()
        return_value = c_lib.load_built_in_system_python(
            system_name.encode("utf-8"),
            ctypes.byref(num_particles),
            ctypes.byref(particle_ids),
            ctypes.byref(x),
            ctypes.byref(v),
            ctypes.byref(m),
            ctypes.byref(G),
        )

        if return_value != 0:
            raise RuntimeError(
                f"Failed to load built-in system '{system_name}' from C library"
            )

        new_system.num_particles = num_particles.value
        new_system.G = G.value

        new_system.particle_ids = np.ctypeslib.as_array(
            particle_ids, shape=(num_particles.value,)
        ).copy()
        new_system.x = np.ctypeslib.as_array(x, shape=(num_particles.value, 3)).copy()
        new_system.v = np.ctypeslib.as_array(v, shape=(num_particles.value, 3)).copy()
        new_system.m = np.ctypeslib.as_array(m, shape=(num_particles.value,)).copy()

        c_lib.free_memory_int32(particle_ids)
        c_lib.free_memory_double(x)
        c_lib.free_memory_double(v)
        c_lib.free_memory_double(m)

        return new_system

    def plot_2d_system(
        self,
        colors: Optional[list[str]] = None,
        labels: Optional[list[str]] = None,
        legend: bool = False,
        xlabel: str = "$x$",
        ylabel: str = "$y$",
        title: Optional[str] = None,
        marker: str = "o",
        markersize: int = 6,
        save_fig: bool = False,
        save_fig_path: Optional[str | Path] = None,
    ) -> None:
        initial_state = np.concatenate(
            [
                self.x.flatten(),
                self.v.flatten(),
            ]
        )[np.newaxis, :]

        plotting.plot_2d_trajectory(
            initial_state,
            colors,
            labels,
            legend,
            xlabel,
            ylabel,
            title,
            marker,
            markersize,
            save_fig,
            save_fig_path,
        )

    def plot_3d_system(
        self,
        colors: Optional[list[str]] = None,
        labels: Optional[list[str]] = None,
        legend: bool = False,
        xlabel: str = "$x$",
        ylabel: str = "$y$",
        zlabel: str = "$z$",
        title: Optional[str] = None,
        marker: str = "o",
        markersize: int = 6,
        save_fig: bool = False,
        save_fig_path: Optional[str | Path] = None,
    ) -> None:
        initial_state = np.concatenate(
            [
                self.x.flatten(),
                self.v.flatten(),
            ]
        )[np.newaxis, :]

        plotting.plot_3d_trajectory(
            initial_state,
            colors,
            labels,
            legend,
            xlabel,
            ylabel,
            zlabel,
            title,
            marker,
            markersize,
            save_fig,
            save_fig_path,
        )
