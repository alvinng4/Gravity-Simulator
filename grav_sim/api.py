"""
Application user interface (API) for gravity simulator

Usage:
    from grav_sim import GravitySimulatorAPI
    gs = GravitySimulatorAPI()
"""

import copy
import ctypes
import sys
import warnings
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np

from . import plotting
from . import utils
from . import parameters
from .simulator import Simulator
from .system import System


class GravitySimulatorAPI:
    """Gravity simulator API"""

    def __init__(self, c_lib_path: Optional[str] = None) -> None:
        """
        Initialize gravity simulator API

        Parameters
        ----------
        c_lib_path : str, optional
            Path to C library, by default None
        """
        if c_lib_path is not None:
            self.c_lib = utils.load_c_lib(Path(c_lib_path))
        else:
            self.c_lib = utils.load_c_lib()
        utils.initialize_c_lib(self.c_lib)

        # System
        self.BUILT_IN_SYSTEMS = System.BUILT_IN_SYSTEMS

        # Plotting
        self.SOLAR_SYSTEM_COLORS = plotting.SOLAR_SYSTEM_COLORS
        self.plot_quantity_against_time = plotting.plot_quantity_against_time
        self.plot_2d_trajectory = plotting.plot_2d_trajectory
        self.plot_3d_trajectory = plotting.plot_3d_trajectory

        # Simulator
        self.simulator = Simulator(c_lib=self.c_lib)
        self.DAYS_PER_YEAR = self.simulator.DAYS_PER_YEAR
        self.launch_simulation = self.simulator.launch_simulation

        # Parameters
        self.AVAILABLE_ACCELERATION_METHODS = parameters.AccelerationParam.AVAILABLE_ACCELERATION_METHODS
        self.AVAILABLE_INTEGRATORS = parameters.IntegratorParam.AVAILABLE_INTEGRATORS
        self.FIXED_STEP_SIZE_INTEGRATORS = parameters.IntegratorParam.FIXED_STEP_SIZE_INTEGRATORS
        self.ADAPTIVE_STEP_SIZE_INTEGRATORS = parameters.IntegratorParam.ADAPTIVE_STEP_SIZE_INTEGRATORS
        self.AVAILABLE_OUTPUT_METHODS = parameters.OutputParam.AVAILABLE_OUTPUT_METHODS
        self.AVAILABLE_OUTPUT_DTYPE = parameters.OutputParam.AVAILABLE_OUTPUT_DTYPE

    def get_new_system(self) -> System:
        """Create a gravitational system

        Returns
        -------
        System object
        """
        return System(c_lib=self.c_lib)

    def load_system(
        self,
        file_path: str | Path,
    ) -> System:
        """Load system from a CSV file

        Parameters
        ----------
        file_path : str
            File path to load the system from
        """
        return System.load_system(self.c_lib, file_path)

    def get_built_in_system(self, system_name: str) -> System:
        """Get a built-in gravitational system

        Parameters
        ----------
        system_name : str
            Name of the built-in system to be loaded.
        """
        return System.get_built_in_system(self.c_lib, system_name)

    @staticmethod
    def get_new_parameters() -> Tuple[
        parameters.AccelerationParam,
        parameters.IntegratorParam,
        parameters.OutputParam,
        parameters.Settings,
    ]:
        """Create new simulation parameters

        Returns
        -------
        Tuple of acceleration, integrator, output, and settings parameters
        """
        acceleration_param = parameters.AccelerationParam()
        integrator_param = parameters.IntegratorParam()
        output_param = parameters.OutputParam()
        settings = parameters.Settings()

        return acceleration_param, integrator_param, output_param, settings

    def days_to_years(self, days: float | np.ndarray) -> float | np.ndarray:
        return days / self.simulator.DAYS_PER_YEAR

    def years_to_days(self, years: float | np.ndarray) -> float | np.ndarray:
        return years * self.simulator.DAYS_PER_YEAR

    @staticmethod
    def read_csv_data(
        output_dir: str | Path,
    ) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Read CSV snapshots from the output directory,
        assuming number of particles and particle_ids
        stays the same

        Parameters
        ----------
        output_dir : str | Path
            Output directory path

        Returns
        -------
        G : float
            Gravitational constant
        time : np.ndarray
            Simulation time of each snapshot
        dt : np.ndarray
            Time step of each snapshot
        particle_ids : np.ndarray
            1D array of Particle IDs
        sol_state : np.ndarray
            3D array of solution state for each snapshot, with shape
            (num_snapshots, num_particles, 7) being
            [m, x, y, z, vx, vy, vz]
        """
        output_dir = Path(output_dir)
        if not output_dir.is_dir():
            raise FileNotFoundError(f"Output directory not found: {output_dir}")

        snapshot_files = sorted(output_dir.glob("snapshot_*.csv"))
        if len(snapshot_files) == 0:
            raise FileNotFoundError(f"No snapshot files found in: {output_dir}")

        G = -1.0
        time = np.zeros(len(snapshot_files), dtype=np.float64)
        dt = np.zeros(len(snapshot_files), dtype=np.float64)

        for i, snapshot_file in enumerate(snapshot_files):
            # Read the metadata
            with open(snapshot_file, "r") as file:
                read_metadata_num_particles = False
                read_metadata_G = False
                read_metadata_time = False
                read_metadata_dt = False
                for line in file:
                    line = line.strip()

                    if line.startswith("#"):
                        if line.startswith("# num_particles"):
                            if i == 0:
                                num_particles = int(line.split(":")[1].strip())
                            elif num_particles != int(line.split(":")[1].strip()):
                                raise ValueError(
                                    f"Number of particles changed from {num_particles} to {int(line.split(':')[1].strip())}"
                                )
                            read_metadata_num_particles = True
                        elif line.startswith("# G"):
                            G = float(line.split(":")[1].strip())
                            read_metadata_G = True
                        elif line.startswith("# time"):
                            time[i] = float(line.split(":")[1].strip())
                            read_metadata_time = True
                        elif line.startswith("# dt"):
                            dt[i] = float(line.split(":")[1].strip())
                            read_metadata_dt = True

                    if (
                        read_metadata_num_particles
                        and read_metadata_G
                        and read_metadata_time
                        and read_metadata_dt
                    ):
                        break

        # Read the data
        particle_ids = np.zeros(num_particles, dtype=np.int32)
        sol_state = np.zeros((len(snapshot_files), num_particles, 7), dtype=np.float64)
        for i, snapshot_file in enumerate(snapshot_files):
            data = np.genfromtxt(snapshot_file, delimiter=",", skip_header=5)
            if i == 0:
                particle_ids = data[:, 0].astype(np.int32)
                particle_ids = np.sort(particle_ids)
                _, num_duplicates = np.unique(particle_ids, return_counts=True)
                if np.any(num_duplicates > 1):
                    raise ValueError(f"Particle IDs are not unique. Particle IDs: {particle_ids}")

            snapshot_particle_ids = data[:, 0].astype(np.int32)

            # Sort the data by particle IDs
            sorted_indices = np.argsort(snapshot_particle_ids)
            data = data[sorted_indices]

            # Check if the particle IDs match
            if not np.array_equal(particle_ids, snapshot_particle_ids):
                raise ValueError(
                    f"Particle IDs do not match in snapshot {i + 1}: {snapshot_file}"
                )

            # Store the data
            sol_state[i, :, :] = data[:, 1:]

        return G, time, dt, particle_ids, sol_state

    @staticmethod
    def delete_snapshots(
        output_dir: str | Path,
    ):
        """Delete all snapshots in the output directory

        Parameters
        ----------
        output_dir : str | Path
            Output directory path
        """
        output_dir = Path(output_dir)
        if not output_dir.is_dir():
            raise FileNotFoundError(f"Output directory not found: {output_dir}")

        snapshot_files = sorted(output_dir.glob("snapshot_*.csv"))
        for snapshot_file in snapshot_files:
            snapshot_file.unlink()

    def compute_energy(
        self: ctypes.CDLL, sol_state: np.ndarray, G: float
    ) -> np.ndarray:
        """Compute the total energy of the system

        Parameters
        ----------
        sol_state : np.ndarray
            3D array of solution state for each snapshot, with shape
            (num_snapshots, num_particles, 7) being
            [m, x, y, z, vx, vy, vz]
        G : float
            Gravitational constant

        Returns
        -------
        energy : np.ndarray
            1D array of total energy for each snapshot
        """
        # Check the dimension and shape of sol_state
        if len(sol_state.shape) != 3:
            raise ValueError("sol_state must be a 3D array")

        if sol_state.shape[2] != 7:
            raise ValueError(
                "sol_state must have shape (num_snapshots, num_particles, 7)"
            )

        # Compute the total energy
        energy = np.zeros(sol_state.shape[0], dtype=np.float64)
        self.c_lib.compute_energy_python(
            energy.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_double(G),
            sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int32(sol_state.shape[0]),
            ctypes.c_int32(sol_state.shape[1]),
        )

        return energy

    def compute_linear_momentum(
        self: ctypes.CDLL, sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute the total linear_momentum of the system

        Parameters
        ----------
        sol_state : np.ndarray
            3D array of solution state for each snapshot, with shape
            (num_snapshots, num_particles, 7) being
            [m, x, y, z, vx, vy, vz]

        Returns
        -------
        linear_momentum : np.ndarray
            1D array of total linear_momentum for each snapshot
        """
        # Check the dimension and shape of sol_state
        if len(sol_state.shape) != 3:
            raise ValueError("sol_state must be a 3D array")

        if sol_state.shape[2] != 7:
            raise ValueError(
                "sol_state must have shape (num_snapshots, num_particles, 7)"
            )

        # Compute the total energy
        linear_momentum = np.zeros(sol_state.shape[0], dtype=np.float64)
        self.c_lib.compute_linear_momentum_python(
            linear_momentum.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int32(sol_state.shape[0]),
            ctypes.c_int32(sol_state.shape[1]),
        )

        return linear_momentum
    
    def compute_angular_momentum(
        self: ctypes.CDLL, sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute the total angular_momentum of the system

        Parameters
        ----------
        sol_state : np.ndarray
            3D array of solution state for each snapshot, with shape
            (num_snapshots, num_particles, 7) being
            [m, x, y, z, vx, vy, vz]

        Returns
        -------
        angular_momentum : np.ndarray
            1D array of total angular_momentum for each snapshot
        """
        # Check the dimension and shape of sol_state
        if len(sol_state.shape) != 3:
            raise ValueError("sol_state must be a 3D array")

        if sol_state.shape[2] != 7:
            raise ValueError(
                "sol_state must have shape (num_snapshots, num_particles, 7)"
            )

        # Compute the total energy
        angular_momentum = np.zeros(sol_state.shape[0], dtype=np.float64)
        self.c_lib.compute_angular_momentum_python(
            angular_momentum.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            sol_state.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int32(sol_state.shape[0]),
            ctypes.c_int32(sol_state.shape[1]),
        )

        return angular_momentum
    
    @staticmethod
    def compute_eccentricity(
        G: float,
        sol_state: np.ndarray,
    ) -> np.ndarray:
        """Compute the eccentricity using the sol_state array,
        assuming that the first object is the central object

        Parameters
        ----------
        G : float
            Gravitational constant
        sol_state : np.ndarray
            Solution state of the system

        Returns
        -------
        np.ndarray
            Eccentricity of the system at each time step,
            with shape (num_snapshots, num_particles - 1)

        Notes
        -----
        - The function assumes that the first object is the central object.
        - C library function is not used here since this can be done with
            purely numpy vectorized operations. Nevertheless, we may consider
            implementing this in C library in the future.
        """
        num_snapshots = sol_state.shape[0]
        m_0 = sol_state[0, 0, 0]
        m = sol_state[0, 1:, 0]
        
        eccentricity = np.zeros(num_snapshots)

        x = sol_state[:, 1:, 1:4].copy() - sol_state[:, 0, 1:4].reshape(-1, 1, 3)
        v = sol_state[:, 1:, 4:7].copy() - sol_state[:, 0, 4:7].reshape(-1, 1, 3)

        denom = G * (m_0 + m)[np.newaxis, :, np.newaxis]
        eccentricity = (
            np.cross(v, np.cross(x, v)) / denom
            - x / np.linalg.norm(x, axis=2)[:, :, np.newaxis]
        )
        eccentricity = np.linalg.norm(eccentricity, axis=2)

        return eccentricity

    @staticmethod
    def compute_inclination(sol_state: np.ndarray) -> np.ndarray:
        """Compute the inclination using the sol_state array,
        assuming that the first object is the central object

        Parameters
        ----------
        sol_state : np.ndarray
            Solution state of the system

        Returns
        -------
        np.ndarray
            Inclination of the system at each time step, 
            with shape (num_snapshots, num_particles - 1)

        Notes
        -----
        - The function assumes that the first object is the central object.
        - C library function is not used here since this can be done with
          purely numpy vectorized operations. Nevertheless, we may consider
          implementing this in C library in the future.
        """
        num_snapshots = sol_state.shape[0]

        inclination = np.zeros(num_snapshots)

        x = sol_state[:, 1:, 1:4].copy() - sol_state[:, 0, 1:4].reshape(-1, 1, 3)
        v = sol_state[:, 1:, 4:7].copy() - sol_state[:, 0, 4:7].reshape(-1, 1, 3)

        unit_angular_momentum_vector = (
            np.cross(x, v) / np.linalg.norm(np.cross(x, v), axis=2)[:, :, np.newaxis]
        )
        unit_z = np.array([0, 0, 1])

        inclination = np.arccos(np.sum(unit_angular_momentum_vector * unit_z, axis=2))

        return inclination

    @staticmethod
    def plot_rel_energy_error(
        sol_energy: np.ndarray,
        sol_time: np.ndarray,
        is_log_y: bool = True,
        title: Optional[str] = None,
        xlabel: Optional[str] = "Time",
        ylabel: Optional[str] = "$(E_0 - E(t)) / E_0$",
        save_fig: bool = False,
        save_fig_path: Optional[str | Path] = None,
    ) -> None:
        if sol_energy[0] == 0.0:
            warnings.warn("The initial energy is zero.")
        rel_energy_error = np.abs((sol_energy - sol_energy[0]) / sol_energy[0])
        plotting.plot_quantity_against_time(
            quantity=rel_energy_error,
            sol_time=sol_time,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            is_log_y=is_log_y,
            save_fig=save_fig,
            save_fig_path=save_fig_path,
        )

    @staticmethod
    def plot_rel_linear_momentum_error(
        sol_linear_momentum: np.ndarray,
        sol_time: np.ndarray,
        is_log_y: bool = True,
        title: Optional[str] = None,
        xlabel: Optional[str] = "Time",
        ylabel: Optional[str] = "Relative linear momentum error",
        save_fig: bool = False,
        save_fig_path: Optional[str | Path] = None,
    ) -> None:
        if sol_linear_momentum[0] == 0.0:
            warnings.warn("The initial linear momentum is zero.")
        rel_linear_momentum_error = np.abs(
            (sol_linear_momentum - sol_linear_momentum[0]) / sol_linear_momentum[0]
        )
        plotting.plot_quantity_against_time(
            quantity=rel_linear_momentum_error,
            sol_time=sol_time,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            is_log_y=is_log_y,
            save_fig=save_fig,
            save_fig_path=save_fig_path,
        )

    @staticmethod
    def plot_rel_angular_momentum_error(
        sol_angular_momentum: np.ndarray,
        sol_time: np.ndarray,
        is_log_y: bool = True,
        title: Optional[str] = None,
        xlabel: Optional[str] = "Time",
        ylabel: Optional[str] = "$(L_0 - L(t)) / L_0$",
        save_fig: bool = False,
        save_fig_path: Optional[str | Path] = None,
    ) -> None:
        if sol_angular_momentum[0] == 0.0:
            warnings.warn("The initial angular momentum is zero.")
        angular_momentum_error = np.abs(
            (sol_angular_momentum - sol_angular_momentum[0]) / sol_angular_momentum[0]
        )
        plotting.plot_quantity_against_time(
            quantity=angular_momentum_error,
            sol_time=sol_time,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            is_log_y=is_log_y,
            save_fig=save_fig,
            save_fig_path=save_fig_path,
        )