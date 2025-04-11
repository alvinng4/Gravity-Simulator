"""
Utility functions for grav_sim API
"""

import ctypes
from pathlib import Path
from typing import Optional

import numpy as np


def load_c_lib(c_lib_path: Optional[str | Path] = None) -> ctypes.CDLL:
    """Load the C dynamic-link library

    Returns
    -------
    c_lib : ctypes.CDLL
        C dynamic-link library object

    Raises
    ------
    OSError
        If the platform is not supported
    FileNotFoundError
        If the C library is not found at the path
    """
    if not c_lib_path:
        search_path = Path(__file__).parent.parent
        c_lib_files = [str(p) for p in search_path.rglob("*libgrav_sim*")]
        if len(c_lib_files) == 0:
            raise FileNotFoundError(f"C library not found from path: {search_path}")

        c_lib_path = c_lib_files[0]

    return ctypes.cdll.LoadLibrary(str(c_lib_path))


def initialize_c_lib(c_lib: ctypes.CDLL) -> None:
    """Initialize C library

    Parameters
    ----------
    c_lib : ctypes.CDLL
        C dynamic-link library object

    Raises
    ------
    ValueError
        If any C functions that should be available are not found in the C library
    """
    c_lib.launch_simulation_python.restype = int
    c_lib.compute_energy_python.restype = None
    c_lib.compute_linear_momentum_python.restype = None
    c_lib.compute_angular_momentum_python.restype = None
    c_lib.load_built_in_system_python.restype = int
    c_lib.free_memory_int32.restype = None
    c_lib.free_memory_double.restype = None
    c_lib.keplerian_to_cartesian_python.restype = None


def keplerian_to_cartesian(
    c_lib: ctypes.CDLL,
    semi_major_axis: float,
    eccentricity: float,
    inclination: float,
    argument_of_periapsis: float,
    longitude_of_ascending_node: float,
    true_anomaly: float,
    total_mass: float,
    G: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert Keplerian elements to Cartesian coordinates

    Parameters
    ----------
    semi_major_axis : float
        Semi-major axis
    eccentricity : float
        Eccentricity
    inclination : float
        Inclination
    raan : float
        Right ascension of ascending node
    arg_periapsis : float
        Argument of periapsis
    true_anomaly : float
        True anomaly
    total_mass : float
        Total mass of the two bodies
    G : float
        Gravitational constant

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Position and velocity vectors in Cartesian coordinates
    """

    x = ctypes.c_double()
    y = ctypes.c_double()
    z = ctypes.c_double()
    v_x = ctypes.c_double()
    v_y = ctypes.c_double()
    v_z = ctypes.c_double()
    c_lib.keplerian_to_cartesian_python(
        ctypes.byref(x),
        ctypes.byref(y),
        ctypes.byref(z),
        ctypes.byref(v_x),
        ctypes.byref(v_y),
        ctypes.byref(v_z),
        ctypes.c_double(semi_major_axis),
        ctypes.c_double(eccentricity),
        ctypes.c_double(inclination),
        ctypes.c_double(argument_of_periapsis),
        ctypes.c_double(longitude_of_ascending_node),
        ctypes.c_double(true_anomaly),
        ctypes.c_double(total_mass),
        ctypes.c_double(G),
    )

    return np.array([x.value, y.value, z.value]), np.array(
        [v_x.value, v_y.value, v_z.value]
    )
