from pathlib import Path


class AccelerationParam:
    """Acceleration parameters for gravity simulation

    Attributes
    ----------
    method : str
        Acceleration method
    softening_length : float
        Softening length for acceleration
    pm_grid_size : int
        Grid size for Particle Mesh method
    opening_angle : float
        Opening angle for Barnes-Hut method
    max_num_particles_per_leaf : int
        Maximum number of particles per leaf for Barnes-Hut method
    """

    AVAILABLE_ACCELERATION_METHODS = [
        "pairwise",
        "massless",
        "barnes_hut",
        "particle_mesh",
    ]
    ACCELERATION_METHODS_TO_ENCODING = {
        "pairwise": 1,
        "massless": 2,
        "barnes_hut": 3,
        "particle_mesh": 4,
    }
    ACCELERATION_ENCODING_TO_METHODS = {
        v: k for k, v in ACCELERATION_METHODS_TO_ENCODING.items()
    }

    def __init__(self) -> None:
        self.method = "pairwise"
        self.opening_angle = 1.0
        self.softening_length = 0.0
        self.pm_grid_size = 128
        self.max_num_particles_per_leaf = 1

    @property
    def method(self) -> str:
        """Get the acceleration method"""
        return self.ACCELERATION_ENCODING_TO_METHODS[self._method]

    @method.setter
    def method(self, value: str) -> None:
        """Set the acceleration method"""
        value = value.lower().strip()
        if value not in self.AVAILABLE_ACCELERATION_METHODS:
            raise ValueError(
                f"Invalid acceleration method: {value}. Available methods: {self.AVAILABLE_ACCELERATION_METHODS}"
            )
        self._method: int = self.ACCELERATION_METHODS_TO_ENCODING[value]

    @property
    def opening_angle(self) -> float:
        """Get the opening angle for Barnes-Hut method"""
        return self._opening_angle

    @opening_angle.setter
    def opening_angle(self, value: float) -> None:
        """Set the opening angle for Barnes-Hut method"""
        value = float(value)
        if value <= 0:
            raise ValueError(f"Invalid opening angle: {value}. Must be positive.")
        self._opening_angle = value

    @property
    def softening_length(self) -> float:
        """Get the softening length"""
        return self._softening_length

    @softening_length.setter
    def softening_length(self, value: float) -> None:
        """Set the softening length"""
        value = float(value)
        if value < 0:
            raise ValueError(
                f"Invalid softening length: {value}. Must be non-negative."
            )
        self._softening_length = value

    @property
    def pm_grid_size(self) -> int:
        """Get the grid size for Particle Mesh method"""
        return self._pm_grid_size

    @pm_grid_size.setter
    def pm_grid_size(self, value: int) -> None:
        """Set the grid size for Particle Mesh method"""
        value = int(value)
        if value <= 0:
            raise ValueError(f"Invalid grid size: {value}. Must be positive.")
        self._pm_grid_size = value

    @property
    def max_num_particles_per_leaf(self) -> int:
        """Get the maximum number of particles per leaf for Barnes-Hut method"""
        return self._max_num_particles_per_leaf

    @max_num_particles_per_leaf.setter
    def max_num_particles_per_leaf(self, value: int) -> None:
        """Set the maximum number of particles per leaf for Barnes-Hut method"""
        value = int(value)
        if value <= 0:
            raise ValueError(
                f"Invalid maximum number of particles per leaf: {value}. Must be positive."
            )
        self._max_num_particles_per_leaf = value


class IntegratorParam:
    """Integrator parameters for gravity simulation

    Attributes
    ----------
    integrator : str
        Integrator
    dt : float
        Time step
    tolerance : float
        Tolerance for adaptive step size integrators
    initial_dt : float
        Initial time step for adaptive step size integrators
    whfast_remove_invalid_particles : bool
        Flag to indicate whether to remove invalid particles in WHFast integrator
    """

    AVAILABLE_INTEGRATORS = [
        "euler",
        "euler_cromer",
        "rk4",
        "leapfrog",
        "rkf45",
        "dopri",
        "dverk",
        "rkf78",
        "ias15",
        "whfast",
    ]
    FIXED_STEP_SIZE_INTEGRATORS = ["euler", "euler_cromer", "rk4", "leapfrog", "whfast"]
    ADAPTIVE_STEP_SIZE_INTEGRATORS = ["rkf45", "dopri", "dverk", "rkf78", "ias15"]

    INTEGRATOR_TO_ENCODING = {
        "euler": 1,
        "euler_cromer": 2,
        "rk4": 3,
        "leapfrog": 4,
        "rkf45": 5,
        "dopri": 6,
        "dverk": 7,
        "rkf78": 8,
        "ias15": 9,
        "whfast": 10,
    }
    ENCODING_TO_INTEGRATORS = {k: v for v, k in INTEGRATOR_TO_ENCODING.items()}

    # Recommended settings for built-in systems with IAS15 integrator
    RECOMMENDED_IAS15_SETTINGS_BUILT_IN_SYSTEMS = {
        # "template": ["tf", "tf unit", "tolerance", "storing_freq"],
        "circular_binary_orbit": [50.0, "days", 1e-9, 1],
        "eccentric_binary_orbit": [2.6, "years", 1e-9, 1],
        "3d_helix": [20.0, "days", 1e-9, 1],
        "sun_earth_moon": [1.0, "years", 1e-9, 1],
        "figure-8": [20.0, "days", 1e-9, 1],
        "pyth-3-body": [70.0, "days", 1e-9, 1],
        "solar_system": [200.0, "years", 1e-9, 1],
        "solar_system_plus": [250.0, "years", 1e-9, 1],
    }

    def __init__(self) -> None:
        self.integrator = "euler"
        self.dt = -1.0
        self.tolerance = -1.0
        self.initial_dt = -1.0
        self.whfast_remove_invalid_particles = True

    @property
    def integrator(self) -> str:
        """Get the integrator"""
        return self.ENCODING_TO_INTEGRATORS[self._integrator]

    @integrator.setter
    def integrator(self, value: str) -> None:
        """Set the integrator"""
        value = value.lower().strip()
        if value not in self.AVAILABLE_INTEGRATORS:
            raise ValueError(
                f"Invalid integrator: {value}. Available integrators: {self.AVAILABLE_INTEGRATORS}"
            )
        self._integrator: int = self.INTEGRATOR_TO_ENCODING[value]

    @property
    def dt(self) -> float:
        """Get the time step size"""
        return self._dt

    @dt.setter
    def dt(self, value: float) -> None:
        """Set the time step size"""
        self._dt = float(value)

    @property
    def tolerance(self) -> float:
        """Get the tolerance for adaptive step size integrators"""
        return self._tolerance

    @tolerance.setter
    def tolerance(self, value: float) -> None:
        """Set the tolerance for adaptive step size integrators"""
        self._tolerance = float(value)

    @property
    def initial_dt(self) -> float:
        """Get the initial time step size for adaptive step size integrators"""
        return self._initial_dt

    @initial_dt.setter
    def initial_dt(self, value: float) -> None:
        """Set the initial time step size for adaptive step size integrators"""
        self._initial_dt = float(value)

    @property
    def whfast_remove_invalid_particles(self) -> bool:
        """Get the flag for removing invalid particles in WHFast integrator"""
        return self._whfast_remove_invalid_particles

    @whfast_remove_invalid_particles.setter
    def whfast_remove_invalid_particles(self, value: bool) -> None:
        """Set the flag for removing invalid particles in WHFast integrator"""
        self._whfast_remove_invalid_particles = bool(value)


class OutputParam:
    """Output parameters for gravity simulation

    Attributes
    ----------
    method : str
        Output method
    output_dir : str
        Output directory
    output_initial : bool
        Flag to indicate whether to output initial state
    output_interval : float
        Output interval
    coordinate_output_dtype : str
        Data type for coordinate output
    velocity_output_dtype : str
        Data type for velocity output
    mass_output_dtype : str
        Data type for mass output
    """

    AVAILABLE_OUTPUT_METHODS = ["disabled", "csv", "hdf5"]
    AVAILABLE_OUTPUT_DTYPE = ["float", "double"]

    OUTPUT_METHODS_TO_ENCODING = {
        "disabled": 1,
        "csv": 2,
        "hdf5": 3,
    }
    OUTPUT_ENCODING_TO_METHODS = {v: k for k, v in OUTPUT_METHODS_TO_ENCODING.items()}

    OUTPUT_DTYPE_TO_ENCODING = {
        "float": 1,
        "double": 2,
    }
    OUTPUT_ENCODING_TO_DTYPE = {v: k for k, v in OUTPUT_DTYPE_TO_ENCODING.items()}

    def __init__(self) -> None:
        self.method = "disabled"
        self.output_dir = "tmp/"
        self.output_initial = True
        self.output_interval = 1.0
        self.coordinate_output_dtype = "double"
        self.velocity_output_dtype = "double"
        self.mass_output_dtype = "double"

    @property
    def method(self) -> str:
        """Get the output method"""
        return self.OUTPUT_ENCODING_TO_METHODS[self._method]

    @method.setter
    def method(self, value: str) -> None:
        """Set the output method"""
        value = value.lower().strip()
        if value not in self.AVAILABLE_OUTPUT_METHODS:
            raise ValueError(
                f"Invalid output method: {value}. Available methods: {self.AVAILABLE_OUTPUT_METHODS}"
            )
        self._method: int = self.OUTPUT_METHODS_TO_ENCODING[value]

    @property
    def output_dir(self) -> str:
        """Get the output directory"""
        return self._output_dir

    @output_dir.setter
    def output_dir(self, value: str) -> None:
        """Set the output directory"""
        value = str(Path(value).resolve())
        if not value.endswith("/"):
            value += "/"
        self._output_dir = value

    @property
    def output_initial(self) -> bool:
        """Get the flag for outputting initial state"""
        return self._output_initial

    @output_initial.setter
    def output_initial(self, value: bool) -> None:
        """Set the flag for outputting initial state"""
        self._output_initial = bool(value)

    @property
    def output_interval(self) -> float:
        """Get the output interval"""
        return self._output_interval

    @output_interval.setter
    def output_interval(self, value: float) -> None:
        """Set the output interval"""
        self._output_interval = float(value)

    @property
    def coordinate_output_dtype(self) -> str:
        """Get the data type for coordinate output"""
        return self.OUTPUT_ENCODING_TO_DTYPE[self._coordinate_output_dtype]

    @coordinate_output_dtype.setter
    def coordinate_output_dtype(self, value: str) -> None:
        """Set the data type for coordinate output"""
        value = value.lower().strip()
        if value not in self.AVAILABLE_OUTPUT_DTYPE:
            raise ValueError(
                f"Invalid coordinate output data type: {value}. Available types: {self.AVAILABLE_OUTPUT_DTYPE}"
            )
        self._coordinate_output_dtype: int = self.OUTPUT_DTYPE_TO_ENCODING[value]

    @property
    def velocity_output_dtype(self) -> str:
        """Get the data type for velocity output"""
        return self.OUTPUT_ENCODING_TO_DTYPE[self._velocity_output_dtype]

    @velocity_output_dtype.setter
    def velocity_output_dtype(self, value: str) -> None:
        """Set the data type for velocity output"""
        value = value.lower().strip()
        if value not in self.AVAILABLE_OUTPUT_DTYPE:
            raise ValueError(
                f"Invalid velocity output data type: {value}. Available types: {self.AVAILABLE_OUTPUT_DTYPE}"
            )
        self._velocity_output_dtype: int = self.OUTPUT_DTYPE_TO_ENCODING[value]

    @property
    def mass_output_dtype(self) -> str:
        """Get the data type for mass output"""
        return self.OUTPUT_ENCODING_TO_DTYPE[self._mass_output_dtype]

    @mass_output_dtype.setter
    def mass_output_dtype(self, value: str) -> None:
        """Set the data type for mass output"""
        value = value.lower().strip()
        if value not in self.AVAILABLE_OUTPUT_DTYPE:
            raise ValueError(
                f"Invalid mass output data type: {value}. Available types: {self.AVAILABLE_OUTPUT_DTYPE}"
            )
        self._mass_output_dtype: int = self.OUTPUT_DTYPE_TO_ENCODING[value]


class Settings:
    """Settings for gravity simulation

    Attributes
    ----------
    verbose : str
        Verbosity level
    enable_progress_bar : bool
        Flag to enable progress bar
    """

    AVAILABLE_VERBOSITY_LEVELS = ["ignore_all", "ignore_info", "normal", "verbose"]
    VERBOSITY_TO_ENCODING = {
        "ignore_all": 0,
        "ignore_info": 1,
        "normal": 2,
        "verbose": 3,
    }
    ENCODING_TO_VERBOSITY = {v: k for k, v in VERBOSITY_TO_ENCODING.items()}

    def __init__(self) -> None:
        self.verbose = "normal"
        self.enable_progress_bar = True

    @property
    def verbose(self) -> str:
        """Get the verbosity level"""
        return self.ENCODING_TO_VERBOSITY[self._verbose]

    @verbose.setter
    def verbose(self, value: str) -> None:
        """Set the verbosity level"""
        value = value.lower().strip()
        if value not in self.AVAILABLE_VERBOSITY_LEVELS:
            raise ValueError(
                f"Invalid verbosity level: {value}. Available levels: {self.AVAILABLE_VERBOSITY_LEVELS}"
            )
        self._verbose: int = self.VERBOSITY_TO_ENCODING[value]

    @property
    def enable_progress_bar(self) -> bool:
        """Get the flag for enabling progress bar"""
        return self._enable_progress_bar

    @enable_progress_bar.setter
    def enable_progress_bar(self, value: bool) -> None:
        """Set the flag for enabling progress bar"""
        self._enable_progress_bar = bool(value)
