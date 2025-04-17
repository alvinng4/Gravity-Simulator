"""
Simulator module for grav_sim
"""

import ctypes
import time
import threading
from queue import Queue

import numpy as np

from .system import System
from .parameters import AccelerationParam, IntegratorParam, OutputParam, Settings


class Simulator:
    DAYS_PER_YEAR = 365.242189

    def __init__(self, c_lib: ctypes.CDLL) -> None:
        self.c_lib = c_lib

    def launch_simulation(
        self,
        system: System,
        acceleration_param: AccelerationParam,
        integrator_param: IntegratorParam,
        output_param: OutputParam,
        settings: Settings,
        tf: float,
    ) -> None:
        """Launch simulation

        Parameters
        ----------
        system : System
            Gravitational system
        acceleration_param : AccelerationParam
            Acceleration parameters
        integrator_param : IntegratorParam
            Integrator parameters
        output_param : OutputParam
            Output parameters
        settings : Settings
            Simulation settings
        tf : float
            Simulation time

        Raises
        -------
        RuntimeError
            If simulation fails
        """
        num_particles = ctypes.c_int32(system.num_particles)
        new_particle_ids = ctypes.POINTER(ctypes.c_int32)()
        new_x = ctypes.POINTER(ctypes.c_double)()
        new_v = ctypes.POINTER(ctypes.c_double)()
        new_m = ctypes.POINTER(ctypes.c_double)()

        is_exit = ctypes.c_bool(False)

        # Create a thread to stop the simulation if user sends KeyboardInterrupt
        def simulation_wrapper(c_lib_launch_simulation_python, return_queue, *args):
            return_queue.put(c_lib_launch_simulation_python(*args))

        queue: Queue = Queue()
        simulation_thread = threading.Thread(
            target=simulation_wrapper,
            args=(
                self.c_lib.launch_simulation_python,
                queue,
                ctypes.byref(num_particles),
                system.particle_ids.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                system.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                system.v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                system.m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.byref(new_particle_ids),
                ctypes.byref(new_x),
                ctypes.byref(new_v),
                ctypes.byref(new_m),
                ctypes.c_double(system.G),
                ctypes.c_int32(integrator_param._integrator),
                ctypes.c_double(integrator_param.dt),
                ctypes.c_double(integrator_param.tolerance),
                ctypes.c_double(integrator_param.initial_dt),
                ctypes.c_bool(integrator_param.whfast_remove_invalid_particles),
                ctypes.c_int32(acceleration_param._method),
                ctypes.c_double(acceleration_param.opening_angle),
                ctypes.c_double(acceleration_param.softening_length),
                ctypes.c_int32(acceleration_param.max_num_particles_per_leaf),
                ctypes.c_int32(output_param._method),
                output_param.output_dir.encode("utf-8"),
                ctypes.c_bool(output_param.output_initial),
                ctypes.c_double(output_param.output_interval),
                ctypes.c_int32(output_param._coordinate_output_dtype),
                ctypes.c_int32(output_param._velocity_output_dtype),
                ctypes.c_int32(output_param._mass_output_dtype),
                ctypes.c_int32(settings._verbose),
                ctypes.c_bool(settings.enable_progress_bar),
                ctypes.byref(is_exit),
                ctypes.c_double(tf),
            ),
        )

        ### Launch simulation ###
        try:
            simulation_thread.start()

            while simulation_thread.is_alive():
                time.sleep(0.05)

        except KeyboardInterrupt:
            print("\nKeyboardInterrupt detected. Stopping simulation...\n")
            is_exit.value = True
            simulation_thread.join()
            raise KeyboardInterrupt

        # End simulation and get return value
        simulation_thread.join()
        ret_val = queue.get()
        if ret_val != 0:
            raise RuntimeError(
                f"Simulation failed. Please check the C library traceback."
            )

        # Update system with new data (just to be safe since 
        # we don't know if the memory will be changed in the C library)
        system.particle_ids = np.ctypeslib.as_array(
            new_particle_ids,
            shape=(num_particles.value,),
        )
        system.x = np.ctypeslib.as_array(
            new_x,
            shape=(num_particles.value, 3),
        )
        system.v = np.ctypeslib.as_array(
            new_v,
            shape=(num_particles.value, 3),
        )
        system.m = np.ctypeslib.as_array(
            new_m,
            shape=(num_particles.value,),
        )
        system.num_particles = num_particles.value
