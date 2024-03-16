import ctypes

import numpy as np
import rich.progress

from common import acceleration

class FIXED_STEP_SIZE_INTEGRATOR:
    def __init__(self, simulator, integrator, objects_count, x, v, m, G, dt, tf, is_c_lib):
        if is_c_lib == True: 
            self.c_lib = simulator.c_lib
            self.simulation_c_lib(integrator, objects_count, x, v, m, G, dt, tf)
        elif is_c_lib == False:
            self.simulation_numpy(integrator, objects_count, x, v, m, G, dt, tf)

    def simulation_c_lib(self, integrator, objects_count, x, v, m, G, dt, tf):
        t = ctypes.c_double(0.0)
        npts = int(np.floor((tf / dt))) + 1
        self.sol_state = np.zeros((npts, objects_count * 3 * 2))
        self.sol_state[0] = np.concatenate(
            (
                np.reshape(x, objects_count * 3),
                np.reshape(v, objects_count * 3),
            )
        )
        self.sol_time = np.linspace(
            0, dt * (npts - 1), npts
        )
        self.sol_dt = np.full(
            shape=(npts), fill_value=f"{dt}", dtype=float
        )
        progress_bar = rich.progress.Progress(
            rich.progress.BarColumn(),
            rich.progress.TextColumn("[green]{task.percentage:>3.0f}%"),
            rich.progress.TextColumn("•"),
            rich.progress.TimeElapsedColumn(),
            rich.progress.TextColumn("•"),
            rich.progress.TimeRemainingColumn(),
        )

        with progress_bar:
            match integrator:
                case "euler":
                    task = progress_bar.add_task("", total=tf)
                    while t.value <= self.sol_time[-1]:
                        self.c_lib.euler(
                            ctypes.c_int(objects_count),
                            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.byref(t),
                            ctypes.c_double(dt),
                            ctypes.c_double(tf),
                            ctypes.c_int(npts),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                            self.sol_state.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                        )
                        progress_bar.update(task, completed=t.value)
                    progress_bar.stop()
                case "euler_cromer":
                    task = progress_bar.add_task("", total=tf)
                    while t.value <= self.sol_time[-1]:
                        self.c_lib.euler_cromer(
                            ctypes.c_int(objects_count),
                            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.byref(t),
                            ctypes.c_double(dt),
                            ctypes.c_double(tf),
                            ctypes.c_int(npts),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                            self.sol_state.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                        )
                        progress_bar.update(task, completed=t.value)
                    progress_bar.stop()
                case "rk4":
                    task = progress_bar.add_task("", total=tf)
                    while t.value <= self.sol_time[-1]:
                        self.c_lib.rk4(
                            ctypes.c_int(objects_count),
                            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.byref(t),
                            ctypes.c_double(dt),
                            ctypes.c_double(tf),
                            ctypes.c_int(npts),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                            self.sol_state.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                        )
                        progress_bar.update(task, completed=t.value)
                    progress_bar.stop()
                case "leapfrog":
                    task = progress_bar.add_task("", total=tf)
                    while t.value <= self.sol_time[-1]:
                        self.c_lib.leapfrog(
                            ctypes.c_int(objects_count),
                            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            v.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.byref(t),
                            ctypes.c_double(dt),
                            ctypes.c_double(tf),
                            ctypes.c_int(npts),
                            m.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            ctypes.c_double(G),
                            self.sol_state.ctypes.data_as(
                                ctypes.POINTER(ctypes.c_double)
                            ),
                        )
                        progress_bar.update(task, completed=t.value)
                    progress_bar.stop() 

    def simulation_numpy(self, integrator, objects_count, x, v, m, G, dt, tf):
        npts = int(np.floor((tf / dt))) + 1
        self.sol_state = np.zeros((npts, objects_count * 3 * 2))
        self.sol_state[0] = np.concatenate(
            (
                np.reshape(x, objects_count * 3),
                np.reshape(v, objects_count * 3),
            )
        )
        self.sol_time = np.linspace(
            0, dt * (npts - 1), npts
        )
        self.sol_dt = np.full(
            shape=(npts), fill_value=f"{dt}", dtype=float
        )
        progress_bar = rich.progress.Progress(
            rich.progress.BarColumn(),
            rich.progress.TextColumn("[green]{task.percentage:>3.0f}%"),
            rich.progress.TextColumn("•"),
            rich.progress.TimeElapsedColumn(),
            rich.progress.TextColumn("•"),
            rich.progress.TimeRemainingColumn(),
        )

        with progress_bar:
            match integrator:
                case "euler":
                    for count in progress_bar.track(range(npts - 1)):
                        a = acceleration(
                            objects_count, x, m, G
                        )
                        x, v = self.euler(x, v, a, dt)
                        self.sol_state[count + 1] = np.concatenate(
                            (
                                np.reshape(x, objects_count * 3),
                                np.reshape(v, objects_count * 3),
                            )
                        )
                case "euler_cromer":
                    for count in progress_bar.track(range(npts - 1)):
                        a = acceleration(
                            objects_count, x, m, G
                        )
                        x, v = self.euler_cromer(x, v, a, dt)
                        self.sol_state[count + 1] = np.concatenate(
                            (
                                np.reshape(x, objects_count * 3),
                                np.reshape(v, objects_count * 3),
                            )
                        )

                case "rk4":
                    for count in progress_bar.track(range(npts - 1)):
                        x, v = self.rk4(
                            objects_count,
                            x,
                            v,
                            m,
                            G,
                            dt,
                        )
                        self.sol_state[count + 1] = np.concatenate(
                            (
                                np.reshape(x, objects_count * 3),
                                np.reshape(v, objects_count * 3),
                            )
                        )

                case "leapfrog":
                    a = acceleration(objects_count, x, m, G)
                    for count in progress_bar.track(range(npts - 1)):
                        x, v, a = self.leapfrog(
                            objects_count,
                            x,
                            v,
                            a,
                            m,
                            dt,
                            G,
                        )
                        self.sol_state[count + 1] = np.concatenate(
                            (
                                np.reshape(x, objects_count * 3),
                                np.reshape(v, objects_count * 3),
                            )
                        )    
                             
    @staticmethod        
    def euler(x, v, a, dt):
        return x + v * dt, v + a * dt

    @staticmethod    
    def euler_cromer(x, v, a, dt):
        v = v + a * dt
        x = x + v * dt
        return x, v

    @staticmethod    
    def rk4(objects_count, x, v, m, G, dt):
        vk1 = acceleration(objects_count, x, m, G)
        xk1 = v

        vk2 = acceleration(objects_count, x + 0.5 * xk1 * dt, m, G)
        xk2 = v + 0.5 * vk1 * dt

        vk3 = acceleration(objects_count, x + 0.5 * xk2 * dt, m, G)
        xk3 = v + 0.5 * vk2 * dt

        vk4 = acceleration(objects_count, x + xk3 * dt, m, G)
        xk4 = v + vk3 * dt

        v = v + dt * (vk1 + 2 * vk2 + 2 * vk3 + vk4) / 6.0
        x = x + dt * (xk1 + 2 * xk2 + 2 * xk3 + xk4) / 6.0

        return x, v

    @staticmethod    
    def leapfrog(objects_count, x, v, a, m, dt, G):
        a_0 = a
        x = x + v * dt + a_0 * 0.5 * dt * dt
        a_1 = acceleration(objects_count, x, m, G)
        v = v + (a_0 + a_1) * 0.5 * dt

        return x, v, a_1


