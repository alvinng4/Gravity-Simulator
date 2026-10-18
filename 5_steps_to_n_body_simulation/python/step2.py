import math
import timeit

import numpy as np

import common

INITIAL_CONDITION = "solar_system"
NUM_REPEATS = 10000


def main() -> None:
    # Get initial conditions
    system, _, _, _ = common.get_initial_conditions(INITIAL_CONDITION)

    ### Benchmark ###
    print("Benchmarking with 10000 repetitions")
    print()

    # Allocate memory
    a = np.zeros((system.num_particles, 3))

    # Acceleration 1
    run_time_1 = np.zeros(NUM_REPEATS)
    for i in range(NUM_REPEATS):
        start = timeit.default_timer()
        acceleration_1(a, system)
        end = timeit.default_timer()
        run_time_1[i] = end - start
    print(
        f"acceleration_1: {run_time_1.mean():.6f} +- {run_time_1.std(ddof=1) / math.sqrt(NUM_REPEATS):.3g} seconds"
    )

    # Acceleration 2
    run_time_2 = np.zeros(NUM_REPEATS)
    for i in range(NUM_REPEATS):
        start = timeit.default_timer()
        acceleration_2(a, system)
        end = timeit.default_timer()
        run_time_2[i] = end - start
    print(
        f"acceleration_2: {run_time_2.mean():.6f} +- {run_time_2.std(ddof=1) / math.sqrt(NUM_REPEATS):.3g} seconds"
    )

    # Acceleration 3
    run_time_3 = np.zeros(NUM_REPEATS)
    for i in range(NUM_REPEATS):
        start = timeit.default_timer()
        acceleration_3(a, system)
        end = timeit.default_timer()
        run_time_3[i] = end - start
    print(
        f"acceleration_3: {run_time_3.mean():.6f} +- {run_time_3.std(ddof=1) / math.sqrt(NUM_REPEATS):.3g} seconds"
    )

    # Acceleration 4
    run_time_4 = np.zeros(NUM_REPEATS)
    for i in range(NUM_REPEATS):
        start = timeit.default_timer()
        acceleration_4(a, system)
        end = timeit.default_timer()
        run_time_4[i] = end - start
    print(
        f"acceleration_4: {run_time_4.mean():.6f} +- {run_time_4.std(ddof=1) / math.sqrt(NUM_REPEATS):.3g} seconds"
    )

    ### Error check ###
    acceleration_1(a, system)
    a_1 = a.copy()
    acceleration_2(a, system)
    a_2 = a.copy()
    acceleration_3(a, system)
    a_3 = a.copy()
    acceleration_4(a, system)
    a_4 = a.copy()

    rel_error_2 = np.sum(np.abs(a_1 - a_2)) / np.sum(np.abs(a_1))
    rel_error_3 = np.sum(np.abs(a_1 - a_3)) / np.sum(np.abs(a_1))
    rel_error_4 = np.sum(np.abs(a_1 - a_4)) / np.sum(np.abs(a_1))

    print()
    print("Error check: (relative difference from acceleration_1)")
    print(f"acceleration_2: {rel_error_2:.3g}")
    print(f"acceleration_3: {rel_error_3:.3g}")
    print(f"acceleration_4: {rel_error_4:.3g}")


def acceleration_1(
    a: np.ndarray,
    system: common.System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Empty acceleration array
    a.fill(0.0)

    # Declare variables
    num_particles = system.num_particles
    x = system.x
    m = system.m
    G = system.G

    # Calculations
    for i in range(num_particles):
        for j in range(num_particles):
            if i == j:
                continue

            R = x[j] - x[i]
            a[i] += G * m[j] * R / (np.linalg.norm(R) ** 3)


def acceleration_2(
    a: np.ndarray,
    system: common.System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Empty acceleration array
    a.fill(0.0)

    # Declare variables
    num_particles = system.num_particles
    x = system.x
    m = system.m
    G = system.G

    # Calculations
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            R = x[j] - x[i]
            temp_value = G * R / (np.linalg.norm(R) ** 3)
            a[i] += temp_value * m[j]
            a[j] -= temp_value * m[i]


def acceleration_3(
    a: np.ndarray,
    system: common.System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Declare variables
    x = system.x
    m = system.m
    G = system.G

    # Compute the displacement vector
    r_ij = x[np.newaxis, :, :] - x[:, np.newaxis, :]

    # Compute the distance
    r_norm = np.linalg.norm(r_ij, axis=2)

    # Compute 1 / r^3
    inv_r_cubed = np.divide(
        1.0, r_norm * r_norm * r_norm, out=np.zeros_like(r_norm), where=r_norm > 0.0
    )

    # Compute the acceleration
    a[:] = G * np.sum(
        r_ij * inv_r_cubed[:, :, np.newaxis] * m[np.newaxis, :, np.newaxis], axis=1
    )


def acceleration_4(
    a: np.ndarray,
    system: common.System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Declare variables
    x = system.x
    m = system.m
    G = system.G

    # Compute the displacement vector
    r_ij = x[np.newaxis, :, :] - x[:, np.newaxis, :]

    # Compute the distance
    r_norm = np.linalg.norm(r_ij, axis=2)

    # Compute 1 / r^3
    inv_r_cubed = np.divide(
        1.0, r_norm * r_norm * r_norm, out=np.zeros_like(r_norm), where=r_norm > 0.0
    )

    # Compute the acceleration
    a[:] = G * np.einsum("ijk,ij,j->ik", r_ij, inv_r_cubed, m)


if __name__ == "__main__":
    main()
