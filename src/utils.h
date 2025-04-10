/**
 * \file utils.h
 * \brief Utility functions
 * 
 * \author Ching-Yin Ng
 */

#ifndef UTILS_H
#define UTILS_H

#include "common.h"
#include "pcg_basic.h"
#include "system.h"


/**
 * \brief Get current time as a decimal number of seconds using clock_gettime(CLOCK_MONOTONIC, )
 * 
 * \return Current time as a decimal number of seconds
 */
double grav_get_current_time(void);

/**
 * \brief Compute the energy of the system at a time step
 * 
 * \param[in] system Pointer to the gravitational system
 * 
 * \return Energy of the system
 */
double compute_energy(const System *__restrict system);

/**
 * \brief Initialize the PCG random number generator
 * 
 * \return Initialized PCG random number generator
 */
pcg32_random_t init_pcg_rng(void);

/**
 * \brief Generate a random number in the range [min, max)
 * 
 * \param min Minimum value of the range
 * \param max Maximum value of the range
 * \param rng Pointer to the PCG random number generator
 */
double grav_randrange(
    const double min,
    const double max,
    pcg32_random_t *rng
);

/**
 * \brief Convert Keplerian elements to Cartesian coordinates
 * 
 * \param[out] x Output position vector
 * \param[out] v Output velocity vector
 * \param[in] semi_major_axis Semi-major axis
 * \param[in] eccentricity Eccentricity
 * \param[in] inclination Inclination
 * \param[in] argument_of_periapsis Argument of periapsis
 * \param[in] longitude_of_ascending_node Longitude of ascending node
 * \param[in] true_anomaly True anomaly
 * \param[in] total_mass Total mass of the system
 * \param[in] G Gravitational constant
 */
void keplerian_to_cartesian(
    double *__restrict x,
    double *__restrict v,
    const double semi_major_axis,
    const double eccentricity,
    const double inclination,
    const double argument_of_periapsis,
    const double longitude_of_ascending_node,
    const double true_anomaly,
    const double total_mass,
    const double G
);

// /**
//  * \brief Compute the energy from solution state
//  * 
//  * \param num_particles Number of particles in the system
//  * \param m Pointer to the mass array
//  * \param G Gravitational constant
//  * \param npts Number of time steps
//  * \param count Pointer to the count variable
//  * \param energy Pointer to the energy array to be updated
//  * \param sol_state Pointer to the solution state array
//  * \param is_exit Pointer to the exit flag
//  */
// void compute_energy_python(
//     const int num_particles,
//     const double *__restrict m,
//     const double G,
//     const int npts,
//     int *__restrict count,
//     double *__restrict energy,
//     const double (*__restrict sol_state)[num_particles * 6],
//     int *__restrict is_exit
// );

// /**
//  * \brief Compute the linear momentum from solution state
//  * 
//  * \param num_particles Number of particles in the system
//  * \param m Pointer to the mass array
//  * \param npts Number of time steps
//  * \param count Pointer to the count variable
//  * \param linear_momentum Pointer to the linear momentum array to be updated
//  * \param sol_state Pointer to the solution state array
//  * \param is_exit Pointer to the exit flag
//  */
// void compute_linear_momentum_python(
//     const int num_particles,
//     const double *__restrict m,
//     const int npts,
//     int *__restrict count,
//     double *__restrict linear_momentum,
//     const double (*__restrict sol_state)[num_particles * 6],
//     int *__restrict is_exit
// );

// /**
//  * \brief Compute the angular momentum from solution state
//  * 
//  * \param num_particles Number of particles in the system
//  * \param m Pointer to the mass array
//  * \param npts Number of time steps
//  * \param count Pointer to the count variable
//  * \param angular_momentum Pointer to the angular momentum array to be updated
//  * \param sol_state Pointer to the solution state array
//  * \param is_exit Pointer to the exit flag
//  */
// void compute_angular_momentum_python(
//     const int num_particles,
//     const double *__restrict m,
//     const int npts,
//     int *__restrict count,
//     double *__restrict angular_momentum,
//     const double (*__restrict sol_state)[num_particles * 6],
//     int *__restrict is_exit
// );

#endif
