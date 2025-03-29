/**
 * \file utils.h
 * \brief Utility functions
 * 
 * \author Ching-Yin Ng
 * \date March 2025
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
 * \brief Free memory allocated for a double array
 * 
 * \param ptr Pointer to the double array
 */
void free_memory_double(double *__restrict ptr);

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

// /**
//  * \brief Compute the energy from solution state
//  * 
//  * \param objects_count Number of objects in the system
//  * \param m Pointer to the mass array
//  * \param G Gravitational constant
//  * \param npts Number of time steps
//  * \param count Pointer to the count variable
//  * \param energy Pointer to the energy array to be updated
//  * \param sol_state Pointer to the solution state array
//  * \param is_exit Pointer to the exit flag
//  */
// void compute_energy_python(
//     const int objects_count,
//     const double *__restrict m,
//     const double G,
//     const int npts,
//     int *__restrict count,
//     double *__restrict energy,
//     const double (*__restrict sol_state)[objects_count * 6],
//     int *__restrict is_exit
// );

// /**
//  * \brief Compute the linear momentum from solution state
//  * 
//  * \param objects_count Number of objects in the system
//  * \param m Pointer to the mass array
//  * \param npts Number of time steps
//  * \param count Pointer to the count variable
//  * \param linear_momentum Pointer to the linear momentum array to be updated
//  * \param sol_state Pointer to the solution state array
//  * \param is_exit Pointer to the exit flag
//  */
// void compute_linear_momentum_python(
//     const int objects_count,
//     const double *__restrict m,
//     const int npts,
//     int *__restrict count,
//     double *__restrict linear_momentum,
//     const double (*__restrict sol_state)[objects_count * 6],
//     int *__restrict is_exit
// );

// /**
//  * \brief Compute the angular momentum from solution state
//  * 
//  * \param objects_count Number of objects in the system
//  * \param m Pointer to the mass array
//  * \param npts Number of time steps
//  * \param count Pointer to the count variable
//  * \param angular_momentum Pointer to the angular momentum array to be updated
//  * \param sol_state Pointer to the solution state array
//  * \param is_exit Pointer to the exit flag
//  */
// void compute_angular_momentum_python(
//     const int objects_count,
//     const double *__restrict m,
//     const int npts,
//     int *__restrict count,
//     double *__restrict angular_momentum,
//     const double (*__restrict sol_state)[objects_count * 6],
//     int *__restrict is_exit
// );

#endif
