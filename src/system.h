/**
 * \file system.h
 * \brief Header file for the system module
 */

#ifndef SYSTEM_H
#define SYSTEM_H

#include "error.h"
#include "settings.h"

typedef struct System
{
    int num_particles;
    int *particle_ids;
    double *x;
    double *v;
    double *m;
    double G;
} System;

typedef struct CosmologicalSystem
{
    int num_particles;
    int *particle_ids;
    double *x;
    double *v;
    double *m;
    double G;
    double h0;
    double omega_m;
    double omega_lambda;
    double omega_k;
    double box_center[3];
    double box_width;
    double unit_mass;
    double unit_length;
    double unit_time;
} CosmologicalSystem;

/**
 * \brief Get a new system structure with uninitialized memory
 * 
 * \return System
 */
System get_new_system(void);

/**
 * \brief Get a new system with initialized memory for the given number of particles
 * 
 * \param[out] system Pointer to the system to be initialized
 * \param[in] num_particles Number of particles
 * 
 * \return ErrorStatus
 */
ErrorStatus get_initialized_system(
    System *__restrict system,
    const int num_particles
);

/**
 * \brief Finalize the system by checking the system members
 * 
 * \param[in, out] system Pointer to the system to be finalized
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_POINTER_ERROR if system or its members are NULL
 * \exception GRAV_VALUE_ERROR if the number of particles is less than 1
 * \exception GRAV_VALUE_ERROR if the gravitational constant is not positive
 */
ErrorStatus finalize_system(System *__restrict system);

/**
 * \brief Get a new cosmological system structure with uninitialized memory
 * 
 * \return CosmologicalSystem
 */
CosmologicalSystem get_new_cosmological_system(void);

/**
 * \brief Get a new cosmological system with initialized memory for the given number of particles
 * 
 * \param[out] system Pointer to the cosmological system to be initialized
 * \param[in] num_particles Number of particles
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if failed to allocate memory
 * \exception GRAV_POINTER_ERROR if system or its members are NULL
 */
ErrorStatus get_initialized_cosmological_system(
    CosmologicalSystem *__restrict system,
    const int num_particles
);

/**
 * \brief Finalize the cosmological system by checking the system members
 * 
 * \param[in, out] system Pointer to the cosmological system to be finalized
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_POINTER_ERROR if system or its members are NULL
 * \exception GRAV_VALUE_ERROR if the number of particles is less than 1
 * \exception GRAV_VALUE_ERROR if the gravitational constant is not positive
 * \exception GRAV_VALUE_ERROR if the Hubble constant is not positive
 * \exception GRAV_VALUE_ERROR if omega_m is not positive
 * \exception GRAV_VALUE_ERROR if omega_lambda is not positive
 * \exception GRAV_VALUE_ERROR if box_width is not positive
 */
ErrorStatus finalize_cosmological_system(CosmologicalSystem *__restrict system);

/**
 * \brief Free the memory allocated for the cosmological system
 * 
 * \param[in] system Pointer to the cosmological system to be freed
 */
void free_cosmological_system(CosmologicalSystem *__restrict system);

/**
 * \brief Free the memory allocated for the system
 * 
 * \param[in] system Pointer to the system to be freed
 */
void free_system(System *__restrict system);

/**
 * \brief Set the boundary condition for the system
 * 
 * \param[in, out] system Pointer to the system
 * \param[in] settings Pointer to the settings
 */
ErrorStatus set_boundary_condition(
    System *__restrict system,
    const Settings *__restrict settings
);

/**
 * \brief Check for invalid indices in a double array
 * 
 * \param[out] has_invalid_idx Pointer to a boolean variable indicating if there are invalid indices
 * \param[out] invalid_idx_array Pointer to an array of invalid indices
 * \param[in] array Pointer to the double array to be checked
 * \param[in] arr_size Size of the array
 */
ErrorStatus check_invalid_idx_double(
    bool *__restrict has_invalid_idx,
    int **invalid_idx_array,
    const double *__restrict array,
    const int arr_size
);

/**
 * \brief Check and remove invalid particles from the system
 * 
 * \param[in, out] system Pointer to the system
 * \param[in] settings Pointer to the settings
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if failed to allocate memory
 * \exception GRAV_POINTER_ERROR if system or its members are NULL
 * \exception Other exceptions if failed to remove particles
 */
ErrorStatus check_and_remove_invalid_particles(
    System *__restrict system,
    const Settings *__restrict settings
);

/**
 * \brief Remove invalid particles from the system
 * 
 * \param[in, out] system Pointer to the system
 * \param[in] remove_idx_list List of indices to be removed
 * \param[in] num_to_remove Number of particles to be removed
 * \param[in] settings Pointer to the settings
 * 
 * \return ErrorStatus
 */
ErrorStatus remove_invalid_particles(
    System *__restrict system,
    const int *__restrict remove_idx_list,
    const int num_to_remove,
    const Settings *__restrict settings
);

/**
 * \brief Remove a list of particles from the system
 * 
 * \param[in, out] system Pointer to the system
 * \param[in] remove_idx_list List of indices to be removed
 * \param[in] num_to_remove Number of particles to be removed
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if failed to reallocate memory for the system
 */
ErrorStatus remove_particles(
    System *__restrict system,
    const int *__restrict remove_idx_list,
    const int num_to_remove
);

/**
 * \brief Remove a list of particles from a double array
 * 
 * \param[out] arr The double array to be modified
 * \param[in] remove_idx_list List of indices to be removed
 * \param[in] num_to_remove Number of particles to be removed
 * \param[in] dim Dimension of the array
 * \param[in] original_size Original size of the array
 * 
 * \return ErrorStatus
 */
ErrorStatus remove_particle_from_double_arr(
    double *__restrict arr,
    const int *__restrict remove_idx_list,
    const int num_to_remove,
    const int dim,
    const int original_size
);

/**
 * \brief Initialize a built-in system with the given name
 * 
 * \param[out] system Pointer to the system to be initialized
 * \param[in] system_name Name of the built-in system
 * \param[in] is_memory_initialized Flag indicating if the memory is already initialized
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if failed to allocate memory
 * \exception GRAV_POINTER_ERROR if system or system_name is NULL
 * \exception GRAV_VALUE_ERROR if system_name is not recognized
 * \exception GRAV_VALUE_ERROR if the initialized memory is less than the required size
 */
ErrorStatus initialize_built_in_system(
    System *__restrict system,
    const char *__restrict system_name,
    const bool is_memory_initialized
);

/**
 * \brief Set the center of mass of the system to zero
 * 
 * \param[in, out] system Pointer to the system
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_POINTER_ERROR if system or its members are NULL
 * \exception GRAV_VALUE_ERROR if total mass is non-positive
 * \exception GRAV_VALUE_ERROR if the center of mass is invalid
 */
ErrorStatus system_set_center_of_mass_zero(System *__restrict system);

/**
 * \brief Set the total momentum of the system to zero
 * 
 * \param[in, out] system Pointer to the system
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_POINTER_ERROR if system or its members are NULL
 * \exception GRAV_VALUE_ERROR if total mass is non-positive
 * \exception GRAV_VALUE_ERROR if the V_CM is invalid
 */
ErrorStatus system_set_total_momentum_zero(System *__restrict system);

/**
 * \brief Sort the system by distance from a primary particle
 * 
 * \param[in, out] system Pointer to the system
 * \param[in] primary_particle_id ID of the primary particle
 * 
 * \return ErrorStatus
 */
ErrorStatus system_sort_by_distance(
    System *__restrict system,
    const int primary_particle_id
);

/**
 * \brief Set periodic boundary conditions for the cosmological system
 * 
 * \param[in, out] system Pointer to the cosmological system
 */
void set_periodic_boundary_conditions(CosmologicalSystem *__restrict system);

#endif
