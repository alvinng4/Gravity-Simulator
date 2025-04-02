/**
 * \file system.h
 * \brief Header file for the system module
 * 
 * \author Ching-Yin Ng
 * \date April 2025
 */

#ifndef SYSTEM_H
#define SYSTEM_H

#include "error.h"

typedef struct System
{
    int objects_count;
    int *particle_ids;
    double *x;
    double *v;
    double *m;
    double G;
} System;

/**
 * \brief Get a new system structure with uninitialized memory
 * 
 * \return System
 */
System get_new_system(void);

/**
 * \brief Get a new system with initialized memory for the given number of objects
 * 
 * \param[out] system Pointer to the system to be initialized
 * \param[in] objects_count Number of objects
 * 
 * \return ErrorStatus
 */
ErrorStatus get_initialized_system(
    System *__restrict system,
    const int objects_count
);

/**
 * \brief Free the memory allocated for the system
 * 
 * \param[in] system Pointer to the system to be freed
 */
void free_system(System *__restrict system);

/**
 * \brief Remove invalid particles from the system
 * 
 * \param[in, out] system Pointer to the system
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if failed to allocate memory
 * \exception GRAV_POINTER_ERROR if system or its members are NULL
 * \exception Other exceptions if failed to remove particles
 */
ErrorStatus remove_invalid_particles(System *__restrict system);

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
 * \brief Initialize a built-in system with the given name
 * 
 * \param[out] system Pointer to the system to be initialized
 * \param[in] system_name Name of the built-in system
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if failed to allocate memory
 * \exception GRAV_POINTER_ERROR if system or system_name is NULL
 * \exception GRAV_VALUE_ERROR if system_name is not recognized
 */
ErrorStatus initialize_built_in_system(
    System *__restrict system,
    const char *__restrict system_name
);

#endif
