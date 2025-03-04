/**
 * \file acceleration.h
 * \author Ching Yin Ng
 * \brief Prototype of acceleration and acceleration-related
 *        functions for gravity-simulator
 */

#ifndef ACCELERATION_H
#define ACCELERATION_H

#include "gravity_sim.h"

#define ACCELERATION_METHOD_PAIRWISE 0
#define ACCELERATION_METHOD_MASSLESS 1
#define ACCELERATION_METHOD_BARNES_HUT 2

#ifdef USE_CUDA
    #define ACCELERATION_METHOD_CUDA_PAIRWISE 100
    #define ACCELERATION_METHOD_CUDA_PAIRWISE_FLOAT 101
    #define ACCELERATION_METHOD_CUDA_BARNES_HUT 102
    #define ACCELERATION_METHOD_CUDA_BARNES_HUT_FLOAT 103
#endif

/**
 * \brief Return acceleration method flag based on the input string
 * 
 * \param acceleration_method Name of the acceleration method
 * 
 * \retval SUCCESS If the acceleration method is recognized
 * \retval ERROR_UNKNOWN_ACCELERATION_METHOD If the acceleration method is not recognized
 */
int get_acceleration_method_flag(
    const char *__restrict acceleration_method,
    uint *__restrict acceleration_method_flag
);

/**
 * \brief Wrapper function for computing acceleration
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \retval SUCCESS If the computation is successful
 * \retval ERROR_UNKNOWN_ACCELERATION_CODE If the acceleration code is not recognized
 */
int acceleration(
    real *__restrict a,
    const System *__restrict system,
    AccelerationParam *__restrict acceleration_param
);

/**
 * \brief Compute acceleration with Barnes-Hut algorithm
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \retval SUCCESS If the computation is successful
 * \retval error code if errors occurred
 */
int acceleration_barnes_hut(
    real *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
);

#endif
