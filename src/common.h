/**
 * \file common.h
 * \brief Common definitions for the gravity simulation
 * 
 * \author Ching-Yin Ng
 */

#ifndef COMMON_H
#define COMMON_H

#include <stdbool.h>
#include <stdint.h>

/* For exporting functions in Windows for dynamic-link library */
#ifdef WIN32DLL_EXPORTS
    #define WIN32DLL_API __declspec(dllexport)
#else
    #define WIN32DLL_API 
#endif

/* Functions that are only used in the same file */
#define IN_FILE static

typedef int32_t int32;
typedef int64_t int64;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef struct IntegratorParam
{
    int integrator;
    double dt;
    double tolerance;
    double initial_dt;
    bool whfast_remove_invalid_particles;
} IntegratorParam;

typedef struct OutputParam
{
    int method;
    char *output_dir;
    bool output_initial;
    double output_interval;
    int coordinate_output_dtype;
    int velocity_output_dtype;
    int mass_output_dtype;

    int output_count_;
} OutputParam;

typedef struct SimulationStatus
{
    double t;
    double dt;
    int64 num_steps;
} SimulationStatus;



#endif
