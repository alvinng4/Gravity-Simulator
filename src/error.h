/**
 * \file error.h
 * \author Ching Yin Ng
 * \brief Error codes and prototypes of error-related functions
 * 
 * This file contains error codes and prototypes of error-related functions for the C library. 
 * \note Not all defined error codes are used. Some are reserved for future use.
 */

#ifndef ERROR_H
#define ERROR_H

/* Negative one for unrecognized error code */
#define ERROR_UNKNOWN_ERROR_CODE -1

/* Zero is used to denote success */
#define SUCCESS 0

// 1 - 99: Reserved for error codes returning to python
#define ERROR_SIMULATION_FAILURE 1
#define ERROR_USER_INTERRUPT 2

// 100 - 499: Error codes reserved for gravity_sim

// Initialize system
#define ERROR_INITIALIZE_SYSTEM_NAME_NULL 100
#define ERROR_UNKNOWN_INITIALIZE_SYSTEM_NAME 101
#define ERROR_INITIALIZE_SYSTEM_MEMORY_NOT_NULL 102
#define ERROR_INITIALIZE_SYSTEM_MEMORY_ALLOC 103

// Check if the system has NULL array
#define ERROR_NULL_SYSTEM_POINTER 104

// 300 - 399: Adaptive step size integrator error (general)
#define ERROR_INITIAL_DT_NEGATIVE 300


/* Acceleration */
// 500 - 1999: Reserved for acceleration error

// 500 - 599: Acceleration error (general)
#define ERROR_UNKNOWN_ACCELERATION_METHOD 500
#define ERROR_UNKNOWN_ACCELERATION_CODE 501

// 600 - 609: Pairwise acceleration error
// 610 - 619: Massless acceleration error
#define ERROR_ACCELERATION_MASSLESS_MEMORY_ALLOC 610

// 620 - 699: Barnes-Hut acceleration error
#define ERROR_BARNES_HUT_MORTON_INDICES_MEMORY_ALLOC 620
#define ERROR_BARNES_HUT_RADIX_SORT_MEMORY_ALLOC 621
#define ERROR_BARNES_HUT_OCTREE_MEMORY_ALLOC 622
#define ERROR_BARNES_HUT_SETUP_NODE_MEMORY_REALLOC 623

// 1000 - 1099: CUDA pairwise acceleration error
#define ERROR_CUDA_PAIRWISE_MEMORY_ALLOC 1000
#define ERROR_CUDA_PAIRWISE_MEMCPY_CPU_TO_GPU 1001
#define ERROR_CUDA_PAIRWISE_MEMCPY_GPU_TO_CPU 1002

// 1100 - 1199: CUDA barnes-hut acceleration error
#define ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC 1100
#define ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU 1101
#define ERROR_CUDA_BARNES_HUT_MEMCPY_GPU_TO_CPU 1102

/* Output storage */
// 2000 - 2099: Storing error (general)
#define ERROR_UNKNOWN_STORING_METHOD 2000
#define ERROR_SOL_OUTPUT_MEMORY_ALLOC 2001
#define ERROR_SOL_SIZE_EXCEED_MEMORY_ALLOC 2002
// #define ERROR_STORE_SOLUTION_STEP_METHOD_DISABLED 2003
#define ERROR_STORE_SOLUTION_STEP_UNKNOWN_METHOD 2004
#define ERROR_SOL_OUTPUT_EXTEND_MEMORY_REALLOC 2005


// 2100 - 2199: Flush error
#define ERROR_FLUSH_FILE_OPEN 2100
#define ERROR_FLUSH_FILE_CLOSE 2101
#define ERROR_FLUSH_FILE_CLOSE_IS_NULL 2102
#define ERROR_FLUSH_FILE_CLOSE_NOT_FLUSH_METHOD 2103

// 2500 - 2999: Settings error

/* Integrator */
// 3000 - 3099: Integrator error (general)
#define ERROR_UNKNOWN_INTEGRATOR_METHOD 3000

// 3100 - 3124: Euler integrator error
#define ERROR_EULER_MEMORY_ALLOC 3100

// 3125 - 3149: Euler-Cromer integrator error
#define ERROR_EULER_CROMER_MEMORY_ALLOC 3125

// 3150 - 3174: RK4 integrator error
#define ERROR_RK4_MEMORY_ALLOC 3150

// 3175 - 3199: LeapFrog integrator error
#define ERROR_LEAPFROG_MEMORY_ALLOC 3175

// 3200 - 3299: Embedded RK integrator error
#define ERROR_UNKNOWN_RK_EMBEDDED_METHOD 3200
#define ERROR_RK_EMBEDDED_BUTCHER_TABLEAUS_UNKNOWN_ORDER 3201
#define ERROR_RK_EMBEDDED_BUTCHER_TABLEAUS_MEMORY_ALLOC 3202
#define ERROR_RK_EMBEDDED_INITIAL_DT_MEMORY_ALLOC 3203
#define ERROR_RK_EMBEDDED_INITIAL_DT_NON_POSITIVE 3204
#define ERROR_RK_EMBEDDED_MEMORY_ALLOC 3205

// 3300 - 3399: IAS15 integrator error
#define ERROR_IAS15_INITIAL_DT_MEMORY_ALLOC 3300
#define ERROR_IAS15_INITIAL_DT_NON_POSITIVE 3301
#define ERROR_IAS15_AUX_MEMORY_ALLOC 3302
#define ERROR_IAS15_MEMORY_ALLOC 3303

// 3400 - 3499: WHFast integrator error
#define ERROR_WHFAST_UNKNOWN_ACCELERATION_METHOD 3400
#define ERROR_WHFAST_MEMORY_ALLOC 3401
#define ERROR_WHFAST_KEPLER_AUTO_REMOVE_MEMORY_ALLOC 3402
#define ERROR_WHFAST_ACC_MASSLESS_MEMORY_ALLOC 3403
#define ERROR_WHFAST_STUMPFF_Z_INFINITE 3404
#define ERROR_WHFAST_STUMPFF_Z_NAN 3405

/* Functions prototypes */
/**
 * \brief Print error message to stderr based on the error code
 * 
 * \param error_code Error code
 * 
 * \return None
 */
void print_error_msg(const int error_code);

#endif