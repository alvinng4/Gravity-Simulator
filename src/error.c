#include <stdio.h>

#include "error.h"
#include "gravity_sim.h"

IN_FILE int get_error_msg(const int error_code, char **error_msg);

WIN32DLL_API void print_error_msg(const int error_code)
{
    char *error_msg = NULL;
    int return_code = get_error_msg(error_code, &error_msg);
    if (return_code == SUCCESS)
    {
        fprintf(stderr, "%s", error_msg);
    }
    else
    {
        fprintf(stderr, "C library error: error detected but error code is not recognized. (Error Code: %d)\n", error_code);
    }
}

IN_FILE int get_error_msg(const int error_code, char **error_msg)
{
    switch (error_code)
    {
        /* Error codes returning to python */
        case ERROR_USER_INTERRUPT:
            *error_msg = "C library error: User interrupt detected.\n";
            return SUCCESS;
        /* Error codes reserved for gravity_sim */
        /* Initialize system */
        case ERROR_INITIALIZE_SYSTEM_NAME_NULL:
            *error_msg = "C library error: System name is NULL in initialize_system().\n";
            return SUCCESS;
        case ERROR_UNKNOWN_INITIALIZE_SYSTEM_NAME:
            *error_msg = "C library error: System name not recognized in initialize_system().\n";
            return SUCCESS;
        case ERROR_INITIALIZE_SYSTEM_MEMORY_NOT_NULL:
            *error_msg = "C library error: Initial memory is not NULL in initialize_system().\n";
            return SUCCESS;
        case ERROR_INITIALIZE_SYSTEM_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in initialize_system().\n";
            return SUCCESS;
        // Check if the system has NULL array
        case ERROR_NULL_SYSTEM_POINTER:
            *error_msg = "C library error: System has NULL pointer in launch_simulation().\n";
            return SUCCESS;

        /* Adaptive step size integrator error (general) */
        case ERROR_INITIAL_DT_NEGATIVE:
            *error_msg = "C library error: Initial time step is negative in _launch_simulation().\n";
            return SUCCESS;

        /* Acceleration */
        case ERROR_UNKNOWN_ACCELERATION_METHOD:
            *error_msg = "C library error: Acceleration method not recognized in get_acceleration_method_flag().\n";
            return SUCCESS;
        case ERROR_UNKNOWN_ACCELERATION_CODE:
            *error_msg = "C library error: Acceleration code not recognized in acceleration().\n";
            return SUCCESS;

        // Massless acceleration error
        case ERROR_ACCELERATION_MASSLESS_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in acceleration_massless().\n";
            return SUCCESS;

        // Barnes-Hut acceleration error
        case ERROR_BARNES_HUT_MORTON_INDICES_MEMORY_ALLOC:
            *error_msg = "C library error: Failed to allocate memory for morton indices in acceleration_barnes_hut().\n";
            return SUCCESS;
        case ERROR_BARNES_HUT_RADIX_SORT_MEMORY_ALLOC:
            *error_msg = "C library error: Failed to allocate memory for radix sort in acceleration_barnes_hut().\n";
            return SUCCESS;
        case ERROR_BARNES_HUT_OCTREE_MEMORY_ALLOC:
            *error_msg = "C library error: Failed to allocate memory for octree in acceleration_barnes_hut().\n";
            return SUCCESS;
        case ERROR_BARNES_HUT_SETUP_NODE_MEMORY_REALLOC:
            *error_msg = "C library error: Failed to reallocate memory for octree internal nodes in acceleration_barnes_hut() for _setup_node().\n";
            return SUCCESS;

        // CUDA pairwise acceleration error
        case ERROR_CUDA_PAIRWISE_MEMORY_ALLOC:
            *error_msg = "C library error: Failed to allocate GPU memory in acceleration_cuda_pairwise().\n";   
            return SUCCESS;
        case ERROR_CUDA_PAIRWISE_MEMCPY_CPU_TO_GPU:
            *error_msg = "C library error: Failed to copy data from CPU to GPU in acceleration_cuda_pairwise().\n";
            return SUCCESS; 
        case ERROR_CUDA_PAIRWISE_MEMCPY_GPU_TO_CPU:
            *error_msg = "C library error: Failed to copy data from GPU to CPU in acceleration_cuda_pairwise().\n";
            return SUCCESS;

        // CUDA Barnes-Hut acceleration error
        case ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC:
            *error_msg = "C library error: Failed to allocate GPU memory in _compute_acceleration() for acceleration_barnes_hut_cuda().\n";   
            return SUCCESS;
        case ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU:
            *error_msg = "C library error: Failed to copy data from CPU to GPU in _compute_acceleration() for acceleration_barnes_hut_cuda().\n";
            return SUCCESS; 
        case ERROR_CUDA_BARNES_HUT_MEMCPY_GPU_TO_CPU:
            *error_msg = "C library error: Failed to copy data from GPU to CPU in _compute_acceleration() for acceleration_barnes_hut_cuda().\n";
            return SUCCESS;

        /* Solution storing (general) */
        // Storing error (general)
        case ERROR_UNKNOWN_STORING_METHOD:
            *error_msg = "C library error: Storing method not recognized in get_storing_method_flag().\n";
            return SUCCESS;
        case ERROR_SOL_OUTPUT_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in allocate_solutions_memory().\n";
            return SUCCESS;
        case ERROR_SOL_SIZE_EXCEED_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in extend_solutions_memory().\n";
            return SUCCESS;
        case ERROR_STORE_SOLUTION_STEP_UNKNOWN_METHOD:
            *error_msg = "C library error: Storing method not recognized in store_solution_step().\n";
            return SUCCESS;
        case ERROR_SOL_OUTPUT_EXTEND_MEMORY_REALLOC:
            *error_msg = "C library error: Memory reallocation failed in extend_solutions_memory().\n";
            return SUCCESS;
        
        // Flush error
        case ERROR_FLUSH_FILE_OPEN:
            *error_msg = "C library error: Failed to open flush file in open_flush_file().\n";
            return SUCCESS;
        case ERROR_FLUSH_FILE_CLOSE:
            *error_msg = "C library error: Failed to close flush file in close_flush_file().\n";
            return SUCCESS;
        case ERROR_FLUSH_FILE_CLOSE_IS_NULL:
            *error_msg = "C library error: Flush file is NULL in close_flush_file().\n";
            return SUCCESS;
        case ERROR_FLUSH_FILE_CLOSE_NOT_FLUSH_METHOD:
            *error_msg = "C library error: close_flush_file() is called but storing method is not flush.\n";
            return SUCCESS;

        /* Integrator */
        case ERROR_UNKNOWN_INTEGRATOR_METHOD:
            *error_msg = "C library error: Integrator method not recognized.\n";
            return SUCCESS;

        // Euler integrator
        case ERROR_EULER_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in euler().\n";
            return SUCCESS;

        // Euler-Cromer integrator
        case ERROR_EULER_CROMER_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in euler_cromer().\n";
            return SUCCESS;

        // RK4 integrator
        case ERROR_RK4_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in rk4().\n";
            return SUCCESS;

        // LeapFrog integrator
        case ERROR_LEAPFROG_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in leapfrog().\n";
            return SUCCESS;

        // Embedded RK integrator
        case ERROR_UNKNOWN_RK_EMBEDDED_METHOD:
            *error_msg = "C library error: Embedded Runge-Kutta method not recognized in get_rk_embedded_order().\n";
            return SUCCESS;
        case ERROR_RK_EMBEDDED_BUTCHER_TABLEAUS_UNKNOWN_ORDER:
            *error_msg = "C library error: Order not recognized in rk_embedded_butcher_tableaus().\n";
            return SUCCESS;
        case ERROR_RK_EMBEDDED_BUTCHER_TABLEAUS_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in rk_embedded_butcher_tableaus().\n";
            return SUCCESS;
        case ERROR_RK_EMBEDDED_INITIAL_DT_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in rk_embedded_initial_dt().\n";
            return SUCCESS;
        case ERROR_RK_EMBEDDED_INITIAL_DT_NON_POSITIVE:
            *error_msg = "C library error: Initial time step is non-positive in rk_embedded_initial_dt().\n";
            return SUCCESS;
        case ERROR_RK_EMBEDDED_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in rk_embedded().\n";
            return SUCCESS;

        // IAS15 integrator
        case ERROR_IAS15_INITIAL_DT_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in ias15 _initial_dt().\n";
            return SUCCESS;
        case ERROR_IAS15_INITIAL_DT_NON_POSITIVE:
            *error_msg = "C library error: Initial time step is non-positive in ias15().\n";
            return SUCCESS;
        case ERROR_IAS15_AUX_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation for auxiliary variables failed in ias15.\n";
            return SUCCESS;
        case ERROR_IAS15_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in ias15().\n";
            return SUCCESS;

        // WHFast integrator
        case ERROR_WHFAST_UNKNOWN_ACCELERATION_METHOD:
            *error_msg = "C library error: Acceleration method not recognized in whfast.\n";
            return SUCCESS;
        case ERROR_WHFAST_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in whfast().\n";
            return SUCCESS;
        case ERROR_WHFAST_KEPLER_AUTO_REMOVE_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in whfast for kepler auto remove.\n";
            return SUCCESS;
        case ERROR_WHFAST_ACC_MASSLESS_MEMORY_ALLOC:
            *error_msg = "C library error: Memory allocation failed in whfast for massless acceleration.\n";
            return SUCCESS;
        case ERROR_WHFAST_STUMPFF_Z_INFINITE:
            *error_msg = "C library error: infinite value detected in whfast stumpff function.\n";
            return SUCCESS;
        case ERROR_WHFAST_STUMPFF_Z_NAN:
            *error_msg = "C library error: NaN value detected in whfast stumpff function.\n";
            return SUCCESS;

        default:
            return ERROR_UNKNOWN_ERROR_CODE;
    }
}
