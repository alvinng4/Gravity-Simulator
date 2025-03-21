#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration.h"
#include "acceleration_barnes_hut.h"
#include "error.h"
#include "gravity_sim.h"

#ifdef USE_CUDA
    #include "acceleration_cuda.cuh"
#endif

/**
 * \brief Pairwise acceleration computation based on Newton's law of gravitational
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \retval SUCCESS If the computation is successful
 */
IN_FILE int acceleration_pairwise(
    real *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
);

/**
 * \brief Pairwise acceleration computation based on Newton's law of gravitational,
 *        ignoring the contribution of massless particles
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \retval SUCCESS If the computation is successful
 * \retval ERROR_ACCELERATION_MASSLESS_MEMORY_ALLOC If failed to allocate memory
 */
IN_FILE int acceleration_massless(
    real *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
);

WIN32DLL_API int get_acceleration_method_flag(
    const char *__restrict acceleration_method,
    uint *__restrict acceleration_method_flag
)
{
    if (strcmp(acceleration_method, "pairwise") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_PAIRWISE;
        return SUCCESS;
    }
    else if (strcmp(acceleration_method, "massless") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_MASSLESS;
        return SUCCESS;
    }
    else if (strcmp(acceleration_method, "barnes_hut") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_BARNES_HUT;
        return SUCCESS;
    }
#ifdef USE_CUDA
    else if (strcmp(acceleration_method, "pairwise_cuda") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_CUDA_PAIRWISE;
        return SUCCESS;
    }
    else if (strcmp(acceleration_method, "pairwise_cuda_float") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_CUDA_PAIRWISE_FLOAT;
        return SUCCESS;
    }
    else if (strcmp(acceleration_method, "barnes_hut_cuda") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_CUDA_BARNES_HUT;
        return SUCCESS;
    }
    else if (strcmp(acceleration_method, "barnes_hut_cuda_float") == 0)
    {
        *acceleration_method_flag = ACCELERATION_METHOD_CUDA_BARNES_HUT_FLOAT;
        return SUCCESS;
    }
#endif
    else
    {
        return ERROR_UNKNOWN_ACCELERATION_METHOD;
    }
}

WIN32DLL_API int acceleration(
    real *__restrict a,
    const System *__restrict system,
    AccelerationParam *__restrict acceleration_param
)
{
    switch (acceleration_param->acceleration_method_flag_)
    {
        case ACCELERATION_METHOD_PAIRWISE:
            return acceleration_pairwise(a, system, acceleration_param);
        case ACCELERATION_METHOD_MASSLESS:
            return acceleration_massless(a, system, acceleration_param);
        case ACCELERATION_METHOD_BARNES_HUT:
            return acceleration_barnes_hut(a, system, acceleration_param);
#ifdef USE_CUDA
        case ACCELERATION_METHOD_CUDA_PAIRWISE:
            return acceleration_pairwise_cuda(a, system, acceleration_param);
        case ACCELERATION_METHOD_CUDA_PAIRWISE_FLOAT:
            return acceleration_pairwise_cuda_float(a, system, acceleration_param);
        case ACCELERATION_METHOD_CUDA_BARNES_HUT:
            return acceleration_barnes_hut_cuda(a, system, acceleration_param);
        case ACCELERATION_METHOD_CUDA_BARNES_HUT_FLOAT:
            return acceleration_barnes_hut_cuda_float(a, system, acceleration_param);
#endif
        default:
            return ERROR_UNKNOWN_ACCELERATION_CODE;
    }
}

IN_FILE int acceleration_pairwise(
    real *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
)
{
    const int objects_count = system->objects_count;
    const real *x = system->x;
    const real *m = system->m;
    const real G = system->G;
    const real softening_length = acceleration_param->softening_length;

    /* Empty the input array */
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Compute the pairwise acceleration */
    for (int i = 0; i < objects_count; i++)
    {
        const real m_i = m[i];
        for (int j = i + 1; j < objects_count; j++)
        {
            real temp_vec[3];
            real R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[i * 3 + 0] - x[j * 3 + 0];
            R[1] = x[i * 3 + 1] - x[j * 3 + 1];
            R[2] = x[i * 3 + 2] - x[j * 3 + 2];
            const real R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            const real temp_value = G / (R_norm * R_norm * R_norm);
            const real m_j = m[j];
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[i * 3 + 0] -= temp_vec[0] * m_j;
            a[i * 3 + 1] -= temp_vec[1] * m_j;
            a[i * 3 + 2] -= temp_vec[2] * m_j;
            a[j * 3 + 0] += temp_vec[0] * m_i;
            a[j * 3 + 1] += temp_vec[1] * m_i;
            a[j * 3 + 2] += temp_vec[2] * m_i;
        }
    }
    
    return SUCCESS;
}

IN_FILE int acceleration_massless(
    real *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
)
{
    const int objects_count = system->objects_count;
    const real *x = system->x;
    const real *m = system->m;
    const real G = system->G;
    const real softening_length = acceleration_param->softening_length;

    /* Empty the input array */
    for (int i = 0; i < objects_count; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Find the numbers of massive and massless objects */
    int massive_objects_count = 0;
    int massless_objects_count = 0;
    for (int i = 0; i < objects_count; i++)
    {
        if (m[i] != 0.0)
        {
            massive_objects_count++;
        }
        else
        {
            massless_objects_count++;
        }
    }

    /* Find the indices of massive and massless objects */
    int *__restrict massive_indices = malloc(massive_objects_count * sizeof(int));
    int *__restrict massless_indices = malloc(massless_objects_count * sizeof(int));
    massive_objects_count = 0;
    massless_objects_count = 0;

    if (massive_indices == NULL || massless_indices == NULL)
    {
        goto malloc_error;
    }

    for (int i = 0; i < objects_count; i++)
    {
        if (m[i] != 0.0)
        {
            massive_indices[massive_objects_count] = i;
            massive_objects_count++;
        }
        else
        {
            massless_indices[massless_objects_count] = i;
            massless_objects_count++;
        }
    }

    /* Pairwise acceleration calculation for massive objects */
    for (int i = 0; i < massive_objects_count; i++)
    {
        const int idx_i = massive_indices[i];
        const real m_i = m[idx_i];
        for (int j = i + 1; j < massive_objects_count; j++)
        {
            const int idx_j = massive_indices[j];
            const real m_j = m[idx_j];
            real temp_vec[3];
            real R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            const real R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            real temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[idx_i * 3 + 0] -= temp_vec[0] * m_j;
            a[idx_i * 3 + 1] -= temp_vec[1] * m_j;
            a[idx_i * 3 + 2] -= temp_vec[2] * m_j;
            a[idx_j * 3 + 0] += temp_vec[0] * m_i;
            a[idx_j * 3 + 1] += temp_vec[1] * m_i;
            a[idx_j * 3 + 2] += temp_vec[2] * m_i;
        }
    }

    /* Acceleration calculation for massless objects due to massive objects */
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = 0; j < massless_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massless_indices[j];
            real R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            real R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            real temp_value = G / (R_norm * R_norm * R_norm);
            a[idx_j * 3 + 0] += temp_value * R[0] * m[i];
            a[idx_j * 3 + 1] += temp_value * R[1] * m[i];
            a[idx_j * 3 + 2] += temp_value * R[2] * m[i];
        }
    }

    free(massive_indices);
    free(massless_indices);

    return SUCCESS;

malloc_error:
    free(massive_indices);
    free(massless_indices);
    return ERROR_ACCELERATION_MASSLESS_MEMORY_ALLOC;
}
