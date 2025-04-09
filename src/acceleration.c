/**
 * \file acceleration.c
 * \brief Functions for computing gravitational acceleration
 * 
 * \author Ching-Yin Ng
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration.h"
#include "common.h"
#include "error.h"
#include "math_functions.h"
#include "system.h"
#include "utils.h"

/**
 * \brief Check the acceleration method
 * 
 * \param acceleration_method Acceleration method
 * 
 * \return ErrorStatus
 */
IN_FILE ErrorStatus check_acceleration_method(const int acceleration_method);

/**
 * \brief Compute acceleration with direct pairwise method
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \return ErrorStatus
 */
IN_FILE ErrorStatus acceleration_pairwise(
    double *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
);

/**
 * \brief Compute acceleration with direct pairwise method,
 *        ignoring the contribution of massless particles
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \return ErrorStatus
 */
IN_FILE ErrorStatus acceleration_massless(
    double *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
);


WIN32DLL_API AccelerationParam get_new_acceleration_param(void)
{
    AccelerationParam acceleration_param = {
        .method = ACCELERATION_METHOD_PAIRWISE,
        .opening_angle = 0.5,
        .softening_length = 0.0,
        .max_num_particles_per_leaf = -1
    };
    return acceleration_param;
}

WIN32DLL_API ErrorStatus finalize_acceleration_param(
    AccelerationParam *__restrict acceleration_param
)
{
    ErrorStatus error_status;

    /* Check the acceleration method */
    error_status = WRAP_TRACEBACK(
        check_acceleration_method(acceleration_param->method)
    );
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check the softening length */
    if (acceleration_param->softening_length < 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Softening length is negative. Got: %.3g",
            acceleration_param->softening_length
        );
    }

    /* Check the opening angle */
    if (
        acceleration_param->method == ACCELERATION_METHOD_BARNES_HUT
        && acceleration_param->opening_angle < 0.0
    )
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Opening angle is negative. Got: %.3g",
            acceleration_param->opening_angle
        );
    }

    /* Check the maximum number of particles per leaf */
    if (acceleration_param->method == ACCELERATION_METHOD_BARNES_HUT)
    {
        if (acceleration_param->max_num_particles_per_leaf == -1)
        {
            acceleration_param->max_num_particles_per_leaf = 1;
        }
        else if (acceleration_param->max_num_particles_per_leaf < 1)
        {
            return WRAP_RAISE_ERROR_FMT(
                GRAV_VALUE_ERROR,
                "Maximum number of particles per leaf must be positive. Got: %d",
                acceleration_param->max_num_particles_per_leaf
            );
        }
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus acceleration(
    double *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
)
{
    switch (acceleration_param->method)
    {
        case ACCELERATION_METHOD_PAIRWISE:
            return acceleration_pairwise(a, system, acceleration_param);
        case ACCELERATION_METHOD_MASSLESS:
            return acceleration_massless(a, system, acceleration_param);
        case ACCELERATION_METHOD_BARNES_HUT:
            return acceleration_barnes_hut(a, system, acceleration_param);
        default:
        {
            return WRAP_RAISE_ERROR_FMT(
                GRAV_VALUE_ERROR,
                "Unknown acceleration method. Got: %d",
                acceleration_param->method
            );
        }
    }
}

IN_FILE ErrorStatus check_acceleration_method(const int acceleration_method)
{
    switch (acceleration_method)
    {
        case ACCELERATION_METHOD_PAIRWISE:
        case ACCELERATION_METHOD_MASSLESS:
        case ACCELERATION_METHOD_BARNES_HUT:
        case ACCELERATION_METHOD_PM:
            break;
        default:
        {
            return WRAP_RAISE_ERROR_FMT(
                GRAV_VALUE_ERROR,
                "Unknown acceleration method. Got: %d",
                acceleration_method
            );
        }
    }

    return make_success_error_status();
}

IN_FILE ErrorStatus acceleration_pairwise(
    double *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
)
{
    const int num_particles = system->num_particles;
    const double *x = system->x;
    const double *m = system->m;
    const double G = system->G;
    const double softening_length = acceleration_param->softening_length;

    /* Empty the input array */
    for (int i = 0; i < num_particles; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Compute the pairwise acceleration */
    for (int i = 0; i < num_particles; i++)
    {
        const double m_i = m[i];
        for (int j = i + 1; j < num_particles; j++)
        {
            // Calculate \vec{R} and its norm
            const double R[3] = {
                x[i * 3 + 0] - x[j * 3 + 0],
                x[i * 3 + 1] - x[j * 3 + 1],
                x[i * 3 + 2] - x[j * 3 + 2]
            };
            const double R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            const double temp_value = G / (R_norm * R_norm * R_norm);
            const double m_j = m[j];
            double temp_vec[3] = {
                temp_value * R[0],
                temp_value * R[1],
                temp_value * R[2]
            };
            a[i * 3 + 0] -= temp_vec[0] * m_j;
            a[i * 3 + 1] -= temp_vec[1] * m_j;
            a[i * 3 + 2] -= temp_vec[2] * m_j;
            a[j * 3 + 0] += temp_vec[0] * m_i;
            a[j * 3 + 1] += temp_vec[1] * m_i;
            a[j * 3 + 2] += temp_vec[2] * m_i;
        }
    }

    return make_success_error_status();
}

IN_FILE ErrorStatus acceleration_massless(
    double *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
)
{
    const int num_particles = system->num_particles;
    const double *x = system->x;
    const double *m = system->m;
    const double G = system->G;
    const double softening_length = acceleration_param->softening_length;

    /* Empty the input array */
    for (int i = 0; i < num_particles; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Find the numbers of massive and massless particles */
    int massive_objects_count = 0;
    int massless_objects_count = 0;
    for (int i = 0; i < num_particles; i++)
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

    /* Find the indices of massive and massless particles */
    int *__restrict massive_indices = malloc(massive_objects_count * sizeof(int));
    int *__restrict massless_indices = malloc(massless_objects_count * sizeof(int));
    massive_objects_count = 0;
    massless_objects_count = 0;

    if (massive_indices == NULL || massless_indices == NULL)
    {
        free(massive_indices);
        free(massless_indices);
        return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for massive and massless indices");
    }

    for (int i = 0; i < num_particles; i++)
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

    /* Pairwise acceleration calculation for massive particles */
    for (int i = 0; i < massive_objects_count; i++)
    {
        const int idx_i = massive_indices[i];
        const double m_i = m[idx_i];
        for (int j = i + 1; j < massive_objects_count; j++)
        {
            const int idx_j = massive_indices[j];
            const double m_j = m[idx_j];
            double temp_vec[3];
            double R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            const double R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            double temp_value = G / (R_norm * R_norm * R_norm);
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

    /* Acceleration calculation for massless particles due to massive particles */
    for (int i = 0; i < massive_objects_count; i++)
    {
        for (int j = 0; j < massless_objects_count; j++)
        {
            int idx_i = massive_indices[i];
            int idx_j = massless_indices[j];
            double R[3];

            // Calculate \vec{R} and its norm
            R[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            R[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            R[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];
            double R_norm = sqrt(
                R[0] * R[0] + 
                R[1] * R[1] + 
                R[2] * R[2] +
                softening_length * softening_length
            );

            // Calculate the acceleration
            double temp_value = G / (R_norm * R_norm * R_norm);
            a[idx_j * 3 + 0] += temp_value * R[0] * m[i];
            a[idx_j * 3 + 1] += temp_value * R[1] * m[i];
            a[idx_j * 3 + 2] += temp_value * R[2] * m[i];
        }
    }

    free(massive_indices);
    free(massless_indices);

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus benchmark_acceleration(
    const System *__restrict system,
    const AccelerationParam *acceleration_params,
    const int num_acceleration_params,
    const int *__restrict num_times_acceleration_param    
)
{
    ErrorStatus error_status;

    double *__restrict reference_a = malloc(
        system->num_particles * 3 * sizeof(double)
    );
    double *__restrict a = malloc(
        system->num_particles * 3 * sizeof(double)
    );
    if (!reference_a || !a)
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for acceleration arrays"
        );
        goto err_malloc;
    }

    fputs("Benchmarking acceleration...\n", stdout);

    for (int i = 0; i < num_acceleration_params; i++)
    {
        const AccelerationParam *acceleration_param = &(acceleration_params[i]);
        const int num_times = num_times_acceleration_param[i];

        if (num_times <= 0)
        {
            printf("Test %d:    Skipped since num_times: %d <= 0\n\n", i, num_times);
            continue;
        }

        double *__restrict run_time = calloc(num_times, sizeof(double));
        double mae = 0.0;

        if (!run_time)
        {
            free(run_time);
            error_status = WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to allocate memory for runtime array"
            );
            goto err_malloc;
        }

        for (int j = 0; j < num_times; j++)
        {
            if (i == 0 && j == 0)
            {
                double start_time = grav_get_current_time();
                error_status = WRAP_TRACEBACK(acceleration(
                    reference_a,
                    system,
                    acceleration_param
                ));
                if (error_status.return_code != GRAV_SUCCESS)
                {
                    return error_status;
                }
                double end_time = grav_get_current_time();
                run_time[j] += (end_time - start_time);
            }
            else
            {
                double start_time = grav_get_current_time();
                error_status = WRAP_TRACEBACK(acceleration(
                    a,
                    system,
                    acceleration_param
                ));
                if (error_status.return_code != GRAV_SUCCESS)
                {
                    return error_status;
                }
                double end_time = grav_get_current_time();
                run_time[j] += (end_time - start_time);
            }

            // Calculate the MAE
            if (i != 0 && j == 0)
            {
                for (int k = 0; k < system->num_particles; k++)
                {
                    const double diff[3] = {
                        reference_a[k * 3 + 0] - a[k * 3 + 0],
                        reference_a[k * 3 + 1] - a[k * 3 + 1],
                        reference_a[k * 3 + 2] - a[k * 3 + 2]
                    };
                    mae += fabs(diff[0]) + fabs(diff[1]) + fabs(diff[2]);
                }
                mae /= system->num_particles;
            }
        }

        printf("Test %d:", i);
        switch(acceleration_param->method)
        {
            case ACCELERATION_METHOD_PAIRWISE:
                fputs("    Method: Pairwise\n", stdout);
                break;
            case ACCELERATION_METHOD_MASSLESS:
                fputs("    Method: Massless\n", stdout);
                break;
            case ACCELERATION_METHOD_BARNES_HUT:
                fputs("    Method: Barnes-Hut\n", stdout);
                break;
            default:
                error_status = WRAP_RAISE_ERROR_FMT(
                    GRAV_VALUE_ERROR,
                    "Unknown acceleration method. Got: %d",
                    acceleration_param->method
                );
                goto err_unknown_acceleration_method;
        }
        
        printf("    Number of times: %d\n", num_times);
        printf("    Avg time: %.3g (+- %.3g) s\n", compute_mean(run_time, num_times), compute_std(run_time, num_times, 1));
        printf("    MAE: %.3g\n", mae);
        printf("\n");

        free(run_time);
    }

    free(reference_a);
    free(a);

    return make_success_error_status();

err_unknown_acceleration_method:
err_malloc:
    free(reference_a);
    free(a);
    return error_status;
}
