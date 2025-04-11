/**
 * \file integrator_cosmology_leapfrog.c
 * \brief Leapfrog integrator for cosmological simulations.
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
#include "cosmology.h"
#include "error.h"
#include "integrator.h"
#include "output.h"
#include "progress_bar.h"
#include "settings.h"
#include "system.h"
#include "utils.h"
#include "math_functions.h"


WIN32DLL_API ErrorStatus leapfrog_cosmology(
    CosmologicalSystem *__restrict system,
    IntegratorParam *__restrict integrator_param,
    AccelerationParam *__restrict acceleration_param,
    OutputParam *__restrict output_param,
    SimulationStatus *__restrict simulation_status,
    Settings *__restrict settings,
    const double a_begin,
    const double a_final,
    const int pm_grid_size
)
{
    /* Declare variables */
    ErrorStatus error_status;

    const int num_particles = system->num_particles;
    double *__restrict x = system->x;
    double *__restrict v = system->v;
    const double h0 = system->h0;
    const double omega_m = system->omega_m;
    const double omega_lambda = system->omega_lambda;
    const double G = compute_G(omega_m, h0);
    system->G = G;

    double dt = integrator_param->dt;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *__restrict output_count_ptr = &(output_param->output_count_);
    const double output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    double *__restrict t_ptr = &(simulation_status->t);
    *t_ptr = a_begin;
    int64 *__restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;

    double da;

    /* Allocate memory */
    double *__restrict momentum = malloc(num_particles * 3 * sizeof(double));
    double *__restrict a = malloc(num_particles * 3 * sizeof(double));

    // Check if memory allocation is successful
    if (!momentum || !a)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for arrays");
        goto err_memory;
    }

    /* Get mean background density */
    const double mean_bkg_density = vec_sum(system->m, num_particles) / ((system->box_width * 2) * (system->box_width * 2) * (system->box_width * 2));

    /* Initial output */
    if (is_output && output_param->output_initial)
    {
        error_status = WRAP_TRACEBACK(output_snapshot_cosmology(
            output_param,
            system,
            integrator_param,
            acceleration_param,
            simulation_status,
            settings
        ));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_initial_output;
        }
    }

    /* Set periodic boundary conditions */
    set_periodic_boundary_conditions(system);

    /* Initialize momentum */
    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            momentum[i * 3 + j] = (*t_ptr) * (*t_ptr) * v[i * 3 + j];
        }
    }

    /* Compute initial acceleration */
    error_status = WRAP_TRACEBACK(acceleration_PM(
        a,
        system,
        acceleration_param,
        mean_bkg_density,
        pm_grid_size,
        *t_ptr
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto err_acceleration;
    }

    /* Main Loop */
    int64 total_num_steps = (int64) ceil((a_final - a_begin) / dt);
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        error_status = WRAP_TRACEBACK(start_progress_bar(&progress_bar_param, total_num_steps));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_start_progress_bar;
        }
    }

    simulation_status->dt = dt;
    *num_steps_ptr = 0;
    while (*num_steps_ptr < total_num_steps)
    {
        /* Check dt overshoot */
        if (*t_ptr + dt > a_final)
        {
            dt = a_final - *t_ptr;
        }
        simulation_status->dt = dt;

        /* Kick (p_1/2) */
        da = compute_da(*t_ptr, h0, omega_m, omega_lambda);
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                momentum[i * 3 + j] -= (0.5 * dt) * a[i * 3 + j] / da;
            }
        }
        *t_ptr += 0.5 * dt;

        /* Drift (x_1) */
        da = compute_da(*t_ptr, h0, omega_m, omega_lambda);
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] += dt * momentum[i * 3 + j] / ((*t_ptr) * (*t_ptr) * da);
            }
        }

        /* Set periodic boundary conditions */
        set_periodic_boundary_conditions(system);

        /* Kick (p_1) */
        error_status = WRAP_TRACEBACK(acceleration_PM(
            a,
            system,
            acceleration_param,
            mean_bkg_density,
            pm_grid_size,
            *t_ptr
        ));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_acceleration;
        }

        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                momentum[i * 3 + j] -= (0.5 * dt) * a[i * 3 + j] / da;
            }
        }

        (*num_steps_ptr)++;
        *t_ptr = a_begin + (*num_steps_ptr) * dt;

        /* Store solution */
        if (is_output && *t_ptr >= next_output_time)
        {
            /* Get velocity from momentum */
            for (int i = 0; i < num_particles; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    v[i * 3 + j] = momentum[i * 3 + j] / ((*t_ptr) * (*t_ptr));
                }
            }
            error_status = WRAP_TRACEBACK(output_snapshot_cosmology(
                output_param,
                system,
                integrator_param,
                acceleration_param,
                simulation_status,
                settings
            ));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                goto err_output;
            }

            next_output_time = (*output_count_ptr) * output_interval;
        }

        if (enable_progress_bar)
        {
            update_progress_bar(&progress_bar_param, *num_steps_ptr, false);
        }

        /* Check exit */
        if (*(settings->is_exit_ptr))
        {
            break;
        }
    }

    if (enable_progress_bar)
    {
        update_progress_bar(&progress_bar_param, *num_steps_ptr, true);
    }

    free(momentum);
    free(a);

    return make_success_error_status();

err_output:
err_acceleration:
err_start_progress_bar:
err_initial_output:
err_memory:
    free(momentum);
    free(a);

    return error_status;
}
