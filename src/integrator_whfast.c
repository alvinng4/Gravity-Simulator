/**
 * \file integrator_whfast.c
 * \brief Function definitions for WHFast integrators
 * 
 * \cite J. Roa, et al. Moving Planets Around: An Introduction to
 *   N-Body Simulations Applied to Exoplanetary Systems*, MIT
 *   Press, 2020
 * 
 * \author Ching-Yin Ng
 * \date March 2025
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
#include "output.h"
#include "progress_bar.h"
#include "settings.h"
#include "system.h"

#define WHFAST_KEPLER_TOL 1e-12
#define WHFAST_KEPLER_MAX_ITER 500

/**
 * \brief Compute the velocity kick
 * 
 * \param[out] jacobi_v Array of Jacobi velocity vectors
 * \param[in] objects_count Number of objects in the system
 * \param[in] a Array of acceleration vectors
 * \param[in] dt Time step of the system
 */
IN_FILE void whfast_kick(
    double *__restrict jacobi_v,
    const int objects_count,
    const double *__restrict a,
    const double dt
);

/** 
 * \brief Compute the position drift
 * 
 * \param[out] jacobi_x Array of Jacobi position vectors
 * \param[out] jacobi_v Array of Jacobi velocity vectors
 * \param[in] system Pointer to the gravitational system
 * \param[in] eta Array of cumulative masses
 * \param[in] dt Time step of the system
 * \param[in] verbose Verbosity level
 *
 * \return ErrorStatus
 * 
 * \exception GRAV_VALUE_ERROR If the input value to the stumpff function is infinite or NaN
 */
IN_FILE ErrorStatus whfast_drift(
    double *__restrict jacobi_x,
    double *__restrict jacobi_v,
    const System *__restrict system,
    const double *__restrict eta,
    const double dt,
    const int verbose
);

/**
 * \brief Transform Cartesian coordinates to Jacobi coordinates
 * 
 * \param jacobi_x Array of Jacobi position vectors to be stored
 * \param jacobi_v Array of Jacobi velocity vectors to be stored
 * \param system Pointer to the gravitational system
 * \param eta Array of cumulative masses
 */
IN_FILE void cartesian_to_jacobi(
    double *__restrict jacobi_x,
    double *__restrict jacobi_v,
    const System *__restrict system,
    const double *__restrict eta
);

/**
 * \brief Transform Jacobi coordinates to Cartesian coordinates
 * 
 * \param system Pointer to the gravitational system
 * \param jacobi_x Array of Jacobi position vectors
 * \param jacobi_v Array of Jacobi velocity vectors
 * \param eta Array of cumulative masses
 */
IN_FILE void jacobi_to_cartesian(
    System *__restrict system,
    const double *__restrict jacobi_x,
    const double *__restrict jacobi_v,
    const double *__restrict eta
);

/**
 * \brief Compute the Stumpff functions c0, c1, c2, and c3 for a given argument z
 * 
 * \param[out] c0 Pointer to store c0
 * \param[out] c1 Pointer to store c1
 * \param[out] c2 Pointer to store c2
 * \param[out] c3 Pointer to store c3
 * \param[in] z Input value
 */
IN_FILE void stumpff_functions(
    double *__restrict c0,
    double *__restrict c1,
    double *__restrict c2,
    double *__restrict c3,
    double z
);

/**
 * \brief Compute the acceleration for the WHFast integrator
 * 
 * \param[out] a Array of acceleration vectors to be stored
 * \param[in] system Pointer to the gravitational system
 * \param[in] jacobi_x Array of Jacobi position vectors
 * \param[in] eta Array of cumulative masses
 * \param[in] acceleration_param Pointer to acceleration parameters
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_VALUE_ERROR If the acceleration method is not supported
 * \exception Errors from the acceleration functions if any error occurs
 */
IN_FILE ErrorStatus whfast_acceleration(
    double *__restrict a,
    const System *__restrict system,
    const double *__restrict jacobi_x,
    const double *__restrict eta,
    const AccelerationParam *__restrict acceleration_param
);

/**
 * \brief Direct pairwise acceleration function for WHFast integrator
 * 
 * \details This is a brute-force pairwise calculation
 *          of gravitational acceleration between all objects,
 *          which is O(n^2) complexity.
 * 
 * \param[out] a Array of acceleration vectors to be stored
 * \param[in] system Pointer to the gravitational system
 * \param[in] jacobi_x Array of Jacobi position vectors
 * \param[in] eta Array of cumulative masses
 * \param[in] acceleration_param Pointer to acceleration parameters
 * 
 * \return ErrorStatus
 */
IN_FILE ErrorStatus whfast_acceleration_pairwise(
    double *__restrict a,
    const System *__restrict system,
    const double *__restrict jacobi_x,
    const double *__restrict eta,
    const AccelerationParam *__restrict acceleration_param
);

/**
 * \brief Acceleration function for WHFast integrator,
 *        separating massive and massless objects
 * 
 * \details This function calculates the gravitational acceleration
 *          between massive and massless objects separately.
 *          This is an O(m^2 + mn) complexity calculation,
 *          where m and n are the number of massive and massless 
 *          objects, respectively.
 * 
 * \param[out] a Array of acceleration vectors to be stored
 * \param[in] system Pointer to the gravitational system
 * \param[in] jacobi_x Array of Jacobi position vectors
 * \param[in] eta Array of cumulative masses
 * \param[in] acceleration_param Pointer to acceleration parameters
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR If memory allocation failed
 */

IN_FILE ErrorStatus whfast_acceleration_massless(
    double *__restrict a,
    const System *__restrict system,
    const double *__restrict jacobi_x,
    const double *__restrict eta,
    const AccelerationParam *__restrict acceleration_param
);

WIN32DLL_API ErrorStatus whfast(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double tf
)
{
    ErrorStatus error_status;

    const int objects_count = system->objects_count;
    double *__restrict m = system->m;

    double dt = integrator_param->dt;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *__restrict output_count_ptr = &(output_param->output_count_);
    const int output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    double *__restrict t_ptr = &(simulation_status->t);
    int64 *__restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;
    const int verbose = settings->verbose;

    /* Allocate memory */
    double *__restrict jacobi_x = calloc(objects_count * 3, sizeof(double));
    double *__restrict jacobi_v = malloc(objects_count * 3 * sizeof(double));
    double *__restrict temp_jacobi_v = malloc(objects_count * 3 * sizeof(double));
    double *__restrict a = malloc(objects_count * 3 * sizeof(double));
    double *__restrict eta = malloc(objects_count * sizeof(double));

    // Check if memory allocation is successful
    if (!jacobi_x || !jacobi_v || !temp_jacobi_v || !a || !eta)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for arrays");
        goto err_memory;
    }

    /* Initial output */
    const int initial_output_offset = (output_param->output_initial) ? 1 : 0;
    if (is_output && output_param->output_initial)
    {
        error_status = WRAP_TRACEBACK(output_snapshot(
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

    /* Initialization */
    eta[0] = m[0];
    for (int i = 1; i < objects_count; i++)
    {
        eta[i] = eta[i - 1] + m[i];
    }
    cartesian_to_jacobi(jacobi_x, jacobi_v, system, eta);
    error_status = whfast_acceleration(a, system, jacobi_x, eta, acceleration_param);
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto err_acceleration;
    }
    whfast_kick(jacobi_v, objects_count, a, 0.5 * dt);

    /* Main Loop */
    int64 total_num_steps = (int64) ceil(tf / dt);
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        start_progress_bar(&progress_bar_param, total_num_steps);
    }

    *t_ptr = 0.0;
    simulation_status->dt = dt;
    *num_steps_ptr = 0;
    while (*num_steps_ptr < total_num_steps)
    {
        /* Check dt overshoot */
        if (*t_ptr + dt > tf)
        {
            dt = tf - *t_ptr;
        }
        simulation_status->dt = dt;

        error_status = whfast_drift(
            jacobi_x,
            jacobi_v,
            system,
            eta,
            dt,
            verbose
        );
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_drift;
        }

        jacobi_to_cartesian(system, jacobi_x, jacobi_v, eta);
        error_status = whfast_acceleration(a, system, jacobi_x, eta, acceleration_param);
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_acceleration;
        }
        whfast_kick(jacobi_v, objects_count, a, dt);

        (*num_steps_ptr)++;
        *t_ptr = (*num_steps_ptr) * dt;

        /* Store solution */
        if (is_output && *t_ptr >= next_output_time)
        {
            // Get v_1 from v_1+1/2
            memcpy(temp_jacobi_v, jacobi_v, objects_count * 3 * sizeof(double));
            whfast_kick(temp_jacobi_v, objects_count, a, -0.5 * dt);
            jacobi_to_cartesian(system, jacobi_x, temp_jacobi_v, eta);
            error_status = WRAP_TRACEBACK(output_snapshot(
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

            next_output_time = (*output_count_ptr - initial_output_offset) * output_interval;
        }

        if (enable_progress_bar)
        {
            update_progress_bar(&progress_bar_param, *num_steps_ptr, false);
        }

        /* Check exit */
        if (settings->is_exit)
        {
            break;
        }
    }

    if (enable_progress_bar)
    {
        update_progress_bar(&progress_bar_param, *num_steps_ptr, true);
    }

    free(jacobi_x);
    free(jacobi_v);
    free(temp_jacobi_v);
    free(a);
    free(eta);

    return make_success_error_status();

err_output:
err_acceleration:
err_initial_output:
err_drift:
err_memory:
    free(eta);
    free(a);
    free(temp_jacobi_v);
    free(jacobi_v);
    free(jacobi_x);

    return error_status;
}

IN_FILE void whfast_kick(
    double *__restrict jacobi_v,
    const int objects_count,
    const double *__restrict a,
    const double dt
)
{
    for (int i = 0; i < objects_count; i++)
    {
        jacobi_v[i * 3 + 0] += a[i * 3 + 0] * dt;
        jacobi_v[i * 3 + 1] += a[i * 3 + 1] * dt;
        jacobi_v[i * 3 + 2] += a[i * 3 + 2] * dt;
    }
}

IN_FILE ErrorStatus whfast_drift(
    double *__restrict jacobi_x,
    double *__restrict jacobi_v,
    const System *__restrict system,
    const double *__restrict eta,
    const double dt,
    const int verbose
)
{
    const int objects_count = system->objects_count;
    const int *__restrict particle_ids = system->particle_ids;
    const double *__restrict m = system->m;
    const double G = system->G;
    for (int i = 1; i < objects_count; i++)
    {
        const double gm = G * m[0] * eta[i] / eta[i - 1];
        const double x[3] = {jacobi_x[i * 3], jacobi_x[i * 3 + 1], jacobi_x[i * 3 + 2]};
        const double v[3] = {jacobi_v[i * 3], jacobi_v[i * 3 + 1], jacobi_v[i * 3 + 2]};

        const double x_norm = vec_norm_3d(x);
        const double v_norm = vec_norm_3d(v);

        // Radial velocity
        const double radial_v = vec_dot_3d(x, v) / x_norm; 

        const double alpha = 2.0 * gm / x_norm - (v_norm * v_norm);

        /* Solve Kepler's equation with Newton-Raphson method */

        // Initial guess
        double s = dt / x_norm;

        // Solve Kepler's equation
        double c0 = 0.0;
        double c1 = 0.0;
        double c2 = 0.0;
        double c3 = 0.0;
        bool is_converged = false;

        for (int j = 0; j < WHFAST_KEPLER_MAX_ITER; j++)
        {
            // Compute Stumpff functions
            const double z = alpha * (s * s);
            if (isnan(z) || isinf(z))
            {
                return WRAP_RAISE_ERROR(
                    GRAV_VALUE_ERROR,
                    "Input value to Stumpff functions is NaN or Inf"
                );
            }
            stumpff_functions(&c0, &c1, &c2, &c3, z);

            // Evaluate Kepler's equation and its derivative
            const double F = (
                x_norm * s * c1
                + x_norm * radial_v * (s * s) * c2
                + gm * (s * s * s) * c3
                - dt
            );
            const double dF = (
                x_norm * c0
                + x_norm * radial_v * s * c1
                + gm * (s * s) * c2
            );

            // Advance step
            const double ds = -F / dF;
            s += ds;

            // Check convergence
            if (fabs(ds) < WHFAST_KEPLER_TOL)
            {
                is_converged = true;
                break;
            }
        }

        // The raidal distance is equal to the derivative of F
        // double r = dF
        const double r = x_norm * c0 + x_norm * radial_v * s * c1 + gm * (s * s) * c2;

        if (!is_converged)
        {
            const double error = (
                x_norm * s * c1
                + x_norm * radial_v * (s * s) * c2
                + gm * (s * s * s) * c3
                - dt
            ) / r;

            /* Print warning message */
            if (verbose >= GRAV_VERBOSITY_IGNORE_INFO)
            {
                const int warning_msg_len = (
                    strlen("Kepler's equation did not converge. Particle id: , error = \n")
                    + snprintf(NULL, 0, "%d", particle_ids[i])
                    + snprintf(NULL, 0, "%23.15g", error)
                    + 1
                );
                char *warning_msg = malloc(warning_msg_len);
                if (!warning_msg)
                {
                    return WRAP_RAISE_ERROR(
                        GRAV_MEMORY_ERROR,
                        "Kepler's equation did not converge and failed to allocate memory for warning message"
                    );
                }

                const int actual_warning_msg_len = snprintf(
                    warning_msg,
                    warning_msg_len,
                    "Warning: Kepler's equation did not converge. "\
                    "Particle id: %d, error = %23.15g\n",
                    particle_ids[i], error
                );

                if (actual_warning_msg_len < 0)
                {
                    free(warning_msg);
                    return WRAP_RAISE_ERROR(
                        GRAV_UNKNOWN_ERROR,
                        "Kepler's equation did not converge and failed to generate warning message"
                    );
                }
                else if (actual_warning_msg_len >= warning_msg_len)
                {
                    free(warning_msg);
                    return WRAP_RAISE_ERROR(
                        GRAV_UNKNOWN_ERROR,
                        "Kepler's equation did not converge and warning message are truncated"
                    );
                }

                WRAP_RAISE_WARNING(warning_msg);
                free(warning_msg);
            }

            // if (kepler_auto_remove && ((fabs(error) > kepler_auto_remove_tol) || isnan(error)))
            // {
            //     kepler_failed_bool_array[i] = true;
            //     *kepler_failed_flag = true;
            // }
        }

        /* Evaluate f and g functions, together with their derivatives */
        const double f = 1.0 - gm * (s * s) * c2 / x_norm;
        const double g = dt - gm * (s * s * s) * c3;

        const double df = -gm * s * c1 / (r * x_norm);
        const double dg = 1.0 - gm * (s * s) * c2 / r; 

        /* Compute position and velocity vectors */
        for (int j = 0; j < 3; j++)
        {
            jacobi_x[i * 3 + j] = f * x[j] + g * v[j];
            jacobi_v[i * 3 + j] = df * x[j] + dg * v[j];
        }
    }

    return make_success_error_status();
}

IN_FILE void cartesian_to_jacobi(
    double *__restrict jacobi_x,
    double *__restrict jacobi_v,
    const System *__restrict system,
    const double *__restrict eta
)
{
    double x_cm[3];
    double v_cm[3];

    const int objects_count = system->objects_count;
    const double *__restrict x = system->x;
    const double *__restrict v = system->v;
    const double *__restrict m = system->m;

    x_cm[0] = m[0] * x[0];
    x_cm[1] = m[0] * x[1];
    x_cm[2] = m[0] * x[2];

    v_cm[0] = m[0] * v[0];
    v_cm[1] = m[0] * v[1];
    v_cm[2] = m[0] * v[2];

    for (int i = 1; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            jacobi_x[i * 3 + j] = x[i * 3 + j] - x_cm[j] / eta[i - 1];
            jacobi_v[i * 3 + j] = v[i * 3 + j] - v_cm[j] / eta[i - 1];
        
            x_cm[j] = x_cm[j] * (1.0 + m[i] / eta[i - 1]) + m[i] * jacobi_x[i * 3 + j];
            v_cm[j] = v_cm[j] * (1.0 + m[i] / eta[i - 1]) + m[i] * jacobi_v[i * 3 + j];
        }
    }

    jacobi_x[0] = x_cm[0] / eta[objects_count - 1];
    jacobi_x[1] = x_cm[1] / eta[objects_count - 1];
    jacobi_x[2] = x_cm[2] / eta[objects_count - 1];

    jacobi_v[0] = v_cm[0] / eta[objects_count - 1];
    jacobi_v[1] = v_cm[1] / eta[objects_count - 1];
    jacobi_v[2] = v_cm[2] / eta[objects_count - 1];
}

IN_FILE void jacobi_to_cartesian(
    System *__restrict system,
    const double *__restrict jacobi_x,
    const double *__restrict jacobi_v,
    const double *__restrict eta
)
{
    double x_cm[3];
    double v_cm[3];

    double *__restrict x = system->x;
    double *__restrict v = system->v;

    const int objects_count = system->objects_count;
    const double *__restrict m = system->m;

    x_cm[0] = eta[objects_count - 1] * jacobi_x[0];
    x_cm[1] = eta[objects_count - 1] * jacobi_x[1];
    x_cm[2] = eta[objects_count - 1] * jacobi_x[2];

    v_cm[0] = eta[objects_count - 1] * jacobi_v[0];
    v_cm[1] = eta[objects_count - 1] * jacobi_v[1];
    v_cm[2] = eta[objects_count - 1] * jacobi_v[2];

    for (int i = (objects_count - 1); i > 0; i--)
    {
        for (int j = 0; j < 3; j++)
        {
            x_cm[j] = (x_cm[j] - m[i] * jacobi_x[i * 3 + j]) / eta[i];
            v_cm[j] = (v_cm[j] - m[i] * jacobi_v[i * 3 + j]) / eta[i];

            x[i * 3 + j] = jacobi_x[i * 3 + j] + x_cm[j];
            v[i * 3 + j] = jacobi_v[i * 3 + j] + v_cm[j];

            x_cm[j] = eta[i - 1] * x_cm[j];
            v_cm[j] = eta[i - 1] * v_cm[j];
        }
    }
    
    x[0] = x_cm[0] / m[0];
    x[1] = x_cm[1] / m[0];
    x[2] = x_cm[2] / m[0];

    v[0] = v_cm[0] / m[0];
    v[1] = v_cm[1] / m[0];
    v[2] = v_cm[2] / m[0];
}

IN_FILE void stumpff_functions(
    double *__restrict c0,
    double *__restrict c1,
    double *__restrict c2,
    double *__restrict c3,
    double z
)
{
    /* Reduce the argument */
    int n = 0;
    while (fabs(z) > 0.1)
    {
        z /= 4.0;
        n++;
    }

    /* Compute stumpff functions */
    double temp_c3 = (
        1.0 - z / 20.0 * (1.0 - z / 42.0 * (1.0 - z / 72.0 * (1.0 - z / 110.0 \
        * (1.0 - z / 156.0 * (1.0 - z / 210.0)))))
    ) / 6.0;
    double temp_c2 = (
        1.0 - z / 12.0 * (1.0 - z / 30.0 * (1.0 - z / 56.0 * (1.0 - z / 90.0 \
        * (1.0 - z / 132.0 * (1.0 - z / 182.0)))))
    ) / 2.0;
    double temp_c1 = 1.0 - z * temp_c3;
    double temp_c0 = 1.0 - z * temp_c2;

    /* Half-angle formulae to recover the actual argument */
    while (n > 0)
    {
        temp_c3 = (temp_c2 + temp_c0 * temp_c3) / 4.0;
        temp_c2 = (temp_c1 * temp_c1) / 2.0;
        temp_c1 = temp_c0 * temp_c1;
        temp_c0 = (2.0 * temp_c0 * temp_c0) - 1.0;
        n--;
    }
    *c3 = temp_c3;
    *c2 = temp_c2;
    *c1 = temp_c1;
    *c0 = temp_c0;
}

IN_FILE ErrorStatus whfast_acceleration(
    double *__restrict a,
    const System *system,
    const double *__restrict jacobi_x,
    const double *__restrict eta,
    const AccelerationParam *acceleration_param
)
{
    switch (acceleration_param->method)
    {
        case ACCELERATION_METHOD_PAIRWISE:
            return whfast_acceleration_pairwise(a, system, jacobi_x, eta, acceleration_param);
        case ACCELERATION_METHOD_MASSLESS:
            return whfast_acceleration_massless(a, system, jacobi_x, eta, acceleration_param);
        default:
            return WRAP_RAISE_ERROR(
                GRAV_VALUE_ERROR,
                "Invalid acceleration method for WHFast integrator. Only pairwise and massless methods are supported."
            );
    }
}

IN_FILE ErrorStatus whfast_acceleration_pairwise(
    double *__restrict a,
    const System *__restrict system,
    const double *__restrict jacobi_x,
    const double *__restrict eta,
    const AccelerationParam *__restrict acceleration_param
)
{
    const int objects_count = system->objects_count;
    const double *__restrict x = system->x;
    const double *__restrict m = system->m;
    const double G = system->G;

    const double softening_length = acceleration_param->softening_length;

    double aux[3];
    double temp_vec[3];
    double temp_vec_norm;
    double temp_vec_norm_cube;
    double temp_jacobi_norm;
    double temp_jacobi_norm_cube;
    double softening_length_cube = softening_length * softening_length * softening_length;
    for (int i = 1; i < objects_count; i++)
    {
        // Calculate x_0i
        temp_vec[0] = x[i * 3 + 0] - x[0];
        temp_vec[1] = x[i * 3 + 1] - x[1];
        temp_vec[2] = x[i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm_3d(temp_vec);
        temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;
        temp_jacobi_norm = vec_norm_3d(&jacobi_x[i * 3]);
        temp_jacobi_norm_cube = (temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm) + softening_length_cube;
        for (int j = 0; j < 3; j++)
        {
            a[i * 3 + j] = G * m[0] * eta[i] / eta[i - 1]
            * (
                jacobi_x[i * 3 + j] / temp_jacobi_norm_cube
                - temp_vec[j] / temp_vec_norm_cube
            );
        }

        for (int j = 1; j < i; j++)
        {
            // Calculate x_ji
            temp_vec[0] = x[i * 3 + 0] - x[j * 3 + 0];
            temp_vec[1] = x[i * 3 + 1] - x[j * 3 + 1];
            temp_vec[2] = x[i * 3 + 2] - x[j * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

            aux[0] += G * m[j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[i * 3 + 0] -= aux[0] * eta[i] / eta[i - 1];
        a[i * 3 + 1] -= aux[1] * eta[i] / eta[i - 1];
        a[i * 3 + 2] -= aux[2] * eta[i] / eta[i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = i + 1; j < objects_count; j++)
        {
            // Calculate x_ij
            temp_vec[0] = x[j * 3 + 0] - x[i * 3 + 0];
            temp_vec[1] = x[j * 3 + 1] - x[i * 3 + 1];
            temp_vec[2] = x[j * 3 + 2] - x[i * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

            aux[0] += G * m[j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[i * 3 + 0] += aux[0];
        a[i * 3 + 1] += aux[1];
        a[i * 3 + 2] += aux[2];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = 0; j < i; j++)
        {
            for (int k = i + 1; k < objects_count; k++)
            {
                // Calculate x_jk
                temp_vec[0] = x[k * 3 + 0] - x[j * 3 + 0];
                temp_vec[1] = x[k * 3 + 1] - x[j * 3 + 1];
                temp_vec[2] = x[k * 3 + 2] - x[j * 3 + 2];

                temp_vec_norm = vec_norm_3d(temp_vec);
                temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

                aux[0] += G * m[j] * m[k] * temp_vec[0] / temp_vec_norm_cube;
                aux[1] += G * m[j] * m[k] * temp_vec[1] / temp_vec_norm_cube;
                aux[2] += G * m[j] * m[k] * temp_vec[2] / temp_vec_norm_cube;
            }
        }
        a[i * 3 + 0] -= aux[0] / eta[i - 1];
        a[i * 3 + 1] -= aux[1] / eta[i - 1];
        a[i * 3 + 2] -= aux[2] / eta[i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;
    }

    return make_success_error_status();
}

IN_FILE ErrorStatus whfast_acceleration_massless(
    double *__restrict a,
    const System *__restrict system,
    const double *__restrict jacobi_x,
    const double *__restrict eta,
    const AccelerationParam *__restrict acceleration_param
)
{
    const int objects_count = system->objects_count;
    const double *__restrict x = system->x;
    const double *__restrict m = system->m;
    const double G = system->G;

    const double softening_length = acceleration_param->softening_length;

    double aux[3];
    double temp_vec[3];
    double temp_vec_norm;
    double temp_vec_norm_cube;
    double temp_jacobi_norm;
    double temp_jacobi_norm_cube;
    double softening_length_cube = softening_length * softening_length * softening_length;

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

    if (!massive_indices || !massless_indices)
    {
        free(massive_indices);
        free(massless_indices);
        return WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for indices of massive and massless objects"
        );
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

    /* Acceleration calculation for massive objects */
    for (int i = 1; i < massive_objects_count; i++)
    {
        int idx_i = massive_indices[i];

        // Calculate x_0i
        temp_vec[0] = x[idx_i * 3 + 0] - x[0];
        temp_vec[1] = x[idx_i * 3 + 1] - x[1];
        temp_vec[2] = x[idx_i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm_3d(temp_vec);
        temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;
        temp_jacobi_norm = vec_norm_3d(&jacobi_x[idx_i * 3]);
        temp_jacobi_norm_cube = (temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm) + softening_length_cube;
        for (int j = 0; j < 3; j++)
        {
            a[idx_i * 3 + j] = G * m[0] * eta[idx_i] / eta[idx_i - 1]
            * (
                jacobi_x[idx_i * 3 + j] / temp_jacobi_norm_cube
                - temp_vec[j] / temp_vec_norm_cube
            );
        }

        for (int j = 1; j < i; j++)
        {
            int idx_j = massive_indices[j];

            // Calculate x_ji
            temp_vec[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            temp_vec[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            temp_vec[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

            aux[0] += G * m[idx_j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[idx_j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[idx_j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[idx_i * 3 + 0] -= aux[0] * eta[idx_i] / eta[idx_i - 1];
        a[idx_i * 3 + 1] -= aux[1] * eta[idx_i] / eta[idx_i - 1];
        a[idx_i * 3 + 2] -= aux[2] * eta[idx_i] / eta[idx_i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = i + 1; j < massive_objects_count; j++)
        {
            int idx_j = massive_indices[j];

            // Calculate x_ij
            temp_vec[0] = x[idx_j * 3 + 0] - x[idx_i * 3 + 0];
            temp_vec[1] = x[idx_j * 3 + 1] - x[idx_i * 3 + 1];
            temp_vec[2] = x[idx_j * 3 + 2] - x[idx_i * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

            aux[0] += G * m[idx_j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[idx_j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[idx_j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[idx_i * 3 + 0] += aux[0];
        a[idx_i * 3 + 1] += aux[1];
        a[idx_i * 3 + 2] += aux[2];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = 0; j < i; j++)
        {
            int idx_j = massive_indices[j];

            for (int k = i + 1; k < massive_objects_count; k++)
            {
                int idx_k = massive_indices[k];

                // Calculate x_jk
                temp_vec[0] = x[idx_k * 3 + 0] - x[idx_j * 3 + 0];
                temp_vec[1] = x[idx_k * 3 + 1] - x[idx_j * 3 + 1];
                temp_vec[2] = x[idx_k * 3 + 2] - x[idx_j * 3 + 2];

                temp_vec_norm = vec_norm_3d(temp_vec);
                temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

                aux[0] += G * m[idx_j] * m[idx_k] * temp_vec[0] / temp_vec_norm_cube;
                aux[1] += G * m[idx_j] * m[idx_k] * temp_vec[1] / temp_vec_norm_cube;
                aux[2] += G * m[idx_j] * m[idx_k] * temp_vec[2] / temp_vec_norm_cube;
            }
        }
        a[idx_i * 3 + 0] -= aux[0] / eta[idx_i - 1];
        a[idx_i * 3 + 1] -= aux[1] / eta[idx_i - 1];
        a[idx_i * 3 + 2] -= aux[2] / eta[idx_i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;
    }

    /* Acceleration calculation for massless objects */
    for (int i = 0; i < massless_objects_count; i++)
    {
        int idx_i = massless_indices[i];
        if (idx_i == 0)
        {
            continue;
        }
        
        // Calculate x_0i
        temp_vec[0] = x[idx_i * 3 + 0] - x[0];
        temp_vec[1] = x[idx_i * 3 + 1] - x[1];
        temp_vec[2] = x[idx_i * 3 + 2] - x[2];

        temp_vec_norm = vec_norm_3d(temp_vec);
        temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;
        temp_jacobi_norm = vec_norm_3d(&jacobi_x[idx_i * 3]);
        temp_jacobi_norm_cube = (temp_jacobi_norm * temp_jacobi_norm * temp_jacobi_norm) + softening_length_cube;
        for (int j = 0; j < 3; j++)
        {
            a[idx_i * 3 + j] = G * m[0]
            * (
                jacobi_x[idx_i * 3 + j] / temp_jacobi_norm_cube
                - temp_vec[j] / temp_vec_norm_cube
            );
        }

        for (int j = 1; j < massive_objects_count; j++)
        {
            int idx_j = massive_indices[j];
            if (idx_j >= idx_i)
            {
                break;
            }

            // Calculate x_ji
            temp_vec[0] = x[idx_i * 3 + 0] - x[idx_j * 3 + 0];
            temp_vec[1] = x[idx_i * 3 + 1] - x[idx_j * 3 + 1];
            temp_vec[2] = x[idx_i * 3 + 2] - x[idx_j * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

            aux[0] += G * m[idx_j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[idx_j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[idx_j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[idx_i * 3 + 0] -= aux[0];
        a[idx_i * 3 + 1] -= aux[1];
        a[idx_i * 3 + 2] -= aux[2];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = 1; j < massive_objects_count; j++)
        {
            int idx_j = massive_indices[j];
            if (idx_j <= idx_i)
            {
                continue;
            }

            // Calculate x_ij
            temp_vec[0] = x[idx_j * 3 + 0] - x[idx_i * 3 + 0];
            temp_vec[1] = x[idx_j * 3 + 1] - x[idx_i * 3 + 1];
            temp_vec[2] = x[idx_j * 3 + 2] - x[idx_i * 3 + 2];

            temp_vec_norm = vec_norm_3d(temp_vec);
            temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

            aux[0] += G * m[idx_j] * temp_vec[0] / temp_vec_norm_cube;
            aux[1] += G * m[idx_j] * temp_vec[1] / temp_vec_norm_cube;
            aux[2] += G * m[idx_j] * temp_vec[2] / temp_vec_norm_cube;
        }
        a[idx_i * 3 + 0] += aux[0];
        a[idx_i * 3 + 1] += aux[1];
        a[idx_i * 3 + 2] += aux[2];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;

        for (int j = 0; j < massive_objects_count; j++)
        {
            int idx_j = massive_indices[j];
            if (idx_j >= idx_i)
            {
                break;
            }

            for (int k = j + 1; k < massive_objects_count; k++)
            {
                int idx_k = massive_indices[k];
                if (idx_k <= idx_i)
                {
                    continue;
                }

                // Calculate x_jk
                temp_vec[0] = x[idx_k * 3 + 0] - x[idx_j * 3 + 0];
                temp_vec[1] = x[idx_k * 3 + 1] - x[idx_j * 3 + 1];
                temp_vec[2] = x[idx_k * 3 + 2] - x[idx_j * 3 + 2];

                temp_vec_norm = vec_norm_3d(temp_vec);
                temp_vec_norm_cube = (temp_vec_norm * temp_vec_norm * temp_vec_norm) + softening_length_cube;

                aux[0] += G * m[idx_j] * m[idx_k] * temp_vec[0] / temp_vec_norm_cube;
                aux[1] += G * m[idx_j] * m[idx_k] * temp_vec[1] / temp_vec_norm_cube;
                aux[2] += G * m[idx_j] * m[idx_k] * temp_vec[2] / temp_vec_norm_cube;
            }
        }
        a[idx_i * 3 + 0] -= aux[0] / eta[idx_i - 1];
        a[idx_i * 3 + 1] -= aux[1] / eta[idx_i - 1];
        a[idx_i * 3 + 2] -= aux[2] / eta[idx_i - 1];

        aux[0] = 0.0;
        aux[1] = 0.0;
        aux[2] = 0.0;
    }

    free(massive_indices);
    free(massless_indices);

    return make_success_error_status();
}
