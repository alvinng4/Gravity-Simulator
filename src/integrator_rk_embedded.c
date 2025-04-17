/**
 * \file integrator_rk_embedded.c
 * \brief Function definitions for Embedded Runge-Kutta integrators
 * 
 * \author Ching-Yin Ng
 * 
 * \ref J. Roa, et al. Moving Planets Around: An Introduction to
 *   N-Body Simulations Applied to Exoplanetary Systems*, MIT
 *   Press, 2020
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration.h"
#include "common.h"
#include "error.h"
#include "integrator.h"
#include "output.h"
#include "progress_bar.h"
#include "settings.h"

/**
 * \brief Get the order of the Embedded RK integrator
 * 
 * \param[out] order Pointer to the order of the integrator
 * \param[in] method Integrator method
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_VALUE_ERROR if the given method is invalid
 */
IN_FILE ErrorStatus get_rk_embedded_order(
    int *restrict order,
    const int method
)
{
    *order = 0;
    switch (method)
    {
        case INTEGRATOR_RKF45:
            *order = 45;
            break;
        case INTEGRATOR_DOPRI:
            *order = 54;
            break;
        case INTEGRATOR_DVERK:
            *order = 65;
            break;
        case INTEGRATOR_RKF78:
            *order = 78;
            break;
        default:
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Invalid integrator method");
    }

    return make_success_error_status();
}

/**
 * \brief Butcher tableaus for Embedded RK integrator
 * 
 * \param[in] order Order of the integrator, must be one of 45 / 54 / 78 / 65
 * \param[out] power Power of the integrator
 * \param[out] power_test Power for error calculation
 * \param[out] coeff Pointer to array of coefficients for the integrator
 * \param[out] len_weights Length of the weights array
 * \param[out] weights Pointer to the array of weights for RK integrator
 * \param[out] weights_test Pointer to the array of weights for error calculation
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_VALUE_ERROR if the given order is invalid
 * \exception GRAV_MEMORY_ERROR if failed to allocate memory for coeff, weights and weights_test
 */
IN_FILE ErrorStatus rk_embedded_butcher_tableaus(
    const int order,
    int *restrict power,
    int *restrict power_test,
    double **coeff,
    int *restrict len_weights,
    double **weights,
    double **weights_test
)
{
    /*  
    *   Select integrator
    *   45) Runge-Kutta-Fehleberg 4(5)
    *   54) Dormand-Prince 5(4)
    *   78) Runge-Kutta-Fehlberg 7(8)
    *   65) Verner's method 6(5), DVERK
    */

    *power = -1;
    *power_test = -1;
    *len_weights = -1;
    *coeff = NULL;
    *weights = NULL;
    *weights_test = NULL;

    switch (order)
    {
        // RUNGE-KUTTA-FEHLBERG 4(5)
        case 45:
            // Order
            *power = 4;
            *power_test = 5;
            // nodes = np.array([1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5])
            *coeff = malloc(25 * sizeof(double));
            *len_weights = 6;
            *weights = malloc(*len_weights * sizeof(double));
            *weights_test = malloc(6 * sizeof(double));
            if (!*coeff || !*weights || !*weights_test)
            {
                free(*coeff);
                free(*weights);
                free(*weights_test);

                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for Butcher tableau");
            }
            memcpy(
                *coeff,
                (double [25]) {
                    1.0L / 4.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    3.0L / 32.0L, 9.0L / 32.0L, 0.0L, 0.0L, 0.0L,
                    1932.0L / 2197.0L, -7200.0L / 2197.0L, 7296.0L / 2197.0L, 0.0L, 0.0L,
                    439.0L / 216.0L, -8.0L, 3680.0L / 513.0L, -845.0L / 4104.0L, 0.0L,
                    -8.0L / 27.0L, 2.0L, -3544.0L / 2565.0L, 1859.0L / 4104.0L, -11.0L / 40.0L
                },
                25 * sizeof(double)
            );
            memcpy(
                *weights,
                (double [6]) {
                    25.0L / 216.0L, 0.0L, 1408.0L / 2565.0L, 2197.0L / 4104.0L, -0.2L, 0.0L
                },
                6 * sizeof(double)
            );
            memcpy(
                *weights_test,
                (double [6]) {
                    16.0L / 135.0L, 0.0L, 6656.0L / 12825.0L, 28561.0L / 56430.0L, -9.0L / 50.0L, 2.0L / 55.0L
                },
                6 * sizeof(double)
            );

            break;

        // DORMAND-PRINCE 5(4)
        case 54:
            // order
            *power = 5;
            *power_test = 4;
            // nodes = np.array([1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0])
            *coeff = malloc(36 * sizeof(double));
            *len_weights = 7;
            *weights = malloc(*len_weights * sizeof(double));
            *weights_test = malloc(7 * sizeof(double));

            if (!*coeff || !*weights || !*weights_test)
            {
                free(*coeff);
                free(*weights);
                free(*weights_test);

                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for Butcher tableau");
            }

            memcpy(
                *coeff,
                (double [36]) {
                    1.0L / 5.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    3.0L / 40.0L, 9.0L / 40.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    44.0L / 45.0L, -56.0L / 15.0L, 32.0L / 9.0L, 0.0L, 0.0L, 0.0L,
                    19372.0L / 6561.0L, -25360.0L / 2187.0L, 64448.0L / 6561.0L, -212.0L / 729.0L, 0.0L, 0.0L,
                    9017.0L / 3168.0L, -355.0L / 33.0L, 46732.0L / 5247.0L, 49.0L / 176.0L, -5103.0L / 18656.0L, 0.0L,
                    35.0L / 384.0L, 0.0L, 500.0L / 1113.0L, 125.0L / 192.0L, -2187.0L / 6784.0L, 11.0L / 84.0L
                },
                36 * sizeof(double)
            );
            memcpy(
                *weights,
                (double [7]) {
                    35.0L / 384.0L, 0.0L, 500.0L / 1113.0L, 125.0L / 192.0L, -2187.0L / 6784.0L, 11.0L / 84.0L, 0.0L
                },
                7 * sizeof(double)
            );
            memcpy(
                *weights_test,
                (double [7]) {
                    5179.0L / 57600.0L, 0.0L, 7571.0L / 16695.0L, 393.0L / 640.0L, -92097.0L / 339200.0L, 187.0L / 2100.0L, 1.0L / 40.0L
                },
                7 * sizeof(double)
            );

            break;

        // RUNGE-KUTTA-FEHLBERG 7(8)
        case 78:
            // Order
            *power = 7;
            *power_test = 8;
            // nodes = np.array(
            //     [
            //         2.0 / 27.0,
            //         1.0 / 9.0,
            //         1.0 / 6.0,
            //         5.0 / 12.0,
            //         1.0 / 2.0,
            //         5.0 / 6.0,
            //         1.0 / 6.0,
            //         2.0 / 3.0,
            //         1.0 / 3.0,
            //         1.0,
            //         0.0,
            //         1.0,
            //     ]
            // )
            
            *coeff = malloc(144 * sizeof(double));
            *len_weights = 13;
            *weights = malloc(*len_weights * sizeof(double));
            *weights_test = malloc(13 * sizeof(double));
            if (!*coeff || !*weights || !*weights_test)
            {
                free(*coeff);
                free(*weights);
                free(*weights_test);

                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for Butcher tableau");
            }

            memcpy(
                *coeff,
                (double [144]) {
                    2.0L / 27.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    1.0L / 36.0L, 1.0L / 12.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    1.0L / 24.0L, 0.0L, 1.0L / 8.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    5.0L / 12.0L, 0.0L, -25.0L / 16.0L, 25.0L / 16.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    1.0L / 20.0L, 0.0L, 0.0L, 1.0L / 4.0L, 1.0L / 5.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    -25.0L / 108.0L, 0.0L, 0.0L, 125.0L / 108.0L, -65.0L / 27.0L, 125.0L / 54.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    31.0L / 300.0L, 0.0L, 0.0L, 0.0L, 61.0L / 225.0L, -2.0L / 9.0L, 13.0L / 900.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    2.0L, 0.0L, 0.0L, -53.0L / 6.0L, 704.0L / 45.0L, -107.0L / 9.0L, 67.0L / 90.0L, 3.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    -91.0L / 108.0L, 0.0L, 0.0L, 23.0L / 108.0L, -976.0L / 135.0L, 311.0L / 54.0L, -19.0L / 60.0L, 17.0L / 6.0L, -1.0L / 12.0L, 0.0L, 0.0L, 0.0L,
                    2383.0L / 4100.0L, 0.0L, 0.0L, -341.0L / 164.0L, 4496.0L / 1025.0L, -301.0L / 82.0L, 2133.0L / 4100.0L, 45.0L / 82.0L, 45.0L / 164.0L, 18.0L / 41.0L, 0.0L, 0.0L,
                    3.0L / 205.0L, 0.0L, 0.0L, 0.0L, 0.0L, -6.0L / 41.0L, -3.0L / 205.0L, -3.0L / 41.0L, 3.0L / 41.0L, 6.0L / 41.0L, 0.0L, 0.0L,
                    -1777.0L / 4100.0L, 0.0L, 0.0L, -341.0L / 164.0L, 4496.0L / 1025.0L, -289.0L / 82.0L, 2193.0L / 4100.0L, 51.0L / 82.0L, 33.0L / 164.0L, 19.0L / 41.0L, 0.0L, 1.0L
                },
                144 * sizeof(double)
            );
            memcpy(
                *weights,
                (double [13]) {
                    41.0L / 840.0L, 0.0L, 0.0L, 0.0L, 0.0L, 34.0L / 105.0L, 9.0L / 35.0L, 9.0L / 35.0L, 9.0L / 280.0L, 9.0L / 280.0L, 41.0L / 840.0L, 0.0L, 0.0L
                },
                13 * sizeof(double)
            );
            memcpy(
                *weights_test,
                (double [13]) {
                    0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 34.0L / 105.0L, 9.0L / 35.0L, 9.0L / 35.0L, 9.0L / 280.0L, 9.0L / 280.0L, 0.0L, 41.0L / 840.0L, 41.0L / 840.0L
                },
                13 * sizeof(double)
            );

            break;

        // VERNER 6(5) DVERK
        case 65:
            // Order
            *power = 6;
            *power_test = 7;
            /* nodes = np.array(
            *     [1.0 / 6.0, 4.0 / 15.0, 2.0 / 3.0, 5.0 / 6.0, 1.0, 1.0 / 15.0, 1.0]
            * )
            */
            *coeff = malloc(49 * sizeof(double));
            *len_weights = 8;
            *weights = malloc(*len_weights * sizeof(double));
            *weights_test = malloc(8 * sizeof(double));
            if (!*coeff || !*weights || !*weights_test)
            {
                free(*coeff);
                free(*weights);
                free(*weights_test);

                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for Butcher tableau");
            }
            memcpy(
                *coeff,
                (double [49]) {
                    1.0L / 6.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    4.0L / 75.0L, 16.0L / 75.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    5.0L / 6.0L, -8.0L / 3.0L, 5.0L / 2.0L, 0.0L, 0.0L, 0.0L, 0.0L,
                    -165.0L / 64.0L, 55.0L / 6.0L, -425.0L / 64.0L, 85.0L / 96.0L, 0.0L, 0.0L, 0.0L,
                    12.0L / 5.0L, -8.0L, 4015.0L / 612.0L, -11.0L / 36.0L, 88.0L / 255.0L, 0.0L, 0.0L,
                    -8263.0L / 15000.0L, 124.0L / 75.0L, -643.0L / 680.0L, -81.0L / 250.0L, 2484.0L / 10625.0L, 0.0L, 0.0L,
                    3501.0L / 1720.0L, -300.0L / 43.0L, 297275.0L / 52632.0L, -319.0L / 2322.0L, 24068.0L / 84065.0L, 0.0L, 3850.0L / 26703.0L
                },
                49 * sizeof(double)
            );
            memcpy(
                *weights,
                (double [8]) {
                    3.0L / 40.0L, 0.0L, 875.0L / 2244.0L, 23.0L / 72.0L, 264.0L / 1955.0L, 0.0L, 125.0L / 11592.0L, 43.0L / 616.0L
                },
                8 * sizeof(double)
            );
            memcpy(
                *weights_test,
                (double [8]) {
                    13.0L / 160.0L, 0.0L, 2375.0L / 5984.0L, 5.0L / 16.0L, 12.0L / 85.0L, 3.0L / 44.0L, 0.0L, 0.0L
                },
                8 * sizeof(double)
            );

            break;

        default:
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Invalid order for Embedded RK integrator");
    }

    return make_success_error_status();
}

/**
 * \brief Calculate initial time step for Embedded Runge-Kutta integrator
 * 
 * \param[out] initial_dt Pointer to the initial time step
 * \param[in] rel_tolerance Relative tolerance
 * \param[in] abs_tolerance Absolute tolerance
 * \param[in] power Power of the integrator
 * \param[in] system Pointer to the system
 * \param[in] acceleration_param Pointer to the acceleration parameters
 * 
 * \note Modified to return dt * 1e-2 since this function gives initial dt thats too large
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if failed to allocate memory for arrays
 * \exception GRAV_VALUE_ERROR if initial_dt is negative
 */
IN_FILE ErrorStatus rk_embedded_initial_dt(
    double *restrict initial_dt,
    const double rel_tolerance,
    const double abs_tolerance,
    const int power,
    const System *system,
    const AccelerationParam *acceleration_param
)
{
    ErrorStatus error_status;
    *initial_dt = -1.0;

    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;

    /* Allocate memory and declare variables */
    double *restrict tolerance_scale_x = malloc(num_particles * 3 * sizeof(double));
    double *restrict tolerance_scale_v = malloc(num_particles * 3 * sizeof(double));
    double *restrict x_1 = malloc(num_particles * 3 * sizeof(double));
    double *restrict v_1 = malloc(num_particles * 3 * sizeof(double));
    double *restrict a_1 = malloc(num_particles * 3 * sizeof(double));
    double *restrict a = malloc(num_particles * 3 * sizeof(double));
    if (
        !tolerance_scale_x ||
        !tolerance_scale_v ||
        !x_1 ||
        !v_1 ||
        !a_1 ||
        !a
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for arrays");
        goto error_memory;
    }

    double sum_0 = 0.0;
    double sum_1 = 0.0;
    double sum_2 = 0.0;
    double d_0;
    double d_1;
    double d_2;
    double dt_0;
    double dt_1;
    double dt;
    System system_1 = {
        .num_particles = num_particles,
        .x = x_1,
        .v = v_1,
        .m = system->m,
        .G = system->G,
    };

    /* Compute acceleration */
    error_status = WRAP_TRACEBACK(acceleration(
        a,
        system,
        acceleration_param
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error_acc;
    }

    /**
     *  tolerance_scale_x = abs_tol + rel_tol * abs(x)
     *  tolerance_scale_v = abs_tol + rel_tol * abs(v)
     */
    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            tolerance_scale_x[i * 3 + j] = abs_tolerance + rel_tolerance * fabs(x[i * 3 + j]);
            tolerance_scale_v[i * 3 + j] = abs_tolerance + rel_tolerance * fabs(v[i * 3 + j]);
        }
    }

    /**
     *  sum_0 = sum(square(x / tolerance_scale_x)) + sum(square(v / tolerance_scale_v))
     *  sum_1 = sum(square(v / tolerance_scale_x)) + sum(square(a / tolerance_scale_x))
     */
    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sum_0 += (x[i * 3 + j] / tolerance_scale_x[i * 3 + j]) * (x[i * 3 + j] / tolerance_scale_x[i * 3 + j]);
            sum_0 += (v[i * 3 + j] / tolerance_scale_v[i * 3 + j]) * (v[i * 3 + j] / tolerance_scale_v[i * 3 + j]);
            sum_1 += (v[i * 3 + j] / tolerance_scale_x[i * 3 + j]) * (v[i * 3 + j] / tolerance_scale_x[i * 3 + j]);
            sum_1 += (a[i * 3 + j] / tolerance_scale_v[i * 3 + j]) * (a[i * 3 + j] / tolerance_scale_v[i * 3 + j]);
        }
    }

    d_0 = sqrt(sum_0 / (num_particles * 6));
    d_1 = sqrt(sum_1 / (num_particles * 6));

    if (d_0 < 1e-5 || d_1 < 1e-5)
    {
        dt_0 = 1e-4;
    }
    else
    {
        dt_0 = d_0 / d_1;
    }

    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            x_1[i * 3 + j] = x[i * 3 + j] + (dt_0 / 100.0) * v[i * 3 + j];
            v_1[i * 3 + j] = v[i * 3 + j] + (dt_0 / 100.0) * a[i * 3 + j];
        }
    }

    error_status = WRAP_TRACEBACK(acceleration(
        a_1,
        &system_1,
        acceleration_param
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error_acc;
    }

    /* Calculate d_2 to measure how much the derivatives have changed */

    /**
     * sum_2 = sum(square((v_1 - v) / tolerance_scale_x)) + sum(square((a_1 - a) / tolerance_scale_v))
     */
    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sum_2 += ((v_1[i * 3 + j] - v[i * 3 + j]) / tolerance_scale_x[i * 3 + j]) * ((v_1[i * 3 + j] - v[i * 3 + j]) / tolerance_scale_x[i * 3 + j]);
            sum_2 += ((a_1[i * 3 + j] - a[i * 3 + j]) / tolerance_scale_v[i * 3 + j]) * ((a_1[i * 3 + j] - a[i * 3 + j]) / tolerance_scale_v[i * 3 + j]);
        }
    }
    d_2 = sqrt(sum_2 / (num_particles * 6)) / dt_0;

    if (fmax(d_1, d_2) <= 1e-15)
    {
        dt_1 = fmax(1e-6L, dt_0 * 1e-3L);
    }
    {
        dt_1 = pow((0.01L / fmax(d_1, d_2)), (1.0L / (1 + power)));
    }
    dt = fmin(100.0L * dt_0, dt_1);

    /* Modified to multiply by 1e-2 */
    *initial_dt = dt * 1e-2;
    if (*initial_dt <= 0.0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initial dt is negative");
        goto error_dt;
    }

    free(tolerance_scale_x);
    free(tolerance_scale_v);
    free(x_1);
    free(v_1);
    free(a_1);
    free(a);

    return make_success_error_status();

error_dt:
error_acc:
error_memory:
    free(tolerance_scale_x);
    free(tolerance_scale_v);
    free(x_1);
    free(v_1);
    free(a_1);
    free(a);
    return error_status;
}

WIN32DLL_API ErrorStatus rk_embedded(
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

    /* Initialization */
    int order;
    error_status = WRAP_TRACEBACK(get_rk_embedded_order(
        &order,
        integrator_param->integrator
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    int power;
    int power_test;
    double *coeff = NULL;
    int len_weights;
    double *weights = NULL;
    double *weights_test = NULL;

    error_status = WRAP_TRACEBACK(rk_embedded_butcher_tableaus(
        order,
        &power,
        &power_test,
        &coeff,
        &len_weights,
        &weights,
        &weights_test
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    const int stages = len_weights;
    const int min_power = power < power_test ? power : power_test;

    double *restrict error_estimation_delta_weights = malloc(len_weights * sizeof(double));
    if (!error_estimation_delta_weights)
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for error estimation delta weights"
        );
        goto error_error_estimation_delta_weights_memory_alloc;
    }

    for (int stage = 0; stage < stages; stage++)
    {
        error_estimation_delta_weights[stage] = weights[stage] - weights_test[stage];
    }

    /* tolerance */
    const double abs_tolerance = integrator_param->tolerance;
    const double rel_tolerance = integrator_param->tolerance;

    /* Safety factors for step-size control */
    const double safety_fac_max = 6.0;
    const double safety_fac_min = 0.33;
    const double safety_fac = pow(0.38, (1.0 / (1.0 + (double) min_power)));

    /* Declare variables */
    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *restrict output_count_ptr = &(output_param->output_count_);
    const double output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    double *restrict t_ptr = &(simulation_status->t);
    int64 *restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;

    double sum;
    double error;
    double dt_new;

    System temp_system = {
        .num_particles = num_particles,
        .x = NULL,
        .v = NULL,
        .m = system->m,
        .G = system->G,
    };

    /* Allocate memory */
    double *restrict v_1 = malloc(num_particles * 3 * sizeof(double));
    double *restrict x_1 = malloc(num_particles * 3 * sizeof(double));
    double *restrict vk = malloc(stages * num_particles * 3 * sizeof(double));
    double *restrict xk = malloc(stages * num_particles * 3 * sizeof(double));    
    double *restrict temp_v = malloc(num_particles * 3 * sizeof(double));
    double *restrict temp_x = malloc(num_particles * 3 * sizeof(double));
    double *restrict error_estimation_delta_v = malloc(num_particles * 3 * sizeof(double));
    double *restrict error_estimation_delta_x = malloc(num_particles * 3 * sizeof(double));
    double *restrict tolerance_scale_v = malloc(num_particles * 3 * sizeof(double));
    double *restrict tolerance_scale_x = malloc(num_particles * 3 * sizeof(double));

    // Compensated summation
    double *restrict x_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict v_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict temp_x_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict temp_v_err_comp_sum = calloc(num_particles * 3, sizeof(double));

    if (
        !v_1 ||
        !x_1 ||
        !vk ||
        !xk ||
        !temp_v ||
        !temp_x ||
        !error_estimation_delta_v ||
        !error_estimation_delta_x ||
        !tolerance_scale_v ||
        !tolerance_scale_x ||
        !x_err_comp_sum ||
        !v_err_comp_sum ||
        !temp_x_err_comp_sum ||
        !temp_v_err_comp_sum
    ) 
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for arrays"
        );
        goto error_memory;
    }

    /* Get initial dt */
    double dt;
    if (integrator_param->initial_dt > 0.0)
    {
        dt = integrator_param->initial_dt;
    }
    else
    {
        error_status = WRAP_TRACEBACK(rk_embedded_initial_dt(
            &dt,
            rel_tolerance,
            abs_tolerance,
            power,
            system,
            acceleration_param
        ));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto error_initial_dt;
        }

        if (dt > tf)
        {
            dt = tf;
        }
    }

    /* Initial output */
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

    /* Main Loop */
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        error_status = WRAP_TRACEBACK(start_progress_bar(&progress_bar_param, tf));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_start_progress_bar;
        }
    }
    
    *t_ptr = 0.0;
    simulation_status->dt = dt;
    *num_steps_ptr = 0;
    while (*t_ptr < tf)
    {
        /* Compute xk and vk */
        error_status = WRAP_TRACEBACK(acceleration(
            vk,
            system,
            acceleration_param
        ));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto acc_error;
        }

        memcpy(xk, v, num_particles * 3 * sizeof(double));
        for (int stage = 1; stage < stages; stage++)
        {
            // Empty temp_v and temp_x
            for (int i = 0; i < num_particles; i++)
            {
                temp_v[i * 3 + 0] = 0.0;
                temp_v[i * 3 + 1] = 0.0;
                temp_v[i * 3 + 2] = 0.0;
                temp_x[i * 3 + 0] = 0.0;
                temp_x[i * 3 + 1] = 0.0;
                temp_x[i * 3 + 2] = 0.0;
            }

            for (int i = 0; i < stage; i++)
            {
                for (int j = 0; j < num_particles; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        temp_v[j * 3 + k] += coeff[(stage - 1) * (stages - 1) + i] * vk[i * num_particles * 3 + j * 3 + k];
                        temp_x[j * 3 + k] += coeff[(stage - 1) * (stages - 1) + i] * xk[i * num_particles * 3 + j * 3 + k];
                    }
                }
            }

            for (int i = 0; i < num_particles; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    temp_v[i * 3 + j] = v[i * 3 + j] + dt * temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                    temp_x[i * 3 + j] = x[i * 3 + j] + dt * temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                }
            }

            temp_system.x = temp_x;
            temp_system.v = temp_v;
            error_status = WRAP_TRACEBACK(acceleration(
                &vk[stage * num_particles * 3],
                &temp_system,
                acceleration_param
            ));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                goto acc_error;
            }
            memcpy(&xk[stage * num_particles * 3], temp_v, num_particles * 3 * sizeof(double));
        }

        // Empty temp_v, temp_x, error_estimation_delta_v, error_estimation_delta_x
        for (int i = 0; i < num_particles; i++)
        {
            temp_v[i * 3 + 0] = 0.0;
            temp_v[i * 3 + 1] = 0.0;
            temp_v[i * 3 + 2] = 0.0;
            temp_x[i * 3 + 0] = 0.0;
            temp_x[i * 3 + 1] = 0.0;
            temp_x[i * 3 + 2] = 0.0;
            error_estimation_delta_v[i * 3 + 0] = 0.0;
            error_estimation_delta_v[i * 3 + 1] = 0.0;
            error_estimation_delta_v[i * 3 + 2] = 0.0;
            error_estimation_delta_x[i * 3 + 0] = 0.0;
            error_estimation_delta_x[i * 3 + 1] = 0.0;
            error_estimation_delta_x[i * 3 + 2] = 0.0;
        }

        // Calculate x_1, v_1 and also delta x, delta v for error estimation
        for(int stage = 0; stage < stages; stage++)
        {
            for (int i = 0; i < num_particles; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    temp_v[i * 3 + j] += weights[stage] * vk[stage * num_particles * 3 + i * 3 + j];
                    temp_x[i * 3 + j] += weights[stage] * xk[stage * num_particles * 3 + i * 3 + j];

                    error_estimation_delta_v[i * 3 + j] += dt * error_estimation_delta_weights[stage] * vk[stage * num_particles * 3 + i * 3 + j];
                    error_estimation_delta_x[i * 3 + j] += dt * error_estimation_delta_weights[stage] * xk[stage * num_particles * 3 + i * 3 + j];
                }
            }
        }

        memcpy(temp_x_err_comp_sum, x_err_comp_sum, num_particles * 3 * sizeof(double));
        memcpy(temp_v_err_comp_sum, v_err_comp_sum, num_particles * 3 * sizeof(double));
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                temp_v_err_comp_sum[i * 3 + j] += dt * temp_v[i * 3 + j];
                temp_x_err_comp_sum[i * 3 + j] += dt * temp_x[i * 3 + j];

                v_1[i * 3 + j] = v[i * 3 + j] + temp_v_err_comp_sum[i * 3 + j];
                x_1[i * 3 + j] = x[i * 3 + j] + temp_x_err_comp_sum[i * 3 + j];

                temp_v_err_comp_sum[i * 3 + j] += v[i * 3 + j] - v_1[i * 3 + j];
                temp_x_err_comp_sum[i * 3 + j] += x[i * 3 + j] - x_1[i * 3 + j];
            }
        }

        /* Error calculation */
        for (int i = 0; i < num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                tolerance_scale_v[i * 3 + j] = abs_tolerance + fmax(fabs(v[i * 3 + j]), fabs(v_1[i * 3 + j])) * rel_tolerance;
                tolerance_scale_x[i * 3 + j] = abs_tolerance + fmax(fabs(x[i * 3 + j]), fabs(x_1[i * 3 + j])) * rel_tolerance;
            }
        }

        // Sum up all the elements of x/tol and v/tol, 
        // square and divide by the total number of elements
        sum = 0.0;
        for (int i = 0; i < num_particles; i++)
        {
            double temp;
            for (int j = 0; j < 3; j++)
            {
                temp = error_estimation_delta_v[i * 3 + j] / tolerance_scale_v[i * 3 + j];
                sum += temp * temp;
                temp = error_estimation_delta_x[i * 3 + j] / tolerance_scale_x[i * 3 + j];
                sum += temp * temp;
            }
        }
        error = sqrt(sum / (num_particles * 3 * 2));

        /* Advance step */
        if (error <= 1 || dt <= tf * 1e-12)
        {
            (*num_steps_ptr)++;
            *t_ptr += dt;

            memcpy(x, x_1, num_particles * 3 * sizeof(double));
            memcpy(v, v_1, num_particles * 3 * sizeof(double));

            memcpy(x_err_comp_sum, temp_x_err_comp_sum, num_particles * 3 * sizeof(double));
            memcpy(v_err_comp_sum, temp_v_err_comp_sum, num_particles * 3 * sizeof(double));

            /* Output */
            if (is_output && *t_ptr >= next_output_time)
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
                    goto err_output;
                }

                next_output_time = (*output_count_ptr) * output_interval;
            }
        }

        if (enable_progress_bar)
        {
            update_progress_bar(&progress_bar_param, *t_ptr, false);
        }

        /* Calculate dt for next step */
        if (error < 1e-10)
        {
            error = 1e-10;  // Prevent error from being too small
        }
        dt_new = dt * safety_fac / pow(error, (1.0 / (1.0 + (double) min_power)));
        if (dt_new > safety_fac_max * dt) 
        {
            dt *= safety_fac_max;
        }
        else if (dt_new < safety_fac_min * dt)
        {
            dt *= safety_fac_min;
        }
        else
        {
            dt = dt_new;
        }

        if (dt_new / tf < 1e-12)
        {
            dt = tf * 1e-12;
        }

        // Correct overshooting
        if (((*t_ptr) < tf) && (((*t_ptr) + dt) > tf))
        {
            dt = tf - (*t_ptr);
        }
        simulation_status->dt = dt;

        /* Check exit */
        if (*(settings->is_exit_ptr))
        {
            break;
        }
    }

    if (enable_progress_bar)
    {
        update_progress_bar(&progress_bar_param, *t_ptr, true);
    }

    /* free memory */
    free(v_1);
    free(x_1);
    free(vk);
    free(xk);
    free(temp_v);
    free(temp_x);
    free(error_estimation_delta_v);
    free(error_estimation_delta_x);
    free(tolerance_scale_v);
    free(tolerance_scale_x);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    free(temp_x_err_comp_sum);
    free(temp_v_err_comp_sum);
    free(error_estimation_delta_weights);
    free(coeff);
    free(weights);
    free(weights_test); 

    return make_success_error_status();

err_output:
acc_error:
err_start_progress_bar:
err_initial_output:
error_initial_dt:
error_memory:
    free(v_1);
    free(x_1);
    free(vk);
    free(xk);
    free(temp_v);
    free(temp_x);
    free(error_estimation_delta_v);
    free(error_estimation_delta_x);
    free(tolerance_scale_v);
    free(tolerance_scale_x);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    free(temp_x_err_comp_sum);
    free(temp_v_err_comp_sum);
error_error_estimation_delta_weights_memory_alloc:
    free(error_estimation_delta_weights);
    free(coeff);
    free(weights);
    free(weights_test); 

    return error_status;
}
