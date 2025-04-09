#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "acceleration.h"
#include "error.h"
#include "system.h"
#include "math_functions.h"

IN_FILE void cloud_in_cell(
    double *__restrict delta,
    const double mean_bkg_density,
    const double *__restrict x,
    const double *__restrict m,
    const int num_particles,
    const int pm_grid_size,
    const double *__restrict box_center,
    const double box_length
)
{
    const int grid_size_2 = pm_grid_size * pm_grid_size;
    const int grid_size_3 = pm_grid_size * pm_grid_size * pm_grid_size; 

    /* Clear the grid */
    for (int i = 0; i < grid_size_3; i++)
    {
        delta[i] = 0.0;
    }

    const double box_width = box_length / 2.0;
    const double cell_length = box_length / pm_grid_size;
    const double inv_cell_length = 1.0 / cell_length;
    
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < num_particles; i++)
    {
        const double x_normalized = (x[i * 3 + 0] - box_center[0] + box_width) * inv_cell_length;
        const double y_normalized = (x[i * 3 + 1] - box_center[1] + box_width) * inv_cell_length;
        const double z_normalized = (x[i * 3 + 2] - box_center[2] + box_width) * inv_cell_length;

        const int n_x = (int) x_normalized;
        const int n_y = (int) y_normalized;
        const int n_z = (int) z_normalized;
        
        const double weight_x = x_normalized - n_x;
        const double weight_y = y_normalized - n_y;
        const double weight_z = z_normalized - n_z;

        const double weight_x_m = 1.0 - weight_x;
        const double weight_y_m = 1.0 - weight_y;
        const double weight_z_m = 1.0 - weight_z;

        delta[
            (n_x       % pm_grid_size) * grid_size_2 +
            (n_y       % pm_grid_size) * pm_grid_size +
            (n_z       % pm_grid_size)
        ] += m[i] * weight_x_m * weight_y_m * weight_z_m;
        delta[
            ((n_x + 1) % pm_grid_size) * grid_size_2 +
            (n_y       % pm_grid_size) * pm_grid_size +
            (n_z       % pm_grid_size)
        ] += m[i] * weight_x * weight_y_m * weight_z_m;
        delta[
            (n_x       % pm_grid_size) * grid_size_2 +
            ((n_y + 1) % pm_grid_size) * pm_grid_size +
            (n_z       % pm_grid_size)
        ] += m[i] * weight_x_m * weight_y * weight_z_m;
        delta[
            (n_x       % pm_grid_size) * grid_size_2 +
            (n_y       % pm_grid_size) * pm_grid_size +
            ((n_z + 1) % pm_grid_size)
        ] += m[i] * weight_x_m * weight_y_m * weight_z;
        delta[
            ((n_x + 1) % pm_grid_size) * grid_size_2 +
            ((n_y + 1) % pm_grid_size) * pm_grid_size +
            (n_z       % pm_grid_size)
        ] += m[i] * weight_x * weight_y * weight_z_m;
        delta[
            ((n_x + 1) % pm_grid_size) * grid_size_2 +
            (n_y       % pm_grid_size) * pm_grid_size +
            ((n_z + 1) % pm_grid_size)
        ] += m[i] * weight_x * weight_y_m * weight_z;
        delta[
            (n_x       % pm_grid_size) * grid_size_2 +
            ((n_y + 1) % pm_grid_size) * pm_grid_size +
            ((n_z + 1) % pm_grid_size)
        ] += m[i] * weight_x_m * weight_y * weight_z;
        delta[
            ((n_x + 1) % pm_grid_size) * grid_size_2 +
            ((n_y + 1) % pm_grid_size) * pm_grid_size +
            ((n_z + 1) % pm_grid_size)
        ] += m[i] * weight_x * weight_y * weight_z;
    }

    const double cell_volume = cell_length * cell_length * cell_length;
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < grid_size_3; i++)
    {
        delta[i] = (delta[i] / (cell_volume * mean_bkg_density)) - 1.0;
    }
}

IN_FILE void get_cloud_in_cell_acceleration(
    double *__restrict a,
    const double *__restrict x,
    const double *__restrict acc_grid,
    const int num_particles,
    const int pm_grid_size,
    const double *__restrict box_center,
    const double box_length
)
{
    const double box_width = box_length / 2.0;
    const double cell_length = box_length / pm_grid_size;
    const double inv_cell_length = 1.0 / cell_length;
    const int grid_size_2 = pm_grid_size * pm_grid_size;
    
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < num_particles; i++)
    {
        const double x_normalized = (x[i * 3 + 0] - box_center[0] + box_width) * inv_cell_length;
        const double y_normalized = (x[i * 3 + 1] - box_center[1] + box_width) * inv_cell_length;
        const double z_normalized = (x[i * 3 + 2] - box_center[2] + box_width) * inv_cell_length;

        const int n_x = (int) x_normalized;
        const int n_y = (int) y_normalized;
        const int n_z = (int) z_normalized;

        const double weight_x = x_normalized - n_x;
        const double weight_y = y_normalized - n_y;
        const double weight_z = z_normalized - n_z;

        const double weight_x_m = 1.0 - weight_x;
        const double weight_y_m = 1.0 - weight_y;
        const double weight_z_m = 1.0 - weight_z;

        for (int j = 0; j < 3; j++)
        {
            a[i * 3 + j] = (
                acc_grid[
                    (
                        (n_x       % pm_grid_size) * grid_size_2 +
                        (n_y       % pm_grid_size) * pm_grid_size +
                        (n_z       % pm_grid_size)
                    ) * 3 + j
                ] * weight_x_m * weight_y_m * weight_z_m
                + acc_grid[
                    (
                        ((n_x + 1) % pm_grid_size) * grid_size_2 +
                        (n_y       % pm_grid_size) * pm_grid_size +
                        (n_z       % pm_grid_size)
                    ) * 3 + j
                ] * weight_x * weight_y_m * weight_z_m
                + acc_grid[
                    (
                        (n_x       % pm_grid_size) * grid_size_2 +
                        ((n_y + 1) % pm_grid_size) * pm_grid_size +
                        (n_z       % pm_grid_size)
                    ) * 3 + j
                ] * weight_x_m * weight_y * weight_z_m
                + acc_grid[
                    (
                        (n_x       % pm_grid_size) * grid_size_2 +
                        (n_y       % pm_grid_size) * pm_grid_size +
                        ((n_z + 1) % pm_grid_size)
                    ) * 3 + j
                ] * weight_x_m * weight_y_m * weight_z
                + acc_grid[
                    (
                        ((n_x + 1) % pm_grid_size) * grid_size_2 +
                        ((n_y + 1) % pm_grid_size) * pm_grid_size +
                        (n_z       % pm_grid_size)
                    ) * 3 + j
                ] * weight_x * weight_y * weight_z_m
                + acc_grid[
                    (
                        ((n_x + 1) % pm_grid_size) * grid_size_2 +
                        (n_y       % pm_grid_size) * pm_grid_size +
                        ((n_z + 1) % pm_grid_size)
                    ) * 3 + j
                ] * weight_x * weight_y_m * weight_z
                + acc_grid[
                    (
                        (n_x       % pm_grid_size) * grid_size_2 +
                        ((n_y + 1) % pm_grid_size) * pm_grid_size +
                        ((n_z + 1) % pm_grid_size)
                    ) * 3 + j
                ] * weight_x_m * weight_y * weight_z
                + acc_grid[
                    (
                        ((n_x + 1) % pm_grid_size) * grid_size_2 +
                        ((n_y + 1) % pm_grid_size) * pm_grid_size +
                        ((n_z + 1) % pm_grid_size)
                    ) * 3 + j
                ] * weight_x * weight_y * weight_z
            );
        }
    }
}

IN_FILE void compute_acceleration_with_gradient(
    double *__restrict acc_grid,
    const double *__restrict phi,
    const int pm_grid_size,
    const double box_length
)
{
    const double cell_length = box_length / pm_grid_size;
    const int grid_size_2 = pm_grid_size * pm_grid_size;
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < pm_grid_size; i++)
    {
        const int i_m_2 = (i - 2 + pm_grid_size) % pm_grid_size;
        const int i_m = (i - 1 + pm_grid_size) % pm_grid_size;
        const int i_p = (i + 1) % pm_grid_size;
        const int i_p_2 = (i + 2) % pm_grid_size;
        for (int j = 0; j < pm_grid_size; j++)
        {
            const int j_m_2 = (j - 2 + pm_grid_size) % pm_grid_size;
            const int j_m = (j - 1 + pm_grid_size) % pm_grid_size;
            const int j_p = (j + 1) % pm_grid_size;
            const int j_p_2 = (j + 2) % pm_grid_size;
            for (int k = 0; k < pm_grid_size; k++)
            {
                const int k_m_2 = (k - 2 + pm_grid_size) % pm_grid_size;
                const int k_m = (k - 1 + pm_grid_size) % pm_grid_size;
                const int k_p = (k + 1) % pm_grid_size;
                const int k_p_2 = (k + 2) % pm_grid_size;

                const int index     = i * grid_size_2 + j * pm_grid_size + k;

                const int index_x_m_2 = i_m_2 * grid_size_2 + j * pm_grid_size + k;
                const int index_x_m   = i_m * grid_size_2 + j * pm_grid_size + k;
                const int index_x_p   = i_p * grid_size_2 + j * pm_grid_size + k;
                const int index_x_p_2 = i_p_2 * grid_size_2 + j * pm_grid_size + k;

                const int index_y_m_2 = i * grid_size_2 + j_m_2 * pm_grid_size + k;
                const int index_y_m   = i * grid_size_2 + j_m * pm_grid_size + k;
                const int index_y_p   = i * grid_size_2 + j_p * pm_grid_size + k;
                const int index_y_p_2 = i * grid_size_2 + j_p_2 * pm_grid_size + k;

                const int index_z_m_2 = i * grid_size_2 + j * pm_grid_size + k_m_2;
                const int index_z_m   = i * grid_size_2 + j * pm_grid_size + k_m;
                const int index_z_p   = i * grid_size_2 + j * pm_grid_size + k_p;
                const int index_z_p_2 = i * grid_size_2 + j * pm_grid_size + k_p_2;

                acc_grid[index * 3 + 0] = -(
                    phi[index_x_m_2] - 8.0 * phi[index_x_m] + 8.0 * phi[index_x_p] - phi[index_x_p_2]
                ) / (12.0 * cell_length);
                acc_grid[index * 3 + 1] = -(
                    phi[index_y_m_2] - 8.0 * phi[index_y_m] + 8.0 * phi[index_y_p] - phi[index_y_p_2]
                ) / (12.0 * cell_length);
                acc_grid[index * 3 + 2] = -(
                    phi[index_z_m_2] - 8.0 * phi[index_z_m] + 8.0 * phi[index_z_p] - phi[index_z_p_2]
                ) / (12.0 * cell_length);
            }
        }
    }
}

WIN32DLL_API ErrorStatus acceleration_PM(
    double *__restrict a,
    const CosmologicalSystem *__restrict system,
    const AccelerationParam *__restrict acceleration_param,
    const double mean_bkg_density,
    const int pm_grid_size,
    const double scale_factor
)
{
    ErrorStatus error_status;
    (void) acceleration_param;

    /* Declare variables */
    const int num_particles = system->num_particles;
    const double *__restrict x = system->x;
    const double *__restrict m = system->m;
    const double G = system->G;

    const double *__restrict box_center = system->box_center;
    const double box_width = system->box_width;
    const double box_length = box_width * 2.0;

    const int grid_size_2 = pm_grid_size * pm_grid_size;
    const int grid_size_3 = pm_grid_size * pm_grid_size * pm_grid_size;

    double *__restrict acc_grid = malloc(grid_size_3 * 3 * sizeof(double));
    double *__restrict delta = fftw_malloc(grid_size_3 * sizeof(double));
    fftw_complex *delta_fourier = fftw_malloc(grid_size_2 * (pm_grid_size / 2 + 1) * sizeof(fftw_complex));
    double *__restrict phi = fftw_malloc(grid_size_3 * sizeof(double));
    if (!acc_grid || !delta || !delta_fourier || !phi)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for particle mesh acceleration");
        goto err_malloc;
    }

    /* Deposit mass onto the grid */
    cloud_in_cell(delta, mean_bkg_density, x, m, num_particles, pm_grid_size, box_center, box_length);

    /* Compute the density perturbation in Fourier space */
    fftw_plan plan_forward = fftw_plan_dft_r2c_3d(pm_grid_size, pm_grid_size, pm_grid_size, delta, delta_fourier, FFTW_ESTIMATE);
    if (!plan_forward)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to create FFTW plan for forward transform");
        goto err_fftw_plan_forward;
    }
    fftw_execute(plan_forward);

    const int half_grid_size = pm_grid_size / 2 + 1;
    const double two_pi_over_box_length_squared = 4.0 * M_PI * M_PI / (box_length * box_length);
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < pm_grid_size; i++)
    {
        const int k_x = (i < half_grid_size) ? i : i - pm_grid_size;
        for (int j = 0; j < pm_grid_size; j++)
        {
            const int k_y = (j < half_grid_size) ? j : j - pm_grid_size;
            for (int k = 0; k < half_grid_size; k++)
            {
                const int k_z = k;
                const double k_sq = (k_x * k_x + k_y * k_y + k_z * k_z) * two_pi_over_box_length_squared;

                const double kernel = k_sq == 0.0 ? 0.0 : 1.0 / k_sq;

                const int index = i * (pm_grid_size * half_grid_size) + j * half_grid_size + k;
                delta_fourier[index][0] *= kernel;
                delta_fourier[index][1] *= kernel;
            }
        }
    }

    /* Compute the potential in real space */
    fftw_plan plan_backward = fftw_plan_dft_c2r_3d(pm_grid_size, pm_grid_size, pm_grid_size, delta_fourier, phi, FFTW_ESTIMATE);
    if (!plan_backward)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to create FFTW plan for backward transform");
        goto err_fftw_plan_backward;
    }
    fftw_execute(plan_backward);

    for (int i = 0; i < grid_size_3; i++)
    {
        phi[i] *= G / (scale_factor * (double) grid_size_3);
    }

    /* Compute the force by taking the gradient of the potential */
    compute_acceleration_with_gradient(
        acc_grid,
        phi, pm_grid_size, box_length
    );

    /* Add the acceleration to the particles */
    get_cloud_in_cell_acceleration(
        a, x, acc_grid,
        num_particles, pm_grid_size, box_center, box_length
    );

    /* Free memory */
    free(acc_grid);
    fftw_free(delta);
    fftw_free(delta_fourier);
    fftw_free(phi);
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);

    return make_success_error_status();

err_fftw_plan_backward:
    fftw_destroy_plan(plan_backward);
err_fftw_plan_forward:
    fftw_destroy_plan(plan_forward);
err_malloc:
    free(acc_grid);
    fftw_free(delta);
    fftw_free(delta_fourier);
    fftw_free(phi);

    return error_status;
}
