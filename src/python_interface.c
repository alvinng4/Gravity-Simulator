#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "grav_sim.h"
#include "math_functions.h"

WIN32DLL_API void free_memory_int32(int32 *__restrict ptr)
{
    free(ptr);
}

WIN32DLL_API void free_memory_double(double *__restrict ptr)
{
    free(ptr);
}

WIN32DLL_API int32 load_built_in_system_python(
    const char *system_name,
    int *num_particles_ptr,
    int32 **particle_ids_ptr,
    double **x_ptr,
    double **v_ptr,
    double **m_ptr,
    double *G_ptr
)
{
    ErrorStatus error_status;

    System system = get_new_system();
    error_status = initialize_built_in_system(
        &system,
        system_name,
        false
    );
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto err;
    }

    *num_particles_ptr = system.num_particles;
    *particle_ids_ptr = system.particle_ids;
    *x_ptr = system.x;
    *v_ptr = system.v;
    *m_ptr = system.m;
    *G_ptr = system.G;

    return GRAV_SUCCESS;

err:
    print_and_free_traceback(&error_status);
    return error_status.return_code;
}

WIN32DLL_API void keplerian_to_cartesian_python(
    double *__restrict x,
    double *__restrict y,
    double *__restrict z,
    double *__restrict v_x,
    double *__restrict v_y,
    double *__restrict v_z,
    const double semi_major_axis,
    const double eccentricity,
    const double inclination,
    const double argument_of_periapsis,
    const double longitude_of_ascending_node,
    const double true_anomaly,
    const double total_mass,
    const double G
)
{
    double x_arr[3];
    double v_arr[3];
    keplerian_to_cartesian(
        x_arr,
        v_arr,
        semi_major_axis,
        eccentricity,
        inclination,
        argument_of_periapsis,
        longitude_of_ascending_node,
        true_anomaly,
        total_mass,
        G
    );
    *x = x_arr[0];
    *y = x_arr[1];
    *z = x_arr[2];
    *v_x = v_arr[0];
    *v_y = v_arr[1];
    *v_z = v_arr[2];
}

WIN32DLL_API int launch_simulation_python(
    int32 *__restrict num_particles,
    int32 *particle_ids,
    double *x,
    double *v,
    double *m,
    int32 **new_particle_ids_ptr,
    double **new_x_ptr,
    double **new_v_ptr,
    double **new_m_ptr,
    const double G,
    const int32 integrator,
    const double dt,
    const double tolerance,
    const double initial_dt,
    const bool whfast_remove_invalid_particles,
    const int32 acceleration_method,
    const double opening_angle,
    const double softening_length,
    const int32 max_num_particles_per_leaf,
    const int32 output_method,
    char *output_dir,
    const bool output_initial,
    const double output_interval,
    const int32 coordinate_output_dtype,
    const int32 velocity_output_dtype,
    const int32 mass_output_dtype,
    const int32 verbose,
    const bool enable_progress_bar,
    bool *is_exit_ptr,
    const double tf
)
{
    /* System */
    System system = get_new_system();
    system.num_particles = *num_particles;
    system.particle_ids = particle_ids;
    system.x = x;
    system.v = v;
    system.m = m;
    system.G = G;

    /* Acceleration parameters */
    AccelerationParam acceleration_param = get_new_acceleration_param();
    acceleration_param.method = acceleration_method;
    acceleration_param.opening_angle = opening_angle;
    acceleration_param.softening_length = softening_length;
    acceleration_param.max_num_particles_per_leaf = max_num_particles_per_leaf;

    /* Integrator parameters */
    IntegratorParam integrator_param = get_new_integrator_param();
    integrator_param.integrator = integrator;
    integrator_param.dt = dt;
    integrator_param.tolerance = tolerance;
    integrator_param.initial_dt = initial_dt;
    integrator_param.whfast_remove_invalid_particles = whfast_remove_invalid_particles;

    /* Output parameters */
    OutputParam output_param = get_new_output_param();
    output_param.method = output_method;
    output_param.output_dir = output_dir;
    output_param.output_interval = output_interval;
    output_param.output_initial = output_initial;
    output_param.coordinate_output_dtype = coordinate_output_dtype;
    output_param.velocity_output_dtype = velocity_output_dtype;
    output_param.mass_output_dtype = mass_output_dtype;

    /* Simulation status */
    SimulationStatus simulation_status;

    /* Settings */
    Settings settings = get_new_settings();
    settings.verbose = verbose;
    settings.enable_progress_bar = enable_progress_bar;
    settings.is_exit_ptr = is_exit_ptr;

    /* Launch simulation */
    ErrorStatus error_status = WRAP_TRACEBACK(launch_simulation(
        &system,
        &integrator_param,
        &acceleration_param,
        &output_param,
        &simulation_status,
        &settings,
        tf
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        print_and_free_traceback(&error_status);
        return 1;
    }

    // Note: no need to free memory since the memory belongs to python
    *new_particle_ids_ptr = system.particle_ids;
    *new_x_ptr = system.x;
    *new_v_ptr = system.v;
    *new_m_ptr = system.m;

    return 0;
}

WIN32DLL_API void compute_energy_python(
    double *__restrict energy,
    const double G,
    const double *__restrict sol_state,
    const int32 num_snapshots,
    const int32 num_particles
)
{
    const int size_snapshot = num_particles * 7;
    for (int n = 0; n < num_snapshots; n++)
    {
        energy[n] = 0.0;
        for (int i = 0; i < num_particles; i++)
        {
            const double m_i = sol_state[n * size_snapshot + i * 7 + 0];
            const double x_i[3] = {
                sol_state[n * size_snapshot + i * 7 + 1],
                sol_state[n * size_snapshot + i * 7 + 2],
                sol_state[n * size_snapshot + i * 7 + 3]
            };
            const double v_i[3] = {
                sol_state[n * size_snapshot + i * 7 + 4],
                sol_state[n * size_snapshot + i * 7 + 5],
                sol_state[n * size_snapshot + i * 7 + 6]
            };

            // KE
            double v_norm = vec_norm_3d(v_i);
            energy[n] += 0.5 * m_i * v_norm * v_norm;

            // PE
            for (int j = i + 1; j < num_particles; j++)
            {
                const double m_j = sol_state[n * size_snapshot + j * 7 + 0];
                const double x_j[3] = {
                    sol_state[n * size_snapshot + j * 7 + 1],
                    sol_state[n * size_snapshot + j * 7 + 2],
                    sol_state[n * size_snapshot + j * 7 + 3]
                };
                const double r_ij[3] = {
                    x_i[0] - x_j[0],
                    x_i[1] - x_j[1],
                    x_i[2] - x_j[2]
                };
                energy[n] -= (G * m_i * m_j / vec_norm_3d(r_ij));
            }
        }
    }
}

WIN32DLL_API void compute_linear_momentum_python(
    double *__restrict linear_momentum,
    const double *__restrict sol_state,
    const int32 num_snapshots,
    const int32 num_particles
)
{
    const int size_snapshot = num_particles * 7;
    for (int n = 0; n < num_snapshots; n++)
    {
        linear_momentum[n] = 0.0;
        for (int i = 0; i < num_particles; i++)
        {
            const double m_i = sol_state[n * size_snapshot + i * 7 + 0];
            const double v_i[3] = {
                sol_state[n * size_snapshot + i * 7 + 4],
                sol_state[n * size_snapshot + i * 7 + 5],
                sol_state[n * size_snapshot + i * 7 + 6]
            };
            
            linear_momentum[n] += m_i * vec_sum_3d(v_i);
        }
    }
}

WIN32DLL_API void compute_angular_momentum_python(
    double *__restrict angular_momentum,
    const double *__restrict sol_state,
    const int32 num_snapshots,
    const int32 num_particles
)
{
    const int size_snapshot = num_particles * 7;
    for (int n = 0; n < num_snapshots; n++)
    {
        angular_momentum[n] = 0.0;
        for (int i = 0; i < num_particles; i++)
        {
            const double m_i = sol_state[n * size_snapshot + i * 7 + 0];
            const double x_i[3] = {
                sol_state[n * size_snapshot + i * 7 + 1],
                sol_state[n * size_snapshot + i * 7 + 2],
                sol_state[n * size_snapshot + i * 7 + 3]
            };
            const double v_i[3] = {
                sol_state[n * size_snapshot + i * 7 + 4],
                sol_state[n * size_snapshot + i * 7 + 5],
                sol_state[n * size_snapshot + i * 7 + 6]
            };

            // L = m * r x v
            double angular_momentum_vec_step[3] = {
                x_i[1] * v_i[2] - x_i[2] * v_i[1],
                x_i[2] * v_i[0] - x_i[0] * v_i[2],
                x_i[0] * v_i[1] - x_i[1] * v_i[0]
            };
            angular_momentum[n] += m_i * vec_sum_3d(angular_momentum_vec_step);
        }
    }
}
