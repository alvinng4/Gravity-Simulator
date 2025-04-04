#include <math.h>
#include <stdio.h>

#include "gravity_sim.h"

#define NUM_PARTICLES 25000
#define DT 180.0
#define TF 5000000 * 365.24
#define OUTPUT_INTERVAL 2500 * 365.24 // Store 2000 snapshots
#define OUTPUT_METHOD OUTPUT_METHOD_HDF5

int main(void)
{
    ErrorStatus error_status;

    /* Initialize system */
    System system;
    error_status = WRAP_TRACEBACK(get_initialized_system(&system, NUM_PARTICLES));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error;
    }

    error_status = WRAP_TRACEBACK(initialize_built_in_system(
        &system,
        "solar_system",
        true
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error;
    }

    /* Sun, Mars, Jupiter, Saturn (1 <- 4, 2 <- 5, 3 <- 6) */
    const int num_massive_objects = 4;
    for (int i = 1, j = 4; i < 4; i++, j++)
    {
        system.m[i] = system.m[j];
        system.x[i * 3 + 0] = system.x[j * 3 + 0];
        system.x[i * 3 + 1] = system.x[j * 3 + 1];
        system.x[i * 3 + 2] = system.x[j * 3 + 2];
        system.v[i * 3 + 0] = system.v[j * 3 + 0];
        system.v[i * 3 + 1] = system.v[j * 3 + 1];
        system.v[i * 3 + 2] = system.v[j * 3 + 2];
    }

    /* Asteroids */
    pcg32_random_t rng = init_pcg_rng();
    for (int i = num_massive_objects; i < NUM_PARTICLES; i++)
    {
        system.m[i] = 0.0;

        const double semi_major_axis = grav_randrange(2.0, 3.35, &rng);
        const double eccentricity = grav_randrange(-0.12, 0.12, &rng);
        const double inclination = grav_randrange(-0.3, 0.3, &rng);
        const double argument_of_periapsis = grav_randrange(0.0, 2.0 * M_PI, &rng);
        const double longitude_of_ascending_node = grav_randrange(0.0, 2.0 * M_PI, &rng);
        const double true_anomaly = grav_randrange(0.0, 2.0 * M_PI, &rng);

        keplerian_to_cartesian(
            &system.x[i * 3],
            &system.v[i * 3],
            semi_major_axis,
            eccentricity,
            inclination,
            argument_of_periapsis,
            longitude_of_ascending_node,
            true_anomaly,
            system.m[0],
            system.G
        );

        system.x[i * 3 + 0] += system.x[0];
        system.x[i * 3 + 1] += system.x[1];
        system.x[i * 3 + 2] += system.x[2];
        system.v[i * 3 + 0] += system.v[0];
        system.v[i * 3 + 1] += system.v[1];
        system.v[i * 3 + 2] += system.v[2];
    }
    error_status = WRAP_TRACEBACK(system_set_center_of_mass_zero(&system));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error;
    }
    error_status = WRAP_TRACEBACK(system_set_total_momentum_zero(&system));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error;
    }

    /* Boundary condition */
    system.boundary_condition = BOUNDARY_COUNDITION_NONE;

    /* Acceleration parameters */
    AccelerationParam acceleration_param = get_new_acceleration_param();
    acceleration_param.method = ACCELERATION_METHOD_MASSLESS;

    /* Integrator parameters */
    IntegratorParam integrator_param = get_new_integrator_param();
    integrator_param.integrator = INTEGRATOR_WHFAST;
    integrator_param.dt = DT;

    /* Output parameters */
    OutputParam output_param = get_new_output_param();
    output_param.method = OUTPUT_METHOD;
    output_param.output_dir = "../snapshots/";
    output_param.output_interval = OUTPUT_INTERVAL;
    output_param.output_initial = true;
    output_param.coordinate_output_dtype = OUTPUT_DTYPE_FLOAT;
    output_param.velocity_output_dtype = OUTPUT_DTYPE_FLOAT;
    output_param.mass_output_dtype = OUTPUT_DTYPE_FLOAT;

    /* Simulation status */
    SimulationStatus simulation_status;

    /* Settings */
    Settings settings = get_new_settings();
    settings.verbose = GRAV_VERBOSITY_VERBOSE;

    printf("Launching simulation...\n");
    error_status = WRAP_TRACEBACK(launch_simulation(
        &system,
        &integrator_param,
        &acceleration_param,
        &output_param,
        &simulation_status,
        &settings,
        TF
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error;
    }
    printf("Done!\n");

    /* Free memory */
    free_system(&system);

    return 0;

error:
    print_and_free_traceback(&error_status);
    return 1;
}
