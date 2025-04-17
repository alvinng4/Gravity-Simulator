#include <stdio.h>

#include "grav_sim.h"

#define NUM_RUNS_PAIRWISE 0
#define NUM_RUNS_BARNES_HUT 1

#define NUM_PARTICLES 10000000
#define X_MAX 1.0
#define X_MIN -1.0

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

    /* Initialize system parameters */
    double *__restrict m = system.m;
    double *__restrict x = system.x;
    pcg32_random_t rng = init_pcg_rng(); // Random number generator
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        m[i] = 1.0;
        x[i * 3 + 0] = grav_randrange(X_MIN, X_MAX, &rng);
        x[i * 3 + 1] = grav_randrange(X_MIN, X_MAX, &rng);
        x[i * 3 + 2] = grav_randrange(X_MIN, X_MAX, &rng);
    }

    /* Acceleration parameters */
    AccelerationParam acceleration_params[2];
    const int num_times_acceleration_param[2] = {NUM_RUNS_PAIRWISE, NUM_RUNS_BARNES_HUT};

    acceleration_params[0] = get_new_acceleration_param();
    acceleration_params[0].method = ACCELERATION_METHOD_PAIRWISE;

    acceleration_params[1] = get_new_acceleration_param();
    acceleration_params[1].method = ACCELERATION_METHOD_BARNES_HUT;
    acceleration_params[1].opening_angle = 0.5;
    acceleration_params[1].max_num_particles_per_leaf = 1;

    error_status = WRAP_TRACEBACK(benchmark_acceleration(
        &system,
        acceleration_params,
        2,
        num_times_acceleration_param
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error;
    }

    /* Free memory */
    free_system(&system);

    return 0;

error:
    print_and_free_traceback(&error_status);
    return 1;
}
