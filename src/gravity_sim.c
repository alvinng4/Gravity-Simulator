#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gravity_sim.h"




ErrorStatus launch_simulation(
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

    /* Check system parameters */
    error_status = WRAP_TRACEBACK(finalize_system(system));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check acceleration parameters */
    error_status = WRAP_TRACEBACK(finalize_acceleration_param(acceleration_param));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check integrator parameters */
    error_status = WRAP_TRACEBACK(finalize_integration_param(integrator_param));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check output parameters */
    error_status = WRAP_TRACEBACK(finalize_output_param(output_param, settings));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check tf */
    if (tf < 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "tf must be non-negative. Got: %g",
            tf
        );
    }

    /* Remove invalid particles */
    if (settings->remove_invalid_particles)
    {
        error_status = WRAP_TRACEBACK(check_and_remove_invalid_particles(system, settings));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            return error_status;
        }
    }

    return WRAP_TRACEBACK(integrator_launch_simulation(
        system,
        integrator_param,
        acceleration_param,
        output_param,
        simulation_status,
        settings,
        tf
    ));
}

ErrorStatus launch_cosmological_simulation(
    CosmologicalSystem *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double a_begin,
    const double a_final,
    const int pm_grid_size
)
{
#ifndef USE_FFTW3
    return WRAP_RAISE_ERROR(
        GRAV_VALUE_ERROR,
        "FFTW3 are required for cosmological simulations"
    );
#else
    ErrorStatus error_status;

    /* Check system parameters */
    error_status = WRAP_TRACEBACK(finalize_cosmological_system(system));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check acceleration parameters */
    error_status = WRAP_TRACEBACK(finalize_acceleration_param(acceleration_param));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check integrator parameters */
    error_status = WRAP_TRACEBACK(finalize_integration_param(integrator_param));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check output parameters */
    error_status = WRAP_TRACEBACK(finalize_output_param(output_param, settings));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check a_final */
    if (a_final < 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "a_final must be non-negative. Got: %g",
            a_final
        );
    }

    /* Remove invalid particles */
    if (settings->remove_invalid_particles)
    {
        System temp_system = get_new_system();
        temp_system.num_particles = system->num_particles;
        temp_system.x = system->x;
        temp_system.v = system->v;
        temp_system.m = system->m;
        temp_system.particle_ids = system->particle_ids;
        error_status = WRAP_TRACEBACK(check_and_remove_invalid_particles(&temp_system, settings));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            return error_status;
        }
    }

    return WRAP_TRACEBACK(leapfrog_cosmology(
        system,
        integrator_param,
        acceleration_param,
        output_param,
        simulation_status,
        settings,
        a_begin,
        a_final,
        pm_grid_size
    ));
#endif
}