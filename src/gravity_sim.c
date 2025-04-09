#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gravity_sim.h"


WIN32DLL_API ErrorStatus launch_simulation(
    System *__restrict system,
    IntegratorParam *__restrict integrator_param,
    AccelerationParam *__restrict acceleration_param,
    OutputParam *__restrict output_param,
    SimulationStatus *__restrict simulation_status,
    Settings *__restrict settings,
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

WIN32DLL_API ErrorStatus launch_cosmological_simulation(
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
#ifndef USE_FFTW3
    (void) system;
    (void) integrator_param;
    (void) acceleration_param;
    (void) output_param;
    (void) simulation_status;
    (void) settings;
    (void) a_begin;
    (void) a_final;
    (void) pm_grid_size;
    return WRAP_RAISE_ERROR(
        GRAV_VALUE_ERROR,
        "FFTW3 are required for cosmological simulations. Please recompile with FFTW3 support."
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
    if (acceleration_param->method != ACCELERATION_METHOD_PM)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Only particle-mesh acceleration is supported for cosmological simulations. Got: %d",
            acceleration_param->method
        );
    }
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