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

    /* Check acceleration parameters */
    error_status = WRAP_TRACEBACK(finalize_acceleration_param(acceleration_param));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check output parameters */
    error_status = WRAP_TRACEBACK(finalize_output_param(output_param));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
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
