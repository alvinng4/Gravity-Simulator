#ifndef GRAVITY_SIM_H
#define GRAVITY_SIM_H

/* Acceleration */
#include "acceleration.h"

/* Common definitions */
#include "common.h"

/* Exception handling */
#include "error.h"

/* Integrator */
#include "integrator.h"

/* Output */
#include "output.h"

/* System */
#include "system.h"

/* Settings */
#include "settings.h"

/* Utils */
#include "utils.h"


/**
 * \brief Main function to launch a simulation.
 * 
 * \param system Pointer to the system.
 * \param integrator_param Pointer to the integrator parameters.
 * \param acceleration_param Pointer to the acceleration parameters.
 * \param output_param Pointer to the output parameters.
 * \param simulation_status Pointer to the simulation status.
 * \param settings Pointer to the settings.
 * \param tf Simulation time.
 * 
 * \return Error status.
 */
ErrorStatus launch_simulation(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double tf
);


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
);






#endif
