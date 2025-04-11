#ifndef GRAV_SIM_H
#define GRAV_SIM_H

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


/* Project version */
#ifndef VERSION_INFO
#define VERSION_INFO "unknown"
#endif

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
    System *__restrict system,
    IntegratorParam *__restrict integrator_param,
    AccelerationParam *__restrict acceleration_param,
    OutputParam *__restrict output_param,
    SimulationStatus *__restrict simulation_status,
    Settings *__restrict settings,
    const double tf
);

/**
 * \brief Main function to launch a cosmological simulation.
 * 
 * \param system Pointer to the system.
 * \param integrator_param Pointer to the integrator parameters.
 * \param acceleration_param Pointer to the acceleration parameters.
 * \param output_param Pointer to the output parameters.
 * \param simulation_status Pointer to the simulation status.
 * \param settings Pointer to the settings.
 * \param a_begin Initial scale factor.
 * \param a_final Final scale factor. 
 * \param pm_grid_size Particle mesh grid size.
 * 
 * \return Error status.
 */
ErrorStatus launch_cosmological_simulation(
    CosmologicalSystem *__restrict system,
    IntegratorParam *__restrict integrator_param,
    AccelerationParam *__restrict acceleration_param,
    OutputParam *__restrict output_param,
    SimulationStatus *__restrict simulation_status,
    Settings *__restrict settings,
    const double a_begin,
    const double a_final,
    const int pm_grid_size
);

/**
 * \brief Get the logo string of grav_sim.
 * 
 * \return Pointer to the logo string.
 */
const char* get_grav_sim_logo_string(void);

/**
 * \brief Print project compilation information.
 */
void print_compilation_info(void);

#endif
