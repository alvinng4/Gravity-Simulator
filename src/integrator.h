/**
 * \file integrator.h
 * \brief Functions for the integrators
 * 
 * \author Ching-Yin NG
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "acceleration.h"
#include "common.h"
#include "error.h"
#include "settings.h"
#include "system.h"

#define INTEGRATOR_EULER 1
#define INTEGRATOR_EULER_CROMER 2
#define INTEGRATOR_RK4 3
#define INTEGRATOR_LEAPFROG 4
#define INTEGRATOR_RKF45 5
#define INTEGRATOR_DOPRI 6
#define INTEGRATOR_DVERK 7
#define INTEGRATOR_RKF78 8
#define INTEGRATOR_IAS15 9
#define INTEGRATOR_WHFAST 10

/**
 * \brief Get a new integrator parameter struct
 * 
 * \return IntegratorParam
 */
IntegratorParam get_new_integrator_param(void);

/**
 * \brief Finalize the integrator parameters
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_POINTER_ERROR If integrator_param is NULL
 * \exception GRAV_VALUE_ERROR If integrator value is invalid
 * \exception GRAV_VALUE_ERROR If integrator is fixed step size integrators and dt is not positive
 * \exception GRAV_VALUE_ERROR If integrator is adaptive step size integrators and tolerance is not positive
 */
ErrorStatus finalize_integration_param(IntegratorParam *__restrict integration_param);

/**
 * \brief Launch the simulation with the specified integrator
 * 
 * \details This function launches the simulation with the specified integrator.
 * This function will be called by the main function. User should not
 * call this function directly as it will bypass the parameter checking.
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param output_param Pointer to the output parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param tf Simulation time
 */
ErrorStatus integrator_launch_simulation(
    System *__restrict system,
    IntegratorParam *__restrict integrator_param,
    AccelerationParam *__restrict acceleration_param,
    OutputParam *__restrict output_param,
    SimulationStatus *__restrict simulation_status,
    Settings *__restrict settings,
    const double tf
);

/**
 * \brief RK Embedded integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \return ErrorStatus
 */
ErrorStatus rk_embedded(
    System *__restrict system,
    IntegratorParam *__restrict integrator_param,
    AccelerationParam *__restrict acceleration_param,
    OutputParam *__restrict output_param,
    SimulationStatus *__restrict simulation_status,
    Settings *__restrict settings,
    const double tf
);

/**
 * \brief IAS15 integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \return ErrorStatus
 */
ErrorStatus ias15(
    System *__restrict system,
    IntegratorParam *__restrict integrator_param,
    AccelerationParam *__restrict acceleration_param,
    OutputParam *__restrict output_param,
    SimulationStatus *__restrict simulation_status,
    Settings *__restrict settings,
    const double tf
);

/**
 * \brief WHFast integrator
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \return ErrorStatus
 */
ErrorStatus whfast(
    System *__restrict system,
    IntegratorParam *__restrict integrator_param,
    AccelerationParam *__restrict acceleration_param,
    OutputParam *__restrict output_param,
    SimulationStatus *__restrict simulation_status,
    Settings *__restrict settings,
    const double tf
);

/**
 * \brief Leapfrog integrator for cosmological simulations.
 * 
 * \param system Pointer to the cosmological system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param output_param Pointer to the output parameters
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param a_begin Initial scale factor
 * \param a_final Final scale factor
 * \param pm_grid_size Size of the PM grid
 * 
 * \return ErrorStatus
 */
ErrorStatus leapfrog_cosmology(
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

#endif
