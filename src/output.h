/**
 * \file output.h
 * \brief Functions for outputting the simulation results
 * 
 * \author Ching-Yin Ng
 * \date March 2025
 */

#ifndef OUTPUT_H
#define OUTPUT_H

/* Output methods */
#define OUTPUT_METHOD_DISABLED 1
#define OUTPUT_METHOD_CSV 2
#define OUTPUT_METHOD_HDF5 3

/* Output data types */
#define OUTPUT_DTYPE_FLOAT 1
#define OUTPUT_DTYPE_DOUBLE 2

#include <stdbool.h>

#include "common.h"
#include "error.h"
#include "settings.h"
#include "system.h"


/**
 * \brief Get a new output parameter structure
 * 
 * \return OutputParam
 */
OutputParam get_new_output_param(void);

/**
 * \brief Finalize and check the output parameters.
 * 
 * \param output_param Pointer to the output parameters.
 * \param[in] settings Pointer to the settings.
 * 
 * \return ErrorStatus
 */
ErrorStatus finalize_output_param(
    OutputParam *restrict output_param,
    const Settings *restrict settings
);

/**
 * \brief Output a snapshot of the simulation.
 * 
 * \param system Pointer to the gravitational system.
 * \param integrator_param Pointer to the integrator parameters.
 * \param acceleration_param Pointer to the acceleration parameters.
 * \param output_param Pointer to the output parameters.
 * \param simulation_status Pointer to the simulation status.
 * \param settings Pointer to the settings.
 * 
 * \return ErrorStatus
 */
ErrorStatus output_snapshot(
    OutputParam *restrict output_param,
    const System *restrict system,
    const IntegratorParam *restrict integrator_param,
    const AccelerationParam *restrict acceleration_param,
    const SimulationStatus *restrict simulation_status,
    const Settings *restrict settings
);

/**
 * \brief Output a snapshot of the cosmological simulation in CSV format.
 * 
 * \param output_param Pointer to the output parameters.
 * \param system Pointer to the gravitational system.
 * \param simulation_status Pointer to the simulation status.
 * \param settings Pointer to the settings.
 * 
 * \return ErrorStatus
 */
ErrorStatus output_snapshot_cosmology(
    OutputParam *restrict output_param,
    const CosmologicalSystem *restrict system,
    const SimulationStatus *restrict simulation_status,
    const Settings *restrict settings
);

#endif
