#include <math.h>
#include <stdio.h>
#include <hdf5.h>

#include "gravity_sim.h"

#define INIT_CONDITION_FILE "../../ics_swift.hdf5"
#define A_FINAL 1.0 // Scale factor at the end of simulation

ErrorStatus read_init_condition(
    CosmologicalSystem *__restrict system,
    double *__restrict softening_length,
    int *__restrict pm_grid_size,
    double *__restrict a_begin
);

int main(void)
{
    /* Acceleration parameters */
    AccelerationParam acceleration_param = get_new_acceleration_param();
    acceleration_param.method = ACCELERATION_METHOD_PM;
    // acceleration_param.opening_angle = 0.5;

    CosmologicalSystem system;
    int pm_grid_size;
    double a_begin;
    ErrorStatus error_status = WRAP_TRACEBACK(read_init_condition(
        &system,
        &acceleration_param.softening_length,
        &pm_grid_size,
        &a_begin
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error;
    }

    /* Integrator parameters */
    IntegratorParam integrator_param = get_new_integrator_param();
    integrator_param.integrator = INTEGRATOR_LEAPFROG;
    integrator_param.dt = (A_FINAL - a_begin) / 1000.0;

    /* Output parameters */
    OutputParam output_param = get_new_output_param();
    output_param.method = OUTPUT_METHOD_HDF5;
    output_param.output_dir = "snapshots/";
    output_param.output_interval = (A_FINAL - a_begin) / 100.0;
    output_param.output_initial = true;
    output_param.coordinate_output_dtype = OUTPUT_DTYPE_FLOAT;
    output_param.velocity_output_dtype = OUTPUT_DTYPE_FLOAT;
    output_param.mass_output_dtype = OUTPUT_DTYPE_FLOAT;

    /* Simulation status */
    SimulationStatus simulation_status;

    /* Settings */
    Settings settings = get_new_settings();
    settings.verbose = GRAV_VERBOSITY_NORMAL;

    error_status = WRAP_TRACEBACK(launch_cosmological_simulation(
        &system,
        &integrator_param,
        &acceleration_param,
        &output_param,
        &simulation_status,
        &settings,
        a_begin,
        A_FINAL,
        pm_grid_size
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        free_cosmological_system(&system);
        goto error;
    }

    /* Free memory */
    free_cosmological_system(&system);

    return 0;

error:
    print_and_free_traceback(&error_status);
    return 1;
}

ErrorStatus read_init_condition(
    CosmologicalSystem *__restrict system,
    double *__restrict softening_length,
    int *__restrict pm_grid_size,
    double *__restrict a_begin
)
{
    ErrorStatus error_status;

    /* Open file */
    hid_t file = H5Fopen(INIT_CONDITION_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file == H5I_INVALID_HID)    
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to open initial condition file");
        goto err_open_file;
    }

    /* Open groups */
    hid_t header = H5Gopen(file, "/Header", H5P_DEFAULT);
    hid_t units = H5Gopen(file, "/Units", H5P_DEFAULT);
    hid_t part_type_1 = H5Gopen(file, "/PartType1", H5P_DEFAULT);
    if (header == H5I_INVALID_HID 
        || units == H5I_INVALID_HID
        || part_type_1 == H5I_INVALID_HID)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to open groups in initial condition file");
        goto err_open_group;
    }

    /* Read header attributes */
    hid_t num_particles_attr = H5Aopen(header, "NumPart_ThisFile", H5P_DEFAULT);
    hid_t box_size_attr = H5Aopen(header, "BoxSize", H5P_DEFAULT);
    hid_t softening_length_attr = H5Aopen(header, "suggested_highressoft", H5P_DEFAULT);
    hid_t grid_size_attr = H5Aopen(header, "suggested_pmgrid", H5P_DEFAULT);
    hid_t omega_m_attr = H5Aopen(header, "Omega0", H5P_DEFAULT);
    hid_t omega_lambda_attr = H5Aopen(header, "OmegaLambda", H5P_DEFAULT);
    hid_t redshift_attr = H5Aopen(header, "Redshift", H5P_DEFAULT);
    hid_t h_attr = H5Aopen(header, "HubbleParam", H5P_DEFAULT);
    hid_t mass_table_attr = H5Aopen(header, "MassTable", H5P_DEFAULT);
    if (
        num_particles_attr == H5I_INVALID_HID
        || box_size_attr == H5I_INVALID_HID
        || softening_length_attr == H5I_INVALID_HID
        || grid_size_attr == H5I_INVALID_HID
        || omega_m_attr == H5I_INVALID_HID
        || omega_lambda_attr == H5I_INVALID_HID
        || redshift_attr == H5I_INVALID_HID
        || h_attr == H5I_INVALID_HID
        || mass_table_attr == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to open attributes in header");
        goto err_header_attr;
    }

    uint32_t temp_num_particles_arr[6];
    if (H5Aread(num_particles_attr, H5T_NATIVE_UINT32, &temp_num_particles_arr) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read number of particles from header");
        goto err_header_attr;
    }
    double box_size;
    if (H5Aread(box_size_attr, H5T_NATIVE_DOUBLE, &box_size) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read box size from header");
        goto err_header_attr;
    }
    if (H5Aread(softening_length_attr, H5T_NATIVE_DOUBLE, softening_length) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read softening length from header");
        goto err_header_attr;
    }
    uint32_t temp_grid_size;
    if (H5Aread(grid_size_attr, H5T_NATIVE_UINT32, &temp_grid_size) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read grid size from header");
        goto err_header_attr;
    }
    *pm_grid_size = temp_grid_size;

    double omega_m;
    if (H5Aread(omega_m_attr, H5T_NATIVE_DOUBLE, &omega_m) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read omega_m from header");
        goto err_header_attr;
    }
    double omega_lambda;
    if (H5Aread(omega_lambda_attr, H5T_NATIVE_DOUBLE, &omega_lambda) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read omega_lambda from header");
        goto err_header_attr;
    }
    if (H5Aread(redshift_attr, H5T_NATIVE_DOUBLE, a_begin) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read redshift from header");
        goto err_header_attr;
    }
    *a_begin = 1.0 / (*a_begin + 1.0);

    double h0;
    if (H5Aread(h_attr, H5T_NATIVE_DOUBLE, &h0) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read Hubble parameter from header");
        goto err_header_attr;
    }
    double temp_mass_arr[6];
    if (H5Aread(mass_table_attr, H5T_NATIVE_DOUBLE, &temp_mass_arr) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read mass table from header");
        goto err_header_attr;
    }
    
    error_status = WRAP_TRACEBACK(get_initialized_cosmological_system(
        system,
        temp_num_particles_arr[1]
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto err_initialize_system;
    }
    system->box_width = box_size / 2.0;
    system->box_center[0] = system->box_width;
    system->box_center[1] = system->box_width;
    system->box_center[2] = system->box_width;
    system->omega_m = omega_m;
    system->omega_lambda = omega_lambda;
    system->h0 = h0;
    for (int i = 0; i < system->num_particles; i++)
    {
        system->m[i] = temp_mass_arr[1];
    }

    /* Read units attributes */
    hid_t unit_length_attr = H5Aopen(units, "Unit length in cgs (U_L)", H5P_DEFAULT);
    hid_t unit_mass_attr = H5Aopen(units, "Unit mass in cgs (U_M)", H5P_DEFAULT);
    hid_t unit_time_attr = H5Aopen(units, "Unit time in cgs (U_t)", H5P_DEFAULT);

    if (
        unit_length_attr == H5I_INVALID_HID
        || unit_mass_attr == H5I_INVALID_HID
        || unit_time_attr == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to open attributes in units");
        goto err_units_attr;
    }

    if (H5Aread(unit_length_attr, H5T_NATIVE_DOUBLE, &(system->unit_length)) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read unit length from units");
        goto err_units_attr;
    }
    if (H5Aread(unit_mass_attr, H5T_NATIVE_DOUBLE, &(system->unit_mass)) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read unit mass from units");
        goto err_units_attr;
    }
    if (H5Aread(unit_time_attr, H5T_NATIVE_DOUBLE, &(system->unit_time)) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read unit time from units");
        goto err_units_attr;
    }

    /* Open particle datasets */
    hid_t part_type_1_masses = H5Dopen(part_type_1, "Masses", H5P_DEFAULT);
    hid_t part_type_1_coordinates = H5Dopen(part_type_1, "Coordinates", H5P_DEFAULT);
    hid_t part_type_1_velocities = H5Dopen(part_type_1, "Velocities", H5P_DEFAULT);
    if (
        part_type_1_masses == H5I_INVALID_HID
        || part_type_1_coordinates == H5I_INVALID_HID
        || part_type_1_velocities == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to open datasets in PartType1");
        goto err_open_particle_dataset;
    }

    /* Read particle data */
    if (H5Dread(part_type_1_masses, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, system->m) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read masses from PartType1");
        goto err_read_particle_data;
    }
    if (H5Dread(part_type_1_coordinates, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, system->x) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read coordinates from PartType1");
        goto err_read_particle_data;
    }
    if (H5Dread(part_type_1_velocities, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, system->v) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read velocities from PartType1");
        goto err_read_particle_data;
    }

    /* Close HDF5 objects */
    H5Dclose(part_type_1_masses);
    H5Dclose(part_type_1_coordinates);
    H5Dclose(part_type_1_velocities);

    H5Aclose(unit_length_attr);
    H5Aclose(unit_mass_attr);
    H5Aclose(unit_time_attr);

    H5Aclose(num_particles_attr);
    H5Aclose(box_size_attr);
    H5Aclose(softening_length_attr);
    H5Aclose(grid_size_attr);
    H5Aclose(omega_m_attr);
    H5Aclose(omega_lambda_attr);
    H5Aclose(redshift_attr);
    H5Aclose(h_attr);
    H5Aclose(mass_table_attr);

    H5Gclose(header);
    H5Gclose(units);
    H5Gclose(part_type_1);

    H5Fclose(file);

    return make_success_error_status();

err_read_particle_data:
err_open_particle_dataset:
    H5Dclose(part_type_1_masses);
    H5Dclose(part_type_1_coordinates);
    H5Dclose(part_type_1_velocities);
err_units_attr:
    H5Aclose(unit_length_attr);
    H5Aclose(unit_mass_attr);
    H5Aclose(unit_time_attr);
    free_cosmological_system(system);
err_initialize_system:
err_header_attr:
    H5Aclose(num_particles_attr);
    H5Aclose(box_size_attr);
    H5Aclose(softening_length_attr);
    H5Aclose(grid_size_attr);
    H5Aclose(omega_m_attr);
    H5Aclose(omega_lambda_attr);
    H5Aclose(redshift_attr);
    H5Aclose(h_attr);
    H5Aclose(mass_table_attr);
err_open_group:
    H5Gclose(header);
    H5Gclose(units);
    H5Gclose(part_type_1);
err_open_file:
    H5Fclose(file);

    return error_status;
}
