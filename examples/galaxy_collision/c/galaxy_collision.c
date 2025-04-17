#include <math.h>
#include <stdio.h>
#include <hdf5.h>

#include "grav_sim.h"

#define INIT_CONDITION_FILE "../../galaxy_collision_IC.hdf5"

#define G_cgs 6.67430e-8  // cm^3 g^-1 s^-2

#define ACC_METHOD ACCELERATION_METHOD_BARNES_HUT
#define OPENING_ANGLE 0.5
#define SOFTENING_LENGTH 0.0

#define TF 4000000000 * 365.24 * 24 * 3600 // 4 billion years

#define INTEGRATOR INTEGRATOR_LEAPFROG
#define DT TF / 2000

#define OUTPUT_METHOD OUTPUT_METHOD_HDF5
#define OUTPUT_INTERVAL TF / 500 // 500 snapshots
#define OUTPUT_DIR "snapshots/"
#define OUTPUT_INITIAL true
#define COORDINATE_OUTPUT_DTYPE OUTPUT_DTYPE_FLOAT
#define VELOCITY_OUTPUT_DTYPE OUTPUT_DTYPE_FLOAT
#define MASS_OUTPUT_DTYPE OUTPUT_DTYPE_FLOAT

ErrorStatus read_init_condition(
    System *__restrict system,
    double *__restrict unit_length,
    double *__restrict unit_mass,
    double *__restrict unit_time
);

int main(void)
{
    /* Acceleration parameters */
    AccelerationParam acceleration_param = get_new_acceleration_param();
    acceleration_param.method = ACCELERATION_METHOD_BARNES_HUT;
    acceleration_param.opening_angle = OPENING_ANGLE;
    acceleration_param.softening_length = SOFTENING_LENGTH;

    System system;
    double unit_length;
    double unit_mass;
    double unit_time;
    ErrorStatus error_status = WRAP_TRACEBACK(read_init_condition(
        &system, &unit_length, &unit_mass, &unit_time
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error;
    }

    /* Convert units */
    system.G = G_cgs * (
        unit_mass
        * unit_time * unit_time
        / (unit_length * unit_length * unit_length)
    );

    /* Integrator parameters */
    IntegratorParam integrator_param = get_new_integrator_param();
    integrator_param.integrator = INTEGRATOR;
    integrator_param.dt = DT / unit_time;

    /* Output parameters */
    OutputParam output_param = get_new_output_param();
    output_param.method = OUTPUT_METHOD;
    output_param.output_dir = OUTPUT_DIR;
    output_param.output_interval = OUTPUT_INTERVAL  / unit_time;
    output_param.output_initial = OUTPUT_INITIAL;
    output_param.coordinate_output_dtype = COORDINATE_OUTPUT_DTYPE;
    output_param.velocity_output_dtype = VELOCITY_OUTPUT_DTYPE;
    output_param.mass_output_dtype = MASS_OUTPUT_DTYPE;

    /* Simulation status */
    SimulationStatus simulation_status;

    /* Settings */
    Settings settings = get_new_settings();
    settings.verbose = GRAV_VERBOSITY_NORMAL;

    error_status = WRAP_TRACEBACK(launch_simulation(
        &system,
        &integrator_param,
        &acceleration_param,
        &output_param,
        &simulation_status,
        &settings,
        TF / unit_time
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        free_system(&system);
        goto error;
    }

    /* Free memory */
    free_system(&system);

    return 0;

error:
    print_and_free_traceback(&error_status);
    return 1;
}

ErrorStatus read_init_condition(
    System *__restrict system,
    double *__restrict unit_length,
    double *__restrict unit_mass,
    double *__restrict unit_time
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
    hid_t part_type_2 = H5Gopen(file, "/PartType2", H5P_DEFAULT);
    if (header == H5I_INVALID_HID 
        || units == H5I_INVALID_HID
        || part_type_1 == H5I_INVALID_HID
        || part_type_2 == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to open groups in initial condition file");
        goto err_open_group;
    }

    /* Read header attributes */
    hid_t num_particles_attr = H5Aopen(header, "NumPart_ThisFile", H5P_DEFAULT);
    if (num_particles_attr == H5I_INVALID_HID)
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
    
    error_status = WRAP_TRACEBACK(get_initialized_system(
        system,
        temp_num_particles_arr[1] + temp_num_particles_arr[2]
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto err_initialize_system;
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

    if (H5Aread(unit_length_attr, H5T_NATIVE_DOUBLE, unit_length) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read unit length from units");
        goto err_units_attr;
    }
    if (H5Aread(unit_mass_attr, H5T_NATIVE_DOUBLE, unit_mass) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read unit mass from units");
        goto err_units_attr;
    }
    if (H5Aread(unit_time_attr, H5T_NATIVE_DOUBLE, unit_time) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read unit time from units");
        goto err_units_attr;
    }

    /* Open particle datasets */
    hid_t part_type_1_masses = H5Dopen(part_type_1, "Masses", H5P_DEFAULT);
    hid_t part_type_1_coordinates = H5Dopen(part_type_1, "Coordinates", H5P_DEFAULT);
    hid_t part_type_1_velocities = H5Dopen(part_type_1, "Velocities", H5P_DEFAULT);
    hid_t part_type_2_masses = H5Dopen(part_type_2, "Masses", H5P_DEFAULT);
    hid_t part_type_2_coordinates = H5Dopen(part_type_2, "Coordinates", H5P_DEFAULT);
    hid_t part_type_2_velocities = H5Dopen(part_type_2, "Velocities", H5P_DEFAULT);
    if (
        part_type_1_masses == H5I_INVALID_HID
        || part_type_1_coordinates == H5I_INVALID_HID
        || part_type_1_velocities == H5I_INVALID_HID
        || part_type_2_masses == H5I_INVALID_HID
        || part_type_2_coordinates == H5I_INVALID_HID
        || part_type_2_velocities == H5I_INVALID_HID
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
    if (H5Dread(part_type_2_masses, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(system->m[temp_num_particles_arr[1]])) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read masses from PartType2");
        goto err_read_particle_data;
    }
    if (H5Dread(part_type_2_coordinates, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(system->x[temp_num_particles_arr[1] * 3])) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read coordinates from PartType2");
        goto err_read_particle_data;
    }
    if (H5Dread(part_type_2_velocities, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(system->v[temp_num_particles_arr[1] * 3])) < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to read velocities from PartType2");
        goto err_read_particle_data;
    }

    /* Close HDF5 objects */
    H5Dclose(part_type_1_masses);
    H5Dclose(part_type_1_coordinates);
    H5Dclose(part_type_1_velocities);
    H5Dclose(part_type_2_masses);
    H5Dclose(part_type_2_coordinates);
    H5Dclose(part_type_2_velocities);

    H5Aclose(unit_length_attr);
    H5Aclose(unit_mass_attr);
    H5Aclose(unit_time_attr);

    H5Aclose(num_particles_attr);

    H5Gclose(header);
    H5Gclose(units);
    H5Gclose(part_type_1);
    H5Gclose(part_type_2);

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
    free_system(system);
err_initialize_system:
err_header_attr:
    H5Aclose(num_particles_attr);
err_open_group:
    H5Gclose(header);
    H5Gclose(units);
    H5Gclose(part_type_1);
    H5Gclose(part_type_2);
err_open_file:
    H5Fclose(file);

    return error_status;
}
