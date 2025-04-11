/**
 * \file output.c
 * \brief Function definitions for simulation output
 * 
 * \author Ching-Yin Ng
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* HDF5 */
#ifdef USE_HDF5
#include <hdf5.h>
#endif

/* For mkdir */
#ifndef _WIN32
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <direct.h>
#include <windows.h>
#endif

#include "acceleration.h"
#include "common.h"
#include "error.h"
#include "integrator.h"
#include "output.h"
#include "system.h"

/**
 * \brief Output a snapshot of the simulation in CSV format.
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
IN_FILE ErrorStatus output_snapshot_csv(
    OutputParam *output_param,
    const System *system,
    const IntegratorParam *integrator_param,
    const AccelerationParam *acceleration_param,
    const SimulationStatus *simulation_status,
    const Settings *settings
);

#ifdef USE_HDF5
/**
 * \brief Output a snapshot of the simulation in HDF5 format.
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
IN_FILE ErrorStatus output_snapshot_hdf5(
    OutputParam *output_param,
    const System *system,
    const IntegratorParam *integrator_param,
    const AccelerationParam *acceleration_param,
    const SimulationStatus *simulation_status,
    const Settings *settings
);

IN_FILE ErrorStatus output_snapshot_cosmology_hdf5(
    OutputParam *output_param,
    const CosmologicalSystem *system,
    const IntegratorParam *integrator_param,
    const AccelerationParam *acceleration_param,
    const SimulationStatus *simulation_status,
    const Settings *settings
);
#endif

OutputParam get_new_output_param(void)
{
    OutputParam output_param = {
        .method = OUTPUT_METHOD_DISABLED,
        .output_dir = NULL,
        .output_initial = false,
        .output_interval = -1.0,
        .coordinate_output_dtype = OUTPUT_DTYPE_DOUBLE,
        .velocity_output_dtype = OUTPUT_DTYPE_DOUBLE,
        .mass_output_dtype = OUTPUT_DTYPE_DOUBLE,
        .output_count_ = 0
    };
    return output_param;
}

IN_FILE ErrorStatus check_output_method(const int output_method)
{
    switch (output_method)
    {
        case OUTPUT_METHOD_DISABLED:
        case OUTPUT_METHOD_CSV:
            break;
        case OUTPUT_METHOD_HDF5:
#ifdef USE_HDF5
            break;
#else
            return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "HDF5 output method is not available");
#endif
        default:
        {
            return WRAP_RAISE_ERROR_FMT(
                GRAV_VALUE_ERROR,
                "Unknown output method. Got: %d",
                output_method
            );
        }
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus finalize_output_param(
    OutputParam *__restrict output_param,
    const Settings *__restrict settings
)
{
    ErrorStatus error_status;

    error_status = WRAP_TRACEBACK(check_output_method(output_param->method));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    if (output_param->method == OUTPUT_METHOD_DISABLED)
    {
        return make_success_error_status();
    }

    /* Check storing interval */
    if (output_param->output_interval <= 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Output interval must be positive. Got: %.17g",
            output_param->output_interval
        );
    }

    /* Check directory path */
    if (!output_param->output_dir)
    {
        // Set to snapshots_YYYYMMDD_HHMMSS if directory path is not set
        size_t path_str_len = strlen("snapshots_YYYYMMDD_HHMMSS/") + 1;
        char *output_dir = malloc(path_str_len * sizeof(char));
        if (!output_dir)
        {
            return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for directory path.");
        }
        const time_t raw_time = time(NULL);
        struct tm *time_info = localtime(&raw_time);
        strftime(output_dir, path_str_len, "snapshots_%Y%m%d_%H%M%S/", time_info);
        output_param->output_dir = output_dir;
    }
    else
    {
        if (output_param->output_dir[strlen(output_param->output_dir) - 1] != '/')
        {
            return WRAP_RAISE_ERROR_FMT(
                GRAV_VALUE_ERROR,
                "Directory path for storing snapshots must end with a trailing slash (\"/\"). Got: \"%s\".",
                output_param->output_dir
            );
        }
    }

    /* Create directory */
#ifdef _WIN32
    if (_mkdir(output_param->output_dir) == -1)
    {
        if (GetFileAttributes(output_param->output_dir) == INVALID_FILE_ATTRIBUTES)
        {
            return WRAP_RAISE_ERROR_FMT(
                GRAV_OS_ERROR,
                "Failed to access path for storing snapshots: \"%s\".",
                output_param->output_dir
            );
        }
        else if (
            (GetFileAttributes(output_param->output_dir) & FILE_ATTRIBUTE_DIRECTORY)
            && (settings->verbose >= GRAV_VERBOSITY_IGNORE_INFO)
        )
        {
            error_status = WRAP_RAISE_WARNING_FMT(
                "Directory for storing snapshots already exists. The files will be overwritten. Directory: \"%s\".",
                output_param->output_dir
            );
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
    }
#else
    struct stat st = {0};
    if (mkdir(output_param->output_dir, 0777) == -1)
    {
        if(stat(output_param->output_dir, &st) != 0)
        {
            return WRAP_RAISE_ERROR_FMT(
                GRAV_OS_ERROR,
                "Failed to access path for storing snapshots: \"%s\".",
                output_param->output_dir
            );
        }

        else if (
            (st.st_mode & S_IFDIR)
            && (settings->verbose >= GRAV_VERBOSITY_IGNORE_INFO)
        )
        {
            error_status = WRAP_RAISE_WARNING_FMT(
                "Directory for storing snapshots already exists. The files will be overwritten. Directory: \"%s\".",
                output_param->output_dir
            );
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
    }
#endif

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus output_snapshot(
    OutputParam *output_param,
    const System *system,
    const IntegratorParam *integrator_param,
    const AccelerationParam *acceleration_param,
    const SimulationStatus *simulation_status,
    const Settings *settings
)
{
    ErrorStatus error_status = make_success_error_status();

    switch (output_param->method)
    {
        case OUTPUT_METHOD_DISABLED:
            break;
        case OUTPUT_METHOD_CSV:
            error_status = WRAP_TRACEBACK(output_snapshot_csv(
                output_param,
                system,
                integrator_param,
                acceleration_param,
                simulation_status,
                settings
            ));
            break;
        case OUTPUT_METHOD_HDF5:
#ifdef USE_HDF5
            error_status = WRAP_TRACEBACK(output_snapshot_hdf5(
                output_param,
                system,
                integrator_param,
                acceleration_param,
                simulation_status,
                settings
            ));
            break;
#else
            error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "HDF5 output method is not available");
            break;
#endif
        default:
            error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Unknown output method");
            break;
    }
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    (output_param->output_count_)++;

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus output_snapshot_cosmology(
    OutputParam *output_param,
    const CosmologicalSystem *system,
    const IntegratorParam *integrator_param,
    const AccelerationParam *acceleration_param,
    const SimulationStatus *simulation_status,
    const Settings *settings
)
{
    ErrorStatus error_status = make_success_error_status();

    switch (output_param->method)
    {
        case OUTPUT_METHOD_DISABLED:
            break;
        case OUTPUT_METHOD_CSV:
            error_status = WRAP_RAISE_ERROR(
                GRAV_VALUE_ERROR,
                "CSV output method is not supported for cosmological simulation."
            );
            break;
        case OUTPUT_METHOD_HDF5:
#ifdef USE_HDF5
            error_status = WRAP_TRACEBACK(output_snapshot_cosmology_hdf5(
                output_param,
                system,
                integrator_param,
                acceleration_param,
                simulation_status,
                settings
            ));
            break;
#else
            error_status = WRAP_RAISE_ERROR(
                GRAV_VALUE_ERROR,
                "HDF5 output method is not available. Please recompile with HDF5 support."
            );
            break;
#endif
        default:
            error_status = WRAP_RAISE_ERROR_FMT(
                GRAV_VALUE_ERROR,
                "Unknown output method. Got: %d",
                output_param->method
            );
            break;
    }
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    (output_param->output_count_)++;

    return make_success_error_status();
}

IN_FILE ErrorStatus output_snapshot_csv(
    OutputParam *__restrict output_param,
    const System *__restrict system,
    const IntegratorParam *__restrict integrator_param,
    const AccelerationParam *__restrict acceleration_param,
    const SimulationStatus *__restrict simulation_status,
    const Settings *__restrict settings
)
{
    ErrorStatus error_status;

    (void) integrator_param;
    (void) acceleration_param;
    (void) settings;

    if (!output_param->output_dir)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Output directory path is NULL.");
        goto err_output_dir_null;
    }

    /* Make file path string */
    const int file_path_length = (
        strlen(output_param->output_dir)
        + snprintf(NULL, 0, "snapshot_%05d.csv", output_param->output_count_)
        + 1  // Null terminator
    );
    char *__restrict file_path = malloc(file_path_length * sizeof(char));
    if (!file_path)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for file path string.");
        goto err_file_path_memory_alloc;
    }
    int actual_file_path_length = snprintf(file_path, file_path_length, "%ssnapshot_%05d.csv", output_param->output_dir, output_param->output_count_);

    if (actual_file_path_length < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Failed to get storing file path string");
        goto err_write_file_path_string;
    }
    else if (actual_file_path_length >= file_path_length)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Storing file path string is truncated.");
        goto err_write_file_path_string;
    }

    /* Open file */
    FILE *file = fopen(file_path, "w");
    if (!file)
    {
        error_status = WRAP_RAISE_ERROR_FMT(
            GRAV_OS_ERROR,
            "Failed to open file for storing snapshots: \"%s\".",
            file_path
        );
        goto err_open_file;
    }

    /* Write metadata */
    fprintf(file, "# num_particles: %d\n", system->num_particles);
    fprintf(file, "# G: %.17g\n", system->G);
    fprintf(file, "# time: %.17g\n", simulation_status->t);
    fprintf(file, "# dt: %.17g\n", simulation_status->dt);

    /* Write header */
    fputs("particle_id,m,x,y,z,vx,vy,vz\n", file);

    /* Write data */
    const int num_particles = system->num_particles;
    const int *__restrict particle_ids = system->particle_ids;
    const double *__restrict x = system->x;
    const double *__restrict v = system->v;
    const double *__restrict m = system->m;

    for (int i = 0; i < num_particles; i++)
    {
        // particle_id
        fprintf(file, "%d,", particle_ids[i]);

        // masses
        if (output_param->mass_output_dtype == OUTPUT_DTYPE_FLOAT)
        {
            fprintf(file, "%.8g,", m[i]);
        }
        else
        {
            fprintf(file, "%.17g,", m[i]);
        }

        // coordinates
        if (output_param->coordinate_output_dtype == OUTPUT_DTYPE_FLOAT)
        {
            fprintf(file, "%.8g,", x[3 * i]);
            fprintf(file, "%.8g,", x[3 * i + 1]);
            fprintf(file, "%.8g,", x[3 * i + 2]);
        }
        else
        {
            fprintf(file, "%.17g,", x[3 * i]);
            fprintf(file, "%.17g,", x[3 * i + 1]);
            fprintf(file, "%.17g,", x[3 * i + 2]);
        }

        // velocities
        if (output_param->velocity_output_dtype == OUTPUT_DTYPE_FLOAT)
        {
            fprintf(file, "%.8g,", v[3 * i]);
            fprintf(file, "%.8g,", v[3 * i + 1]);
            fprintf(file, "%.8g\n", v[3 * i + 2]);
        }
        else
        {
            fprintf(file, "%.17g,", v[3 * i]);
            fprintf(file, "%.17g,", v[3 * i + 1]);
            fprintf(file, "%.17g\n", v[3 * i + 2]);
        }
    }

    fflush(file);
    fclose(file);

    return make_success_error_status();

err_open_file:
err_write_file_path_string:
err_file_path_memory_alloc:
    free(file_path);
err_output_dir_null:
    return error_status;
}

#ifdef USE_HDF5
IN_FILE ErrorStatus output_snapshot_hdf5(
    OutputParam *__restrict output_param,
    const System *__restrict system,
    const IntegratorParam *__restrict integrator_param,
    const AccelerationParam *__restrict acceleration_param,
    const SimulationStatus *__restrict simulation_status,
    const Settings *__restrict settings
)
{
    ErrorStatus error_status;

    (void) integrator_param;
    (void) acceleration_param;
    (void) settings;

    if (!output_param->output_dir)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Output directory path is NULL.");
        goto err_output_dir_null;
    }

    /* Declare variables */
    const int num_particles = system->num_particles;
    const int *__restrict particle_ids = system->particle_ids;
    const double *__restrict x = system->x;
    const double *__restrict v = system->v;
    const double *__restrict m = system->m;

    /* Make file path string */
    const int file_path_length = (
        strlen(output_param->output_dir)
        + snprintf(NULL, 0, "snapshot_%05d.hdf5", output_param->output_count_)
        + 1  // Null terminator
    );
    char *__restrict file_path = malloc(file_path_length * sizeof(char));
    if (!file_path)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for file path string.");
        goto err_file_path_memory_alloc;
    }
    int actual_file_path_length = snprintf(file_path, file_path_length, "%ssnapshot_%05d.hdf5", output_param->output_dir, output_param->output_count_);

    if (actual_file_path_length < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Failed to get storing file path string");
        goto err_write_file_path_string;
    }
    else if (actual_file_path_length >= file_path_length)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Storing file path string is truncated.");
        goto err_write_file_path_string;
    }

    /* Create HDF5 file */
    hid_t file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (!file)
    {
        error_status = WRAP_RAISE_ERROR_FMT(
            GRAV_OS_ERROR,
            "Failed to create HDF5 file for storing snapshots: \"%s\".",
            file_path
        );
        goto err_create_hdf5_file;
    }

    /* Create group */
    hid_t header_group = H5Gcreate(file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t part_type_0_group = H5Gcreate(file, "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (header_group == H5I_INVALID_HID || part_type_0_group == H5I_INVALID_HID)
    {
        H5Gclose(header_group);
        H5Gclose(part_type_0_group);
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 groups.");
        goto err_create_hdf5_group;
    }

    /* Create dataspaces */
    hsize_t dims_1d_1[1] = {1};
    hsize_t dims_1d_objects_count[1] = {num_particles};
    hsize_t dims_3d_objects_count[2] = {num_particles, 3};
    hid_t dataspace_1d_1 = H5Screate_simple(1, dims_1d_1, NULL);
    hid_t dataspace_1d_objects_count = H5Screate_simple(1, dims_1d_objects_count, NULL);
    hid_t dataspace_3d_objects_count = H5Screate_simple(2, dims_3d_objects_count, NULL);
    if (
        dataspace_1d_1 == H5I_INVALID_HID
        || dataspace_1d_objects_count == H5I_INVALID_HID
        || dataspace_3d_objects_count == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 dataspace.");
        goto err_create_hdf5_dataspace;
    }
    
    /* Create attributes for header */
    hid_t header_attr_num_files_per_snapshot = H5Acreate(header_group, "NumFilesPerSnapshot", H5T_NATIVE_INT, dataspace_1d_1, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_num_part_this_file = H5Acreate(header_group, "NumPart_ThisFile", H5T_NATIVE_INT, dataspace_1d_1, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_num_part_total = H5Acreate(header_group, "NumPart_Total", H5T_NATIVE_INT, dataspace_1d_1, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_time = H5Acreate(header_group, "Time", H5T_NATIVE_DOUBLE, dataspace_1d_1, H5P_DEFAULT, H5P_DEFAULT);
    if (
        header_attr_num_files_per_snapshot == H5I_INVALID_HID
        || header_attr_num_part_this_file == H5I_INVALID_HID
        || header_attr_num_part_total == H5I_INVALID_HID
        || header_attr_time == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 attribute for header.");
        goto err_create_hdf5_header_attr;
    }

    /* Create datasets for PartType0 */
    hid_t part_type_0_dataset_part_ids = H5Dcreate(part_type_0_group, "ParticleIDs", H5T_NATIVE_INT, dataspace_1d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t part_type_0_dataset_masses;
    hid_t part_type_0_dataset_coordinates;
    hid_t part_type_0_dataset_velocities;
    
    if (output_param->mass_output_dtype == OUTPUT_DTYPE_FLOAT)
    {
        part_type_0_dataset_masses = H5Dcreate(part_type_0_group, "Masses", H5T_NATIVE_FLOAT, dataspace_1d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        part_type_0_dataset_masses = H5Dcreate(part_type_0_group, "Masses", H5T_NATIVE_DOUBLE, dataspace_1d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (output_param->coordinate_output_dtype == OUTPUT_DTYPE_FLOAT)
    {
        part_type_0_dataset_coordinates = H5Dcreate(part_type_0_group, "Coordinates", H5T_NATIVE_FLOAT, dataspace_3d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        part_type_0_dataset_coordinates = H5Dcreate(part_type_0_group, "Coordinates", H5T_NATIVE_DOUBLE, dataspace_3d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (output_param->velocity_output_dtype == OUTPUT_DTYPE_FLOAT)
    {
        part_type_0_dataset_velocities = H5Dcreate(part_type_0_group, "Velocities", H5T_NATIVE_FLOAT, dataspace_3d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        part_type_0_dataset_velocities = H5Dcreate(part_type_0_group, "Velocities", H5T_NATIVE_DOUBLE, dataspace_3d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (
        part_type_0_dataset_part_ids == H5I_INVALID_HID
        || part_type_0_dataset_masses == H5I_INVALID_HID
        || part_type_0_dataset_coordinates == H5I_INVALID_HID
        || part_type_0_dataset_velocities == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 datasets.");
        goto err_create_hdf5_datasets;
    }
    
    /* Write attributes for header */
    const int num_files_per_snapshot = 1;
    H5Awrite(header_attr_num_files_per_snapshot, H5T_NATIVE_INT, &num_files_per_snapshot);
    H5Awrite(header_attr_num_part_this_file, H5T_NATIVE_INT, &num_particles);
    H5Awrite(header_attr_num_part_total, H5T_NATIVE_INT, &num_particles);
    H5Awrite(header_attr_time, H5T_NATIVE_DOUBLE, &simulation_status->t);

    /* Write data to HDF5 dataset */
    H5Dwrite(part_type_0_dataset_part_ids, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, particle_ids);
    H5Dwrite(part_type_0_dataset_masses, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, m);
    H5Dwrite(part_type_0_dataset_coordinates, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
    H5Dwrite(part_type_0_dataset_velocities, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);

    /* Close HDF5 objects */
    H5Dclose(part_type_0_dataset_part_ids);
    H5Dclose(part_type_0_dataset_masses);
    H5Dclose(part_type_0_dataset_coordinates);
    H5Dclose(part_type_0_dataset_velocities);

    H5Aclose(header_attr_num_files_per_snapshot);
    H5Aclose(header_attr_num_part_this_file);
    H5Aclose(header_attr_num_part_total);
    H5Aclose(header_attr_time);

    H5Sclose(dataspace_1d_1);
    H5Sclose(dataspace_1d_objects_count);
    H5Sclose(dataspace_3d_objects_count);

    H5Gclose(header_group);
    H5Gclose(part_type_0_group);

    H5Fclose(file);

    return make_success_error_status();

    H5Dclose(part_type_0_dataset_part_ids);
    H5Dclose(part_type_0_dataset_masses);
    H5Dclose(part_type_0_dataset_coordinates);
    H5Dclose(part_type_0_dataset_velocities);
err_create_hdf5_datasets:
    H5Aclose(header_attr_num_files_per_snapshot);
    H5Aclose(header_attr_num_part_this_file);
    H5Aclose(header_attr_num_part_total);
    H5Aclose(header_attr_time);
err_create_hdf5_header_attr:
    H5Sclose(dataspace_1d_1);
    H5Sclose(dataspace_1d_objects_count);
    H5Sclose(dataspace_3d_objects_count);
err_create_hdf5_dataspace:
    H5Gclose(header_group);
    H5Gclose(part_type_0_group);
err_create_hdf5_group:
    H5Fclose(file);
err_create_hdf5_file:
err_write_file_path_string:
err_file_path_memory_alloc:
    free(file_path);
err_output_dir_null:
    return error_status;
}

IN_FILE ErrorStatus output_snapshot_cosmology_hdf5(
    OutputParam *__restrict output_param,
    const CosmologicalSystem *__restrict system,
    const IntegratorParam *__restrict integrator_param,
    const AccelerationParam *__restrict acceleration_param,
    const SimulationStatus *__restrict simulation_status,
    const Settings *__restrict settings
)
{
    ErrorStatus error_status;

    (void) integrator_param;
    (void) acceleration_param;
    (void) settings;

    if (!output_param->output_dir)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Output directory path is NULL.");
        goto err_output_dir_null;
    }

    /* Declare variables */
    const int num_particles = system->num_particles;
    const int *__restrict particle_ids = system->particle_ids;
    const double *__restrict x = system->x;
    const double *__restrict v = system->v;
    const double *__restrict m = system->m;

    /* Make file path string */
    const int file_path_length = (
        strlen(output_param->output_dir)
        + snprintf(NULL, 0, "snapshot_%05d.hdf5", output_param->output_count_)
        + 1  // Null terminator
    );
    char *__restrict file_path = malloc(file_path_length * sizeof(char));
    if (!file_path)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for file path string.");
        goto err_file_path_memory_alloc;
    }
    int actual_file_path_length = snprintf(file_path, file_path_length, "%ssnapshot_%05d.hdf5", output_param->output_dir, output_param->output_count_);

    if (actual_file_path_length < 0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Failed to get storing file path string");
        goto err_write_file_path_string;
    }
    else if (actual_file_path_length >= file_path_length)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Storing file path string is truncated.");
        goto err_write_file_path_string;
    }

    /* Create HDF5 file */
    hid_t file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (!file)
    {
        error_status = WRAP_RAISE_ERROR_FMT(
            GRAV_OS_ERROR,
            "Failed to create HDF5 snapshot file: \"%s\".",
            file_path
        );
        goto err_create_hdf5_file;
    }

    /* Create group */
    hid_t header_group = H5Gcreate(file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t part_type_0_group = H5Gcreate(file, "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t units_group = H5Gcreate(file, "/Units", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (header_group == H5I_INVALID_HID || part_type_0_group == H5I_INVALID_HID || units_group == H5I_INVALID_HID)
    {
        H5Gclose(header_group);
        H5Gclose(part_type_0_group);
        H5Gclose(units_group);
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 groups.");
        goto err_create_hdf5_group;
    }

    /* Create dataspaces */
    hsize_t dims_1d_1[1] = {1};
    hsize_t dims_1d_objects_count[1] = {num_particles};
    hsize_t dims_3d_objects_count[2] = {num_particles, 3};
    hid_t dataspace_scaler = H5Screate(H5S_SCALAR);
    hid_t dataspace_1d_1 = H5Screate_simple(1, dims_1d_1, NULL);
    hid_t dataspace_1d_objects_count = H5Screate_simple(1, dims_1d_objects_count, NULL);
    hid_t dataspace_3d_objects_count = H5Screate_simple(2, dims_3d_objects_count, NULL);
    if (
        dataspace_scaler == H5I_INVALID_HID
        || dataspace_1d_1 == H5I_INVALID_HID
        || dataspace_1d_objects_count == H5I_INVALID_HID
        || dataspace_3d_objects_count == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 dataspace.");
        goto err_create_hdf5_dataspace;
    }
    
    /* Create attributes for header */
    hid_t header_attr_num_files_per_snapshot = H5Acreate(header_group, "NumFilesPerSnapshot", H5T_NATIVE_INT, dataspace_1d_1, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_num_part_this_file = H5Acreate(header_group, "NumPart_ThisFile", H5T_NATIVE_INT, dataspace_1d_1, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_num_part_total = H5Acreate(header_group, "NumPart_Total", H5T_NATIVE_INT, dataspace_1d_1, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_time = H5Acreate(header_group, "Time", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_redshift = H5Acreate(header_group, "Redshift", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_omega_m = H5Acreate(header_group, "Omega0", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    hid_t header_attr_omega_lambda = H5Acreate(header_group, "OmegaLambda", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    if (
        header_attr_num_files_per_snapshot == H5I_INVALID_HID
        || header_attr_num_part_this_file == H5I_INVALID_HID
        || header_attr_num_part_total == H5I_INVALID_HID
        || header_attr_time == H5I_INVALID_HID
        || header_attr_redshift == H5I_INVALID_HID
        || header_attr_omega_m == H5I_INVALID_HID
        || header_attr_omega_lambda == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 attribute for header.");
        goto err_create_hdf5_header_attr;
    }

    /* Create datasets for PartType0 */
    hid_t part_type_0_dataset_part_ids = H5Dcreate(part_type_0_group, "ParticleIDs", H5T_NATIVE_INT, dataspace_1d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t part_type_0_dataset_masses;
    hid_t part_type_0_dataset_coordinates;
    hid_t part_type_0_dataset_velocities;
    
    if (output_param->mass_output_dtype == OUTPUT_DTYPE_FLOAT)
    {
        part_type_0_dataset_masses = H5Dcreate(part_type_0_group, "Masses", H5T_NATIVE_FLOAT, dataspace_1d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        part_type_0_dataset_masses = H5Dcreate(part_type_0_group, "Masses", H5T_NATIVE_DOUBLE, dataspace_1d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (output_param->coordinate_output_dtype == OUTPUT_DTYPE_FLOAT)
    {
        part_type_0_dataset_coordinates = H5Dcreate(part_type_0_group, "Coordinates", H5T_NATIVE_FLOAT, dataspace_3d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        part_type_0_dataset_coordinates = H5Dcreate(part_type_0_group, "Coordinates", H5T_NATIVE_DOUBLE, dataspace_3d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (output_param->velocity_output_dtype == OUTPUT_DTYPE_FLOAT)
    {
        part_type_0_dataset_velocities = H5Dcreate(part_type_0_group, "Velocities", H5T_NATIVE_FLOAT, dataspace_3d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        part_type_0_dataset_velocities = H5Dcreate(part_type_0_group, "Velocities", H5T_NATIVE_DOUBLE, dataspace_3d_objects_count, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (
        part_type_0_dataset_part_ids == H5I_INVALID_HID
        || part_type_0_dataset_masses == H5I_INVALID_HID
        || part_type_0_dataset_coordinates == H5I_INVALID_HID
        || part_type_0_dataset_velocities == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 datasets.");
        goto err_create_hdf5_datasets;
    }
    
    /* Create attributes for units */
    hid_t units_attr_current_unit = H5Acreate(units_group, "Unit current in cgs (U_I)", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    hid_t units_attr_length_unit = H5Acreate(units_group, "Unit length in cgs (U_L)", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    hid_t units_attr_mass_unit = H5Acreate(units_group, "Unit mass in cgs (U_M)", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    hid_t units_attr_time_unit = H5Acreate(units_group, "Unit time in cgs (U_t)", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    hid_t units_attr_temperature_unit = H5Acreate(units_group, "Unit temperature in cgs (U_T)", H5T_NATIVE_DOUBLE, dataspace_scaler, H5P_DEFAULT, H5P_DEFAULT);
    if (
        units_attr_current_unit == H5I_INVALID_HID
        || units_attr_length_unit == H5I_INVALID_HID
        || units_attr_mass_unit == H5I_INVALID_HID
        || units_attr_time_unit == H5I_INVALID_HID
        || units_attr_temperature_unit == H5I_INVALID_HID
    )
    {
        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, "Failed to create HDF5 attribute for units.");
        goto err_create_hdf5_units_attr;
    }
    
    /* Write attributes for header */
    const int num_files_per_snapshot = 1;
    H5Awrite(header_attr_num_files_per_snapshot, H5T_NATIVE_INT, &num_files_per_snapshot);
    H5Awrite(header_attr_num_part_this_file, H5T_NATIVE_INT, &num_particles);
    H5Awrite(header_attr_num_part_total, H5T_NATIVE_INT, &num_particles);
    H5Awrite(header_attr_time, H5T_NATIVE_DOUBLE, &(simulation_status->t));
    const double z = 1.0 / simulation_status->t - 1.0;
    H5Awrite(header_attr_redshift, H5T_NATIVE_DOUBLE, &z);
    H5Awrite(header_attr_omega_m, H5T_NATIVE_DOUBLE, &(system->omega_m));
    H5Awrite(header_attr_omega_lambda, H5T_NATIVE_DOUBLE, &(system->omega_lambda));

    /* Write data to HDF5 dataset */
    H5Dwrite(part_type_0_dataset_part_ids, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, particle_ids);
    H5Dwrite(part_type_0_dataset_masses, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, m);
    H5Dwrite(part_type_0_dataset_coordinates, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
    H5Dwrite(part_type_0_dataset_velocities, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);

    /* Write attributes for units */
    double unit_current = 1.0;
    double unit_temperature = 1.0;
    H5Awrite(units_attr_current_unit, H5T_NATIVE_DOUBLE, &unit_current);
    H5Awrite(units_attr_length_unit, H5T_NATIVE_DOUBLE, &(system->unit_length));
    H5Awrite(units_attr_mass_unit, H5T_NATIVE_DOUBLE, &(system->unit_mass));
    H5Awrite(units_attr_time_unit, H5T_NATIVE_DOUBLE, &(system->unit_time));
    H5Awrite(units_attr_temperature_unit, H5T_NATIVE_DOUBLE, &unit_temperature);

    /* Close HDF5 objects */
    H5Dclose(part_type_0_dataset_part_ids);
    H5Dclose(part_type_0_dataset_masses);
    H5Dclose(part_type_0_dataset_coordinates);
    H5Dclose(part_type_0_dataset_velocities);

    H5Aclose(units_attr_current_unit);
    H5Aclose(units_attr_length_unit);
    H5Aclose(units_attr_mass_unit);
    H5Aclose(units_attr_time_unit);
    H5Aclose(units_attr_temperature_unit);

    H5Aclose(header_attr_num_files_per_snapshot);
    H5Aclose(header_attr_num_part_this_file);
    H5Aclose(header_attr_num_part_total);
    H5Aclose(header_attr_time);
    H5Aclose(header_attr_redshift);
    H5Aclose(header_attr_omega_m);
    H5Aclose(header_attr_omega_lambda);

    H5Sclose(dataspace_scaler);
    H5Sclose(dataspace_1d_1);
    H5Sclose(dataspace_1d_objects_count);
    H5Sclose(dataspace_3d_objects_count);

    H5Gclose(header_group);
    H5Gclose(part_type_0_group);
    H5Gclose(units_group);

    H5Fclose(file);

    return make_success_error_status();

    H5Aclose(units_attr_current_unit);
    H5Aclose(units_attr_length_unit);
    H5Aclose(units_attr_mass_unit);
    H5Aclose(units_attr_time_unit);
    H5Aclose(units_attr_temperature_unit);
err_create_hdf5_units_attr:
    H5Dclose(part_type_0_dataset_part_ids);
    H5Dclose(part_type_0_dataset_masses);
    H5Dclose(part_type_0_dataset_coordinates);
    H5Dclose(part_type_0_dataset_velocities);
err_create_hdf5_datasets:
    H5Aclose(header_attr_num_files_per_snapshot);
    H5Aclose(header_attr_num_part_this_file);
    H5Aclose(header_attr_num_part_total);
    H5Aclose(header_attr_time);
    H5Aclose(header_attr_redshift);
    H5Aclose(header_attr_omega_m);
    H5Aclose(header_attr_omega_lambda);
err_create_hdf5_header_attr:
    H5Sclose(dataspace_1d_1);
    H5Sclose(dataspace_1d_objects_count);
    H5Sclose(dataspace_3d_objects_count);
err_create_hdf5_dataspace:
    H5Gclose(header_group);
    H5Gclose(part_type_0_group);
    H5Gclose(units_group);
err_create_hdf5_group:
    H5Fclose(file);
err_create_hdf5_file:
err_write_file_path_string:
err_file_path_memory_alloc:
    free(file_path);
err_output_dir_null:
    return error_status;
}
#endif
