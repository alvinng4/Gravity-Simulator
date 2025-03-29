/**
 * \file output.c
 * \brief Function definitions for simulation output
 * 
 * \author Ching-Yin Ng
 * \date March 2025
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
            const int error_msg_len = (
                strlen("Unknown output method. Got: ")
                + snprintf(NULL, 0, "%d", output_method)
                + 1  // Null terminator
            );
            char *error_msg = malloc(error_msg_len * sizeof(char));
            if (!error_msg)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR, "Unknown output method and failed to allocate memory for error message"
                );
            }

            const int actual_error_msg_len = snprintf(
                error_msg,
                error_msg_len,
                "Unknown output method. Got: %d",
                output_method
            );

            if (actual_error_msg_len < 0)
            {
                free(error_msg);
                return WRAP_RAISE_ERROR(GRAV_UNKNOWN_ERROR, "Unknown output method and failed to generate error message");
            }
            else if (actual_error_msg_len >= error_msg_len)
            {
                free(error_msg);
                return WRAP_RAISE_ERROR(GRAV_UNKNOWN_ERROR, "Unknown output method and error message are truncated");
            }

            ErrorStatus error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, error_msg);
            free(error_msg);
            return error_status;
        }
    }

    return make_success_error_status();
}

ErrorStatus finalize_output_param(OutputParam *__restrict output_param)
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
        const int error_msg_len = (
            strlen("Output interval must be positive. Got: ")
            + snprintf(NULL, 0, "%.17g", output_param->output_interval)
            + 1  // Null terminator
        );
        char *error_msg = malloc(error_msg_len * sizeof(char));
        if (!error_msg)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR, "Output interval is negative and failed to allocate memory for error message"
            );
        }

        const int actual_error_msg_len = snprintf(
            error_msg,
            error_msg_len,
            "Output interval must be positive. Got: %.17g",
            output_param->output_interval
        );

        if (actual_error_msg_len < 0)
        {
            free(error_msg);
            return WRAP_RAISE_ERROR(GRAV_UNKNOWN_ERROR, "Output interval is negative and failed to generate error message");
        }
        else if (actual_error_msg_len >= error_msg_len)
        {
            free(error_msg);
            return WRAP_RAISE_ERROR(GRAV_UNKNOWN_ERROR, "Output interval is negative and error message are truncated");
        }

        ErrorStatus error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, error_msg);
        free(error_msg);
        return error_status;
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
            size_t buffer_size = (
                strlen("Directory path for storing snapshots must end with a trailing slash (\"/\"). Got: \"\".")
                + strlen(output_param->output_dir)
                + 1  // Null terminator
            );
            char *error_message = malloc(buffer_size * sizeof(char));
            if (!error_message)
            {
                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for error message.");
            }
            snprintf(
                error_message,
                buffer_size,
                "Directory path for storing snapshots must end with a trailing slash (\"/\"). Got: \"%s\".",
                output_param->output_dir
            );
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, error_message);
        }
    }

    /* Create directory */
#ifdef _WIN32
    if (_mkdir(output_param->output_dir) == -1)
    {
        if (GetFileAttributes(output_param->output_dir) == INVALID_FILE_ATTRIBUTES)
        {
            size_t buffer_size = (
                strlen("Failed to access path for storing snapshots: \"\".")
                + strlen(output_param->output_dir)
                + 1  // Null terminator
            );
            char *error_message = malloc(buffer_size * sizeof(char));
            if (!error_message)
            {
                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for error message.");
            }
            snprintf(
                error_message,
                buffer_size,
                "Failed to access path for storing snapshots: \"%s\".",
                output_param->output_dir
            );
            return WRAP_RAISE_ERROR(GRAV_OS_ERROR, error_message);
        }
        else if (GetFileAttributes(output_param->output_dir) & FILE_ATTRIBUTE_DIRECTORY)
        {
            int buffer_size = (
                strlen("Directory for storing snapshots already exists. The files will be overwritten. Directory: \"\".")
                + strlen(output_param->output_dir)
                + 1  // Null terminator
            );
            char *warning_message = malloc(buffer_size * sizeof(char));
            if (!warning_message)
            {
                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for warning message.");
            }
            const int actual_warning_message_length = snprintf(
                warning_message,
                buffer_size,
                "Directory for storing snapshots already exists. The files will be overwritten. Directory: \"%s\".",
                output_param->output_dir
            );
            if (actual_warning_message_length < 0)
            {
                return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Failed to get warning message string.");
            }
            else if (actual_warning_message_length >= buffer_size)
            {
                return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Warning message string is truncated.");
            }
            WRAP_RAISE_WARNING(warning_message);
        }
    }
#else
    struct stat st = {0};
    if (mkdir(output_param->output_dir, 0777) == -1)
    {
        if(stat(output_param->output_dir, &st) != 0)
        {
            size_t buffer_size = (
                strlen("Failed to access path for storing snapshots: \"\".")
                + strlen(output_param->output_dir)
                + 1  // Null terminator
            );
            char *error_message = malloc(buffer_size * sizeof(char));
            if (!error_message)
            {
                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for error message.");
            }
            snprintf(
                error_message,
                buffer_size,
                "Failed to access path for storing snapshots: \"%s\".",
                output_param->output_dir
            );
            return WRAP_RAISE_ERROR(GRAV_OS_ERROR, error_message);
        }

        else if (st.st_mode & S_IFDIR)
        {
            int buffer_size = (
                strlen("Directory for storing snapshots already exists. The files will be overwritten. Directory: \"\".")
                + strlen(output_param->output_dir)
                + 1  // Null terminator
            );
            char *warning_message = malloc(buffer_size * sizeof(char));
            if (!warning_message)
            {
                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for warning message.");
            }
            const int actual_warning_message_length = snprintf(
                warning_message,
                buffer_size,
                "Directory for storing snapshots already exists. The files will be overwritten. Directory: \"%s\".",
                output_param->output_dir
            );
            if (actual_warning_message_length < 0)
            {
                return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Failed to get warning message string.");
            }
            else if (actual_warning_message_length >= buffer_size)
            {
                return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Warning message string is truncated.");
            }
            WRAP_RAISE_WARNING(warning_message);
        }
    }
#endif

    return make_success_error_status();
}

ErrorStatus output_snapshot(
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

IN_FILE ErrorStatus output_snapshot_csv(
    OutputParam *output_param,
    const System *system,
    const IntegratorParam *integrator_param,
    const AccelerationParam *acceleration_param,
    const SimulationStatus *simulation_status,
    const Settings *settings
)
{
    ErrorStatus error_status;

    (void) integrator_param;
    (void) acceleration_param;
    (void) simulation_status;
    (void) settings;

    if (!output_param->output_dir)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Output directory path is NULL.");
        goto err_output_dir_null;
    }

    /* Make file path string */
    const int file_path_length = (
        strlen(output_param->output_dir)
        + snprintf(NULL, 0, "snapshot_%d.csv", output_param->output_count_)
        + 1  // Null terminator
    );
    char *__restrict file_path = malloc(file_path_length * sizeof(char));
    if (!file_path)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for file path string.");
        goto err_file_path_memory_alloc;
    }
    int actual_file_path_length = snprintf(file_path, file_path_length, "%ssnapshot_%d.csv", output_param->output_dir, output_param->output_count_);

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
        const int error_msg_len = (
            strlen("Failed to open file for storing snapshots: \"\".")
            + strlen(file_path)
            + 1  // Null terminator
        );
        char *error_msg = malloc(error_msg_len * sizeof(char));
        if (!error_msg)
        {
            error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for error message.");
            goto err_open_file;
        }

        int actual_error_msg_len = snprintf(
            error_msg,
            error_msg_len,
            "Failed to open file for storing snapshots: \"%s\".",
            file_path
        );

        if (actual_error_msg_len < 0)
        {
            error_status = WRAP_RAISE_ERROR(GRAV_UNKNOWN_ERROR, "Failed to open file for storing snapshots and failed to generate error message.");
            goto err_open_file;
        }
        else if (actual_error_msg_len >= error_msg_len)
        {
            error_status = WRAP_RAISE_ERROR(GRAV_UNKNOWN_ERROR, "Failed to open file for storing snapshots and error message is truncated.");
            goto err_open_file;
        }

        error_status = WRAP_RAISE_ERROR(GRAV_OS_ERROR, error_msg);
        free(error_msg);
        goto err_open_file;
    }

    /* Write header */
    fputs("particle_id,m,x,y,z,vx,vy,vz\n", file);

    /* Write data */
    const int objects_count = system->objects_count;
    const int *__restrict particle_ids = system->particle_ids;
    const double *__restrict x = system->x;
    const double *__restrict v = system->v;
    const double *__restrict m = system->m;

    for (int i = 0; i < objects_count; i++)
    {
        // particle_id
        fprintf(file, "%d,", particle_ids[i]);

        // masses
        if (output_param->mass_output_dtype == OUTPUT_DTYPE_FLOAT)
        {
            fprintf(file, "%.8g,", (float) m[i]);
        }
        else
        {
            fprintf(file, "%.17g,", m[i]);
        }

        // coordinates
        if (output_param->coordinate_output_dtype == OUTPUT_DTYPE_FLOAT)
        {
            fprintf(file, "%.8g,", (float) x[3 * i]);
            fprintf(file, "%.8g,", (float) x[3 * i + 1]);
            fprintf(file, "%.8g,", (float) x[3 * i + 2]);
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
            fprintf(file, "%.8g,", (float) v[3 * i]);
            fprintf(file, "%.8g,", (float) v[3 * i + 1]);
            fprintf(file, "%.8g\n", (float) v[3 * i + 2]);
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

