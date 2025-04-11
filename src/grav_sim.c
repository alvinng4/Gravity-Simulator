#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grav_sim.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_HDF5
#include <hdf5.h>
#endif

#ifdef USE_FFTW3
#include <fftw3.h>
#endif

/**
 * \brief Print simulation information.
 * 
 * \param system Pointer to the system.
 * \param integrator_param Pointer to the integrator parameters.
 * \param acceleration_param Pointer to the acceleration parameters.
 * \param output_param Pointer to the output parameters.
 * \param settings Pointer to the settings.
 * \param tf Simulation time.
 */
IN_FILE void print_simulation_info(
    const System *__restrict system,
    const IntegratorParam *__restrict integrator_param,
    const AccelerationParam *__restrict acceleration_param,
    const OutputParam *__restrict output_param,
    const Settings *__restrict settings,
    const double tf
);

#ifdef USE_FFTW3
/**
 * \brief Print cosmological simulation information.
 * 
 * \param system Pointer to the cosmological system.
 * \param integrator_param Pointer to the integrator parameters.
 * \param acceleration_param Pointer to the acceleration parameters.
 * \param output_param Pointer to the output parameters.
 * \param settings Pointer to the settings.
 * \param a_begin Initial scale factor.
 * \param a_final Final scale factor.
 * \param pm_grid_size Size of the PM grid.
 */
IN_FILE void print_cosmological_simulation_info(
    const CosmologicalSystem *__restrict system,
    const IntegratorParam *__restrict integrator_param,
    const AccelerationParam *__restrict acceleration_param,
    const OutputParam *__restrict output_param,
    const Settings *__restrict settings,
    const double a_begin,
    const double a_final,
    const int pm_grid_size
);
#endif

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

    /* Print message */
    if (settings->verbose >= GRAV_VERBOSITY_NORMAL)
    {
        print_compilation_info();
        print_simulation_info(
            system,
            integrator_param,
            acceleration_param,
            output_param,
            settings,
            tf
        );
    }

    if (settings->verbose >= GRAV_VERBOSITY_IGNORE_INFO)
    {
        fputs("Launching simulation...\n", stdout);
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
    if (output_param->method == OUTPUT_METHOD_CSV)
    {
        return WRAP_RAISE_ERROR(
            GRAV_VALUE_ERROR,
            "CSV output is not supported for cosmological simulations. Please use HDF5 output format."
        );
    }
    error_status = WRAP_TRACEBACK(finalize_output_param(output_param, settings));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Check a_begin and a_final */
    if (a_begin <= 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "a_begin must be positive. Got: %g",
            a_begin
        );
    }
    if (a_final <= 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "a_final must be positive. Got: %g",
            a_final
        );
    }
    if (a_final < a_begin)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "a_final must be greater or equal to a_begin. Got: %g, %g",
            a_begin,
            a_final
        );
    }

    /* Check pm_grid_size */
    if (pm_grid_size <= 0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "pm_grid_size must be positive. Got: %d",
            pm_grid_size
        );
    }

    if (settings->verbose >= GRAV_VERBOSITY_NORMAL)
    {
        print_compilation_info();
        print_cosmological_simulation_info(
            system,
            integrator_param,
            acceleration_param,
            output_param,
            settings,
            a_begin,
            a_final,
            pm_grid_size
        );
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

WIN32DLL_API const char* get_grav_sim_logo_string(void)
{
    const char *logo = (
        "                                              __                   \n"
        "    __   _ __    __     __  __           ____/\\_\\    ___ ___       \n"
        "  /'_ `\\/\\`'__\\/'__`\\  /\\ \\/\\ \\         /',__\\/\\ \\ /' __` __`\\     \n"
        " /\\ \\L\\ \\ \\ \\//\\ \\L\\.\\_\\ \\ \\_/ |       /\\__, `\\ \\ \\/\\ \\/\\ \\/\\ \\    \n"
        " \\ \\____ \\ \\_\\\\ \\__/.\\_\\\\ \\___/        \\/\\____/\\ \\_\\ \\_\\ \\_\\ \\_\\   \n"
        "  \\/___L\\ \\/_/ \\/__/\\/_/ \\/__/   _______\\/___/  \\/_/\\/_/\\/_/\\/_/   \n"
        "    /\\____/                     /\\______\\                          \n"
        "    \\_/__/                      \\/______/                          \n"
    );

    return logo;
}

IN_FILE void print_simulation_info(
    const System *__restrict system,
    const IntegratorParam *__restrict integrator_param,
    const AccelerationParam *__restrict acceleration_param,
    const OutputParam *__restrict output_param,
    const Settings *__restrict settings,
    const double tf
)
{
    const char *__restrict new_line = "\n";
    const char *__restrict straight_line = "-----------------------------------------------------------------\n";
    
    fputs("Simulation parameters:\n", stdout);

    /* TF */
    printf("  tf: %g\n", tf);
    fputs(new_line, stdout);

    /* System */
    fputs("System:\n", stdout);
    printf("  Number of particles: %d\n", system->num_particles);
    printf("  Gravitational constant: %g\n", system->G);
    fputs(new_line, stdout);

    /* Integrator */
    fputs("Integrator parameters:\n", stdout);

    // Integrator name
    switch(integrator_param->integrator)
    {
        case INTEGRATOR_EULER:
            fputs("  Integrator: Euler\n", stdout);
            break;
        case INTEGRATOR_EULER_CROMER:
            fputs("  Integrator: Euler-Cromer\n", stdout);
            break;
        case INTEGRATOR_RK4:
            fputs("  Integrator: RK4\n", stdout);
            break;
        case INTEGRATOR_LEAPFROG:
            fputs("  Integrator: Leapfrog\n", stdout);
            break;
        case INTEGRATOR_RKF45:
            fputs("  Integrator: RKF4(5)\n", stdout);
            break;
        case INTEGRATOR_DOPRI:
            fputs("  Integrator: DOPRI\n", stdout);
            break;
        case INTEGRATOR_DVERK:
            fputs("  Integrator: DVERK\n", stdout);
            break;
        case INTEGRATOR_RKF78:
            fputs("  Integrator: RKF7(8)\n", stdout);
            break;
        case INTEGRATOR_IAS15:
            fputs("  Integrator: IAS15\n", stdout);
            break;
        case INTEGRATOR_WHFAST:
            fputs("  Integrator: WHFAST\n", stdout);
            break;
        default:
            fputs("  Integrator: Unknown\n", stdout);
            break;
    }

    // Time step / tolerance
    switch (integrator_param->integrator)
    {
        case INTEGRATOR_EULER:
        case INTEGRATOR_EULER_CROMER:
        case INTEGRATOR_RK4:
        case INTEGRATOR_LEAPFROG:
        case INTEGRATOR_WHFAST:
            printf("  dt: %g\n", integrator_param->dt);
            break;
        case INTEGRATOR_RKF45:
        case INTEGRATOR_DOPRI:
        case INTEGRATOR_DVERK:
        case INTEGRATOR_RKF78:
        case INTEGRATOR_IAS15:
            printf("  Tolerance: %g\n", integrator_param->tolerance);
            printf("  Initial dt (if applicable): %g\n", integrator_param->initial_dt);
            break;
    }

    // WHFast remove invalid particles
    if (integrator_param->integrator == INTEGRATOR_WHFAST)
    {
        printf("  WHFast remove invalid particles: %s\n", integrator_param->whfast_remove_invalid_particles ? "true" : "false");
    }
    fputs(new_line, stdout);
    
    /* Acceleration */
    fputs("Acceleration parameters:\n", stdout);

    // Method
    switch (acceleration_param->method)
    {
        case ACCELERATION_METHOD_PAIRWISE:
            fputs("  Acceleration method: Pairwise\n", stdout);
            break;
        case ACCELERATION_METHOD_MASSLESS:
            fputs("  Acceleration method: Massless\n", stdout);
            break;
        case ACCELERATION_METHOD_BARNES_HUT:
            fputs("  Acceleration method: Barnes-Hut\n", stdout);
            break;
        case ACCELERATION_METHOD_PM:
            fputs("  Acceleration method: Particle-mesh\n", stdout);
            break;
        default:
            fputs("  Acceleration method: Unknown\n", stdout);
            break;
    }

    // Softening length
    printf("  Softening length: %g\n", acceleration_param->softening_length);

    // Opening angle and Max number of particles per leaf
    if (acceleration_param->method == ACCELERATION_METHOD_BARNES_HUT)
    {
        printf("  Opening angle: %g\n", acceleration_param->opening_angle);
        printf("  Max number of particles per leaf: %d\n", acceleration_param->max_num_particles_per_leaf);
    }
    fputs(new_line, stdout);

    /* Output */
    fputs("Output parameters:\n", stdout);

    // Method
    switch (output_param->method)
    {
        case OUTPUT_METHOD_DISABLED:
            fputs("  Output method: Disabled\n", stdout);
            break;
        case OUTPUT_METHOD_CSV:
            fputs("  Output method: CSV\n", stdout);
            break;
        case OUTPUT_METHOD_HDF5:
            fputs("  Output method: HDF5\n", stdout);
            break;
        default:
            fputs("  Output method: Unknown\n", stdout);
            break;
    }

    // Output directory
    printf("  Output directory: %s\n", output_param->output_dir);

    // Output initial
    printf("  Output initial condition: %s\n", output_param->output_initial ? "true" : "false");

    // Output interval
    printf("  Output interval: %g\n", output_param->output_interval);

    // Output data types
    switch (output_param->coordinate_output_dtype)
    {
        case OUTPUT_DTYPE_FLOAT:
            fputs("  Coordinate output data type: float\n", stdout);
            break;
        case OUTPUT_DTYPE_DOUBLE:
            fputs("  Coordinate output data type: double\n", stdout);
            break;
        default:
            fputs("  Coordinate output data type: Unknown\n", stdout);
            break;
    }
    switch (output_param->velocity_output_dtype)
    {
        case OUTPUT_DTYPE_FLOAT:
            fputs("  Velocity output data type: float\n", stdout);
            break;
        case OUTPUT_DTYPE_DOUBLE:
            fputs("  Velocity output data type: double\n", stdout);
            break;
        default:
            fputs("  Velocity output data type: Unknown\n", stdout);
            break;
    }
    switch (output_param->mass_output_dtype)
    {
        case OUTPUT_DTYPE_FLOAT:
            fputs("  Mass output data type: float\n", stdout);
            break;
        case OUTPUT_DTYPE_DOUBLE:
            fputs("  Mass output data type: double\n", stdout);
            break;
        default:
            fputs("  Mass output data type: Unknown\n", stdout);
            break;
    }
    fputs(new_line, stdout);

    /* Settings */
    fputs("Settings:\n", stdout);
    
    // Verbose
    switch (settings->verbose)
    {
        case GRAV_VERBOSITY_IGNORE_ALL:
            fputs("  Verbose: Ignore all\n", stdout);
            break;
        case GRAV_VERBOSITY_IGNORE_INFO:
            fputs("  Verbose: Ignore info\n", stdout);
            break;
        case GRAV_VERBOSITY_NORMAL:
            fputs("  Verbose: Normal\n", stdout);
            break;
        case GRAV_VERBOSITY_VERBOSE:
            fputs("  Verbose: Verbose\n", stdout);
            break;
        default:
            fputs("  Verbose: Unknown\n", stdout);
            break;
    }

    printf("  Enable progress bar: %s\n", settings->enable_progress_bar ? "true" : "false");
    fputs(straight_line, stdout);
}

#ifdef USE_FFTW3
IN_FILE void print_cosmological_simulation_info(
    const CosmologicalSystem *__restrict system,
    const IntegratorParam *__restrict integrator_param,
    const AccelerationParam *__restrict acceleration_param,
    const OutputParam *__restrict output_param,
    const Settings *__restrict settings,
    const double a_begin,
    const double a_final,
    const int pm_grid_size
)
{
    const char *__restrict new_line = "\n";
    const char *__restrict straight_line = "-----------------------------------------------------------------\n";
    
    fputs("Simulation parameters:\n", stdout);

    /* System */
    fputs("System:\n", stdout);
    printf("  Number of particles: %d\n", system->num_particles);
    printf("  Hubble parameter: %g\n", system->h0);
    printf("  Omega_m: %g\n", system->omega_m);
    printf("  Omega_lambda: %g\n", system->omega_lambda);
    printf("  Omega k: %g\n", system->omega_k);
    printf("  Box center: (%g, %g, %g)\n", system->box_center[0], system->box_center[1], system->box_center[2]);
    printf("  Box width: %g\n", system->box_width);
    printf("  Unit mass in cgs: %g\n", system->unit_mass);
    printf("  Unit length in cgs: %g\n", system->unit_length);
    printf("  Unit time in cgs: %g\n", system->unit_time);

    printf("  a_begin: %g\n", a_begin);
    printf("  a_final: %g\n", a_final);

    fputs(new_line, stdout);

    /* Integrator */
    fputs("Integrator parameters:\n", stdout);
    fputs("  Integrator: Leapfrog\n", stdout);
    printf("  dt: %g\n", integrator_param->dt);
    fputs(new_line, stdout);
    
    /* Acceleration */
    fputs("Acceleration parameters:\n", stdout);

    // Method
    switch (acceleration_param->method)
    {
        case ACCELERATION_METHOD_PAIRWISE:
            fputs("  Acceleration method: Pairwise\n", stdout);
            break;
        case ACCELERATION_METHOD_MASSLESS:
            fputs("  Acceleration method: Massless\n", stdout);
            break;
        case ACCELERATION_METHOD_BARNES_HUT:
            fputs("  Acceleration method: Barnes-Hut\n", stdout);
            break;
        case ACCELERATION_METHOD_PM:
            fputs("  Acceleration method: Particle-mesh\n", stdout);
            break;
        default:
            fputs("  Acceleration method: Unknown\n", stdout);
            break;
    }

    printf("  Particle-mesh grid size: %d\n", pm_grid_size);

    // // Softening length
    // printf("  Softening length: %g\n", acceleration_param->softening_length);

    // // Opening angle and Max number of particles per leaf
    // if (acceleration_param->method == ACCELERATION_METHOD_BARNES_HUT)
    // {
    //     printf("  Opening angle: %g\n", acceleration_param->opening_angle);
    //     printf("  Max number of particles per leaf: %d\n", acceleration_param->max_num_particles_per_leaf);
    // }
    fputs(new_line, stdout);

    /* Output */
    fputs("Output parameters:\n", stdout);

    // Method
    switch (output_param->method)
    {
        case OUTPUT_METHOD_CSV:
            fputs("  Output method: CSV\n", stdout);
            break;
        case OUTPUT_METHOD_HDF5:
            fputs("  Output method: HDF5\n", stdout);
            break;
        default:
            fputs("  Output method: Unknown\n", stdout);
            break;
    }

    // Output directory
    printf("  Output directory: %s\n", output_param->output_dir);

    // Output initial
    printf("  Output initial condition: %s\n", output_param->output_initial ? "true" : "false");

    // Output interval
    printf("  Output interval: %g\n", output_param->output_interval);

    // Output data types
    switch (output_param->coordinate_output_dtype)
    {
        case OUTPUT_DTYPE_FLOAT:
            fputs("  Coordinate output data type: float\n", stdout);
            break;
        case OUTPUT_DTYPE_DOUBLE:
            fputs("  Coordinate output data type: double\n", stdout);
            break;
        default:
            fputs("  Coordinate output data type: Unknown\n", stdout);
            break;
    }
    switch (output_param->velocity_output_dtype)
    {
        case OUTPUT_DTYPE_FLOAT:
            fputs("  Velocity output data type: float\n", stdout);
            break;
        case OUTPUT_DTYPE_DOUBLE:
            fputs("  Velocity output data type: double\n", stdout);
            break;
        default:
            fputs("  Velocity output data type: Unknown\n", stdout);
            break;
    }
    switch (output_param->mass_output_dtype)
    {
        case OUTPUT_DTYPE_FLOAT:
            fputs("  Mass output data type: float\n", stdout);
            break;
        case OUTPUT_DTYPE_DOUBLE:
            fputs("  Mass output data type: double\n", stdout);
            break;
        default:
            fputs("  Mass output data type: Unknown\n", stdout);
            break;
    }
    fputs(new_line, stdout);

    /* Settings */
    fputs("Settings:\n", stdout);
    
    // Verbose
    switch (settings->verbose)
    {
        case GRAV_VERBOSITY_IGNORE_ALL:
            fputs("  Verbose: Ignore all\n", stdout);
            break;
        case GRAV_VERBOSITY_IGNORE_INFO:
            fputs("  Verbose: Ignore info\n", stdout);
            break;
        case GRAV_VERBOSITY_NORMAL:
            fputs("  Verbose: Normal\n", stdout);
            break;
        case GRAV_VERBOSITY_VERBOSE:
            fputs("  Verbose: Verbose\n", stdout);
            break;
        default:
            fputs("  Verbose: Unknown\n", stdout);
            break;
    }

    printf("  Enable progress bar: %s\n", settings->enable_progress_bar ? "true" : "false");
    fputs(straight_line, stdout);
}
#endif

WIN32DLL_API void print_compilation_info(void)
{
    const char *new_line = "\n";
    const char *straight_line = "-----------------------------------------------------------------\n";

    fputs(straight_line, stdout);
    fputs(get_grav_sim_logo_string(), stdout);
    fputs(new_line, stdout);
    fputs(new_line, stdout);

    printf("grav_sim version %s\n", VERSION_INFO);
    fputs(new_line, stdout);

    /* Print OS information */
#ifdef _WIN32
    fputs("Operating System: Windows\n", stdout);
#elif __APPLE__
    fputs("Operating System: MacOS\n", stdout);
#elif __linux__
    fputs("Operating System: Linux\n", stdout);
#else
    fputs("Operating System: Unknown\n", stdout);
#endif

    /* Print compilation information */
    fputs("Compilation Info:\n", stdout);

    /* OpenMP */
#ifdef USE_OPENMP
    fputs("  Compiled with OpenMP: true\n", stdout);
#else
    fputs("  Compiled with OpenMP: false\n", stdout);
#endif

    /* HDF5 */
#ifdef USE_HDF5
    fputs("  Compiled with HDF5: true\n", stdout);
    #ifdef H5_VERS_MAJOR
        printf("    Version: %d.%d.%d\n", H5_VERS_MAJOR, H5_VERS_MINOR, H5_VERS_RELEASE);
    #endif
#else
    fputs("  Compiled with HDF5: false\n", stdout);
#endif

    /* FFTW3 */
#ifdef USE_FFTW3
    fputs("  Compiled with FFTW3: true\n", stdout);
    printf("    Version: %s\n", fftw_version);;
#else
    fputs("  Compiled with FFTW3: false\n", stdout);
#endif

    fputs(new_line, stdout);

    /* Build and compiler info */
    printf("Build time: %s %s\n", __DATE__, __TIME__);

#ifdef _MSC_VER
    printf("Compiler: MSVC (version: %d)\n", _MSC_VER);
#elif defined(__clang__)
    printf("Compiler: Clang (version: %d)\n", __clang_major__);
#elif defined(__GNUC__)
    printf("Compiler: GCC (version: %d)\n", __GNUC__);
#endif
    fputs(straight_line, stdout);
}
