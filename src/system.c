/**
 * \file system.c
 * \brief System module for grav_sim
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "cosmology.h"
#include "error.h"
#include "math_functions.h"
#include "system.h"
#include "settings.h"


typedef struct
{
    int index;
    double distance;
} HelperSystemSortByDistanceStruct;


WIN32DLL_API System get_new_system(void)
{
    System system;
    system.num_particles = 0;
    system.particle_ids = NULL;
    system.x = NULL;
    system.v = NULL;
    system.m = NULL;
    system.G = 0.000295912208284119496676630; // Default value for AU^3 day^-2 M_sun^-1
    return system;
}

WIN32DLL_API ErrorStatus get_initialized_system(
    System *__restrict system,
    const int num_particles
)
{
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    *system = get_new_system();
    system->num_particles = num_particles;
    system->particle_ids = malloc(num_particles * sizeof(int));
    system->x = calloc(num_particles * 3, sizeof(double));
    system->v = calloc(num_particles * 3, sizeof(double));
    system->m = calloc(num_particles, sizeof(double));
    system->G = 0.000295912208284119496676630; // Default value for AU^3 d^-2

    if (!system->particle_ids || !system->x || !system->v || !system->m)
    {
        free_system(system);
        return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for system");
    }

    for (int i = 0; i < num_particles; i++)
    {
        (system->particle_ids)[i] = i;
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus finalize_system(System *__restrict system)
{
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    if (system->num_particles <= 0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Number of particles must be positive. Got: %d",
            system->num_particles
        );
    }
    if (!system->particle_ids)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System array particle_ids is NULL");
    }
    if (!system->x)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System array x is NULL");
    }
    if (!system->v)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System array v is NULL");
    }
    if (!system->m)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System array m is NULL");
    }
    if (system->G <= 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Gravitational constant G must be positive. Got: %g",
            system->G
        );
    }

    return make_success_error_status();
}

WIN32DLL_API void free_system(System *__restrict system)
{
    free(system->particle_ids);
    free(system->x);
    free(system->v);
    free(system->m);
}

WIN32DLL_API CosmologicalSystem get_new_cosmological_system(void)
{
    CosmologicalSystem system;
    system.num_particles = 0;
    system.particle_ids = NULL;
    system.x = NULL;
    system.v = NULL;
    system.m = NULL;
    system.G = -1.0;
    system.h0 = -1.0;
    system.omega_m = -1.0;
    system.omega_lambda = -1.0;
    system.box_center[0] = 0.0;
    system.box_center[1] = 0.0;
    system.box_center[2] = 0.0;
    system.box_width = -1.0;
    system.unit_mass = 1.0;
    system.unit_length = 1.0;
    system.unit_time = 1.0;

    return system;
}

WIN32DLL_API ErrorStatus get_initialized_cosmological_system(
    CosmologicalSystem *__restrict system,
    const int num_particles
)
{
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    *system = get_new_cosmological_system();
    system->num_particles = num_particles;
    system->particle_ids = malloc(num_particles * sizeof(int));
    system->x = calloc(num_particles * 3, sizeof(double));
    system->v = calloc(num_particles * 3, sizeof(double));
    system->m = calloc(num_particles, sizeof(double));

    if (!system->particle_ids || !system->x || !system->v || !system->m)
    {
        free_cosmological_system(system);
        return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for system");
    }

    for (int i = 0; i < num_particles; i++)
    {
        (system->particle_ids)[i] = i;
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus finalize_cosmological_system(CosmologicalSystem *__restrict system)
{
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    if (system->num_particles <= 0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Number of particles must be positive. Got: %d",
            system->num_particles
        );
    }
    if (!system->particle_ids)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System array particle_ids is NULL");
    }
    if (!system->x)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System array x is NULL");
    }
    if (!system->v)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System array v is NULL");
    }
    if (!system->m)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System array m is NULL");
    }
    if (system->h0 <= 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Hubble constant h0 must be positive. Got: %g",
            system->h0
        );
    }

    if (system->omega_m <= 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "omega_m must be positive. Got: %g",
            system->omega_m
        );
    }

    if (system->omega_lambda <= 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "omega_lambda must be positive. Got: %g",
            system->omega_lambda
        );
    }

    system->G = compute_G(system->omega_m, system->h0);

    if (system->box_width <= 0.0)
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "Box width must be positive. Got: %g",
            system->box_width
        );
    }

    return make_success_error_status();
}

WIN32DLL_API void free_cosmological_system(CosmologicalSystem *__restrict system)
{
    free(system->particle_ids);
    free(system->x);
    free(system->v);
    free(system->m);
}

WIN32DLL_API ErrorStatus check_invalid_idx_double(
    bool *__restrict has_invalid_idx,
    int **invalid_idx_array,
    const double *__restrict array,
    const int arr_size
)
{
    if (!array)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Array is NULL");
    }
    if (!has_invalid_idx)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "has_invalid_idx is NULL");
    }

    int invalid_count = 0;
    int buffer_size = 10;
    int *__restrict invalid_particle_idx = malloc(buffer_size * sizeof(int));
    if (!invalid_particle_idx)
    {
        return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for invalid particle index");
    }

    for (int i = 0; i < arr_size; i++)
    {
        if (!isfinite(array[i]))
        {
            invalid_particle_idx[invalid_count] = i;
            invalid_count++;
        }

        if (invalid_count >= buffer_size)
        {
            buffer_size *= 2;
            int *__restrict new_invalid_particle_idx = realloc(invalid_particle_idx, buffer_size * sizeof(int));
            if (!new_invalid_particle_idx)
            {
                free(invalid_particle_idx);
                return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to reallocate memory for invalid particle index");
            }
            invalid_particle_idx = new_invalid_particle_idx;
        }
    }

    if (invalid_count == 0)
    {
        free(invalid_particle_idx);
        *has_invalid_idx = false;
        return make_success_error_status();
    }
    else
    {
        *has_invalid_idx = true;
        *invalid_idx_array = invalid_particle_idx;
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus check_and_remove_invalid_particles(
    System *__restrict system,
    const Settings *__restrict settings
)
{
    ErrorStatus error_status;
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    /* Declare variables */
    const int num_particles = system->num_particles;
    double *__restrict x = system->x;
    double *__restrict v = system->v;
    double *__restrict m = system->m;
    int *__restrict particle_ids = system->particle_ids;

    if (!x || !v || !m || !particle_ids)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System members are NULL");
    }

    int invalid_count = 0;
    int buffer_size = 10;
    int *__restrict invalid_particle_idx = malloc(buffer_size * sizeof(int));
    if (!invalid_particle_idx)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to allocate memory for invalid particle index");
        goto err_memory;
    }

    for (int i = 0; i < num_particles; i++)
    {
        if (
            !isfinite(x[i * 3 + 0]) ||
            !isfinite(x[i * 3 + 1]) ||
            !isfinite(x[i * 3 + 2]) ||
            !isfinite(v[i * 3 + 0]) ||
            !isfinite(v[i * 3 + 1]) ||
            !isfinite(v[i * 3 + 2]) ||
            !isfinite(m[i])            
        )
        {
            invalid_particle_idx[invalid_count] = i;
            invalid_count++;
        }

        if (invalid_count >= buffer_size)
        {
            buffer_size *= 2;
            int *__restrict new_invalid_particle_idx = realloc(invalid_particle_idx, buffer_size * sizeof(int));
            if (!new_invalid_particle_idx)
            {
                error_status = WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to reallocate memory for invalid particle index");
                goto err_memory;
            }
            invalid_particle_idx = new_invalid_particle_idx;
        }
    }

    if (invalid_count != 0)
    {
        error_status = WRAP_TRACEBACK(remove_invalid_particles(
            system,
            invalid_particle_idx,
            invalid_count,
            settings
        ));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_remove_particles;
        }
    }

    free(invalid_particle_idx);

    return make_success_error_status();

err_remove_particles:
err_memory:
    free(invalid_particle_idx);
    return error_status;
}

WIN32DLL_API ErrorStatus remove_invalid_particles(
    System *__restrict system,
    const int *__restrict remove_idx_list,
    const int num_to_remove,
    const Settings *__restrict settings
)
{
    if (num_to_remove == 0)
    {
        return make_success_error_status();
    }

    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }
    if (!remove_idx_list)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Remove index list is NULL");
    }
    if (!settings)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Settings is NULL");
    }

    if (num_to_remove < 0)
    {
        return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Number of particles to remove must be positive");
    }

    if (settings->verbose >= GRAV_VERBOSITY_VERBOSE)
    {
        fprintf(stderr, "remove_invalid_particles: Removing %d invalid particles. Particle IDs: [%d", num_to_remove, remove_idx_list[0]);
        for (int i = 0; i < num_to_remove; i++)
        {
            fprintf(stderr, ", %d", remove_idx_list[i]);
        }
        fputs("]\n", stderr);
    }

    return WRAP_TRACEBACK(remove_particles(
        system,
        remove_idx_list,
        num_to_remove
    ));
}

WIN32DLL_API ErrorStatus remove_particles(
    System *__restrict system,
    const int *__restrict remove_idx_list,
    const int num_to_remove
)
{
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    const int num_particles = system->num_particles;
    double *__restrict x = system->x;
    double *__restrict v = system->v;
    double *__restrict m = system->m;
    int *__restrict particle_ids = system->particle_ids;

    if (!x || !v || !m || !particle_ids)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System members are NULL");
    }

    for (int i = 0, last_shifted_idx = remove_idx_list[i]; i < num_to_remove; i++)
    {
        const int idx = remove_idx_list[i];
        int idx_next;
        if (i == num_to_remove - 1)
        {
            idx_next = num_particles;
        }
        else
        {
            idx_next = remove_idx_list[i + 1];
        }

        for (int j = 0; j < idx_next - idx - 1; j++)    
        { 
            const int from_idx = last_shifted_idx + i + 1;
            particle_ids[last_shifted_idx] = particle_ids[from_idx];
            x[last_shifted_idx * 3 + 0] = x[from_idx * 3 + 0];
            x[last_shifted_idx * 3 + 1] = x[from_idx * 3 + 1];
            x[last_shifted_idx * 3 + 2] = x[from_idx * 3 + 2];
            v[last_shifted_idx * 3 + 0] = v[from_idx * 3 + 0];
            v[last_shifted_idx * 3 + 1] = v[from_idx * 3 + 1];
            v[last_shifted_idx * 3 + 2] = v[from_idx * 3 + 2];
            m[last_shifted_idx] = m[from_idx];

            last_shifted_idx++;
        }
    }

    system->num_particles -= num_to_remove;

    /* Realloc */
    // size_t new_size = (size_t) system->num_particles;
    // int *new_particle_ids = realloc(system->particle_ids, new_size * sizeof(int));
    // if (!new_particle_ids)
    // {
    //     return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to reallocate memory for particle ids");
    // }
    // system->particle_ids = new_particle_ids;

    // double *new_x = realloc(system->x, new_size * 3 * sizeof(double));
    // if (!new_x)
    // {
    //     return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to reallocate memory for x");
    // }
    // system->x = new_x;

    // double *new_v = realloc(system->v, new_size * 3 * sizeof(double));
    // if (!new_v)
    // {
    //     return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to reallocate memory for v");
    // }
    // system->v = new_v;

    // double *new_m = realloc(system->m, new_size * sizeof(double));
    // if (!new_m)
    // {
    //     return WRAP_RAISE_ERROR(GRAV_MEMORY_ERROR, "Failed to reallocate memory for m");
    // }
    // system->m = new_m;

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus remove_particle_from_double_arr(
    double *__restrict arr,
    const int *__restrict remove_idx_list,
    const int num_to_remove,
    const int dim,
    const int original_size
)
{
    if (!arr)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Array is NULL");
    }
    if (!remove_idx_list)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Remove index list is NULL");
    }
    if (num_to_remove <= 0)
    {
        return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Number of particles to remove must be positive");
    }

    for (int i = 0, last_shifted_idx = remove_idx_list[i]; i < num_to_remove; i++)
    {
        const int idx = remove_idx_list[i];
        int idx_next;
        if (i == num_to_remove - 1)
        {
            idx_next = original_size;
        }
        else
        {
            idx_next = remove_idx_list[i + 1];
        }

        for (int j = 0; j < idx_next - idx - 1; j++)    
        { 
            const int from_idx = last_shifted_idx + i + 1;
            for (int d = 0; d < dim; d++)
            {
                arr[last_shifted_idx * dim + d] = arr[from_idx * dim + d];
            }

            last_shifted_idx++;
        }
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus initialize_built_in_system(
    System *__restrict system,
    const char *__restrict system_name,
    const bool is_memory_initialized
)
{
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }
    if (!system_name)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System name is NULL");
    }

    /* Pre-defined constants */

    // Conversion factor from km^3 s^-2 to AU^3 d^-2
    // double CONVERSION_FACTOR = ((double) 86400.0L * 86400.0L) / (149597870.7L * 149597870.7L * 149597870.7L);

    // GM values (km^3 s^-2)
    // ref: https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf
    const double GM_SI_SUN = 132712440041.279419L;
    const double GM_SI_MERCURY = 22031.868551L;
    const double GM_SI_VENUS = 324858.592000L;
    const double GM_SI_EARTH = 398600.435507L;
    const double GM_SI_MARS = 42828.375816L;
    const double GM_SI_JUPITER = 126712764.100000L;
    const double GM_SI_SATURN = 37940584.841800L;
    const double GM_SI_URANUS = 5794556.400000L;
    const double GM_SI_NEPTUNE = 6836527.100580L;
    const double GM_SI_MOON = 4902.800118L;
    const double GM_SI_PLUTO = 975.500000L;
    const double GM_SI_CERES = 62.62890L;
    const double GM_SI_VESTA = 17.288245L;

    // GM values (AU^3 d^-2)
    // const double GM_SUN = 132712440041.279419L * CONVERSION_FACTOR;
    // const double GM_MERCURY = 22031.868551L * CONVERSION_FACTOR;
    // const double GM_VENUS = 324858.592000L * CONVERSION_FACTOR;
    // const double GM_EARTH = 398600.435507L * CONVERSION_FACTOR;
    // const double GM_MARS = 42828.375816L * CONVERSION_FACTOR;
    // const double GM_JUPITER = 126712764.100000L * CONVERSION_FACTOR;
    // const double GM_SATURN = 37940584.841800L * CONVERSION_FACTOR;
    // const double GM_URANUS = 5794556.400000L * CONVERSION_FACTOR;
    // const double GM_NEPTUNE = 6836527.100580L * CONVERSION_FACTOR;
    // const double GM_MOON = 4902.800118L * CONVERSION_FACTOR;
    // const double GM_PLUTO = 975.500000L * CONVERSION_FACTOR;
    // const double GM_CERES = 62.62890L * CONVERSION_FACTOR;
    // const double GM_VESTA = 17.288245L * CONVERSION_FACTOR;

    // Solar system masses (M_sun^-1)
    const double MASS_SUN = 1.0;
    const double MASS_MERCURY = GM_SI_MERCURY / GM_SI_SUN;
    const double MASS_VENUS = GM_SI_VENUS / GM_SI_SUN;
    const double MASS_EARTH = GM_SI_EARTH / GM_SI_SUN;
    const double MASS_MARS = GM_SI_MARS / GM_SI_SUN;
    const double MASS_JUPITER = GM_SI_JUPITER / GM_SI_SUN;
    const double MASS_SATURN = GM_SI_SATURN / GM_SI_SUN;
    const double MASS_URANUS = GM_SI_URANUS / GM_SI_SUN;
    const double MASS_NEPTUNE = GM_SI_NEPTUNE / GM_SI_SUN;
    const double MASS_MOON = GM_SI_MOON / GM_SI_SUN;
    const double MASS_PLUTO = GM_SI_PLUTO / GM_SI_SUN;
    const double MASS_CERES = GM_SI_CERES / GM_SI_SUN;
    const double MASS_VESTA = GM_SI_VESTA / GM_SI_SUN;

    // Gravitational constant (kg^-1 m^3 s^-2):
    // const double G_SI = 6.67430e-11;
    // Gravitational constant (M_sun^-1 AU^3 d^-2):
    // const double G_default = GM_SUN;

    /*
    * Solar system position and velocities data
    * Units: AU-D
    * Coordinate center: Solar System Barycenter
    * Data dated on A.D. 2024-Jan-01 00:00:00.0000 TDB
    * Computational data generated by NASA JPL Horizons System https://ssd.jpl.nasa.gov/horizons/
    */
    const double POS_SUN[3] = {-7.967955691533730e-03L, -2.906227441573178e-03L, 2.103054301547123e-04L};
    const double POS_MERCURY[3] = {-2.825983269538632e-01L, 1.974559795958082e-01L, 4.177433558063677e-02L};
    const double POS_VENUS[3] = {-7.232103701666379e-01L, -7.948302026312400e-02L, 4.042871428174315e-02L};
    const double POS_EARTH[3] = {-1.738192017257054e-01L, 9.663245550235138e-01L, 1.553901854897183e-04L};
    const double POS_MARS[3] = {-3.013262392582653e-01L, -1.454029331393295e00L, -2.300531433991428e-02L};
    const double POS_JUPITER[3] = {3.485202469657674e00L, 3.552136904413157e00L, -9.271035442798399e-02L};
    const double POS_SATURN[3] = {8.988104223143450e00L, -3.719064854634689e00L, -2.931937777323593e-01L};
    const double POS_URANUS[3] = {1.226302417897505e01L, 1.529738792480545e01L, -1.020549026883563e-01L};
    const double POS_NEPTUNE[3] = {2.983501460984741e01L, -1.793812957956852e00L, -6.506401132254588e-01L};
    const double POS_MOON[3] = {-1.762788124769829e-01L, 9.674377513177153e-01L, 3.236901585768862e-04L};
    const double POS_PLUTO[3] = {1.720200478843485e01L, -3.034155683573043e01L, -1.729127607100611e00L};
    const double POS_CERES[3] = {-1.103880510367569e00L, -2.533340440444230e00L, 1.220283937721780e-01L};
    const double POS_VESTA[3] = {-8.092549658731499e-02L, 2.558381434460076e00L, -6.695836142398572e-02L};

    const double VEL_SUN[3] = {4.875094764261564e-06L, -7.057133213976680e-06L, -4.573453713094512e-08L};
    const double VEL_MERCURY[3] = {-2.232165900189702e-02L, -2.157207103176252e-02L, 2.855193410495743e-04L};
    const double VEL_VENUS[3] = {2.034068201002341e-03L, -2.020828626592994e-02L, -3.945639843855159e-04L};
    const double VEL_EARTH[3] = {-1.723001232538228e-02L, -2.967721342618870e-03L, 6.382125383116755e-07L};
    const double VEL_MARS[3] = {1.424832259345280e-02L, -1.579236181580905e-03L, -3.823722796161561e-04L};
    const double VEL_JUPITER[3] = {-5.470970658852281e-03L, 5.642487338479145e-03L, 9.896190602066252e-05L};
    const double VEL_SATURN[3] = {1.822013845554067e-03L, 5.143470425888054e-03L, -1.617235904887937e-04L};
    const double VEL_URANUS[3] = {-3.097615358317413e-03L, 2.276781932345769e-03L, 4.860433222241686e-05L};
    const double VEL_NEPTUNE[3] = {1.676536611817232e-04L, 3.152098732861913e-03L, -6.877501095688201e-05L};
    const double VEL_MOON[3] = {-1.746667306153906e-02L, -3.473438277358121e-03L, -3.359028758606074e-05L};
    const double VEL_PLUTO[3] = {2.802810313667557e-03L, 8.492056438614633e-04L, -9.060790113327894e-04L};
    const double VEL_CERES[3] = {8.978653480111301e-03L, -4.873256528198994e-03L, -1.807162046049230e-03L};
    const double VEL_VESTA[3] = {-1.017876585480054e-02L, -5.452367109338154e-04L, 1.255870551153315e-03L};

    // Pre-defined systems
    if (strcmp(system_name, "circular_binary_orbit") == 0) 
    {
        const int system_num_particles = 2;
        if (!is_memory_initialized)
        {
            ErrorStatus error_status = WRAP_TRACEBACK(get_initialized_system(system, system_num_particles));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
        else if (system->num_particles < system_num_particles)
        {
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initialized system is not big enough for the built-in system");
        }

        (system->x)[0] = 1.0;
        (system->x)[1] = 0.0;
        (system->x)[2] = 0.0;

        (system->x)[3] = -1.0;
        (system->x)[4] = 0.0;
        (system->x)[5] = 0.0; 

        (system->v)[0] = 0.0;
        (system->v)[1] = 0.5;
        (system->v)[2] = 0.0;

        (system->v)[3] = 0.0;
        (system->v)[4] = -0.5;
        (system->v)[5] = 0.0;
        
        (system->m)[0] = 1.0 / system->G;
        (system->m)[1] = 1.0 / system->G;
    }
    else if (strcmp(system_name, "eccentric_binary_orbit") == 0) 
    {
        const int system_num_particles = 2;
        if (!is_memory_initialized)
        {
            ErrorStatus error_status = WRAP_TRACEBACK(get_initialized_system(system, system_num_particles));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
        else if (system->num_particles < system_num_particles)
        {
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initialized system is not big enough for the built-in system");
        }


        (system->x)[0] = 1.0;
        (system->x)[1] = 0.0;
        (system->x)[2] = 0.0;

        (system->x)[3] = -1.25;
        (system->x)[4] = 0.0;
        (system->x)[5] = 0.0; 

        (system->v)[0] = 0.0;
        (system->v)[1] = 0.5;
        (system->v)[2] = 0.0;

        (system->v)[3] = 0.0;
        (system->v)[4] = -0.625;
        (system->v)[5] = 0.0;
        
        (system->m)[0] = 1.0 / system->G;
        (system->m)[1] = 0.8 / system->G;
    }
    else if (strcmp(system_name, "3d_helix") == 0) 
    {
        const int system_num_particles = 3;
        if (!is_memory_initialized)
        {
            ErrorStatus error_status = WRAP_TRACEBACK(get_initialized_system(system, system_num_particles));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
        else if (system->num_particles < system_num_particles)
        {
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initialized system is not big enough for the built-in system");
        }

        (system->x)[0] = 0.0;
        (system->x)[1] = 0.0;
        (system->x)[2] = -1.0;

        (system->x)[3] = -sqrt(3.0) / 2.0;
        (system->x)[4] = 0.0;
        (system->x)[5] = 0.5;

        (system->x)[6] = sqrt(3.0) / 2.0;
        (system->x)[7] = 0.0;
        (system->x)[8] = 0.5;

        double v0 = sqrt(1.0 / sqrt(3));

        (system->v)[0] = -v0;
        (system->v)[1] = 0.5;
        (system->v)[2] = 0.0;

        (system->v)[3] = 0.5 * v0;
        (system->v)[4] = 0.5;
        (system->v)[5] = (sqrt(3.0) / 2.0) * v0;

        (system->v)[6] = 0.5 * v0;
        (system->v)[7] = 0.5;
        (system->v)[8] = -(sqrt(3.0) / 2.0) * v0;
        
        (system->m)[0] = 1.0 / system->G;
        (system->m)[1] = 1.0 / system->G;
        (system->m)[2] = 1.0 / system->G;
    }
    else if (strcmp(system_name, "sun_earth_moon") == 0) 
    {
        const int system_num_particles = 3;
        if (!is_memory_initialized)
        {
            ErrorStatus error_status = WRAP_TRACEBACK(get_initialized_system(system, system_num_particles));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
        else if (system->num_particles < system_num_particles)
        {
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initialized system is not big enough for the built-in system");
        }

        (system->m)[0] = MASS_SUN;
        (system->m)[1] = MASS_EARTH;
        (system->m)[2] = MASS_MOON;
        
        double R_CM[3];
        double V_CM[3];
        const double M = (system->m)[0] + (system->m)[1] + (system->m)[2];
        for (int i = 0; i < 3; i++)
        {
            R_CM[i] = 1 / M * ((system->m)[0] * POS_SUN[i] + (system->m)[1] * POS_EARTH[i] + (system->m)[2] * POS_MOON[i]);
            V_CM[i] = 1 / M * ((system->m)[0] * VEL_SUN[i] + (system->m)[1] * VEL_EARTH[i] + (system->m)[2] * VEL_MOON[i]);
        }

        for (int i = 0; i < 3; i++)
        {
            (system->x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
            (system->x)[1 * 3 + i] = POS_EARTH[i] - R_CM[i];
            (system->x)[2 * 3 + i] = POS_MOON[i] - R_CM[i];

            (system->v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
            (system->v)[1 * 3 + i] = VEL_EARTH[i] - V_CM[i];
            (system->v)[2 * 3 + i] = VEL_MOON[i] - V_CM[i];
        }
    }
    else if (strcmp(system_name, "figure-8") == 0) 
    {
        const int system_num_particles = 3;
        if (!is_memory_initialized)
        {
            ErrorStatus error_status = WRAP_TRACEBACK(get_initialized_system(system, system_num_particles));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
        else if (system->num_particles < system_num_particles)
        {
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initialized system is not big enough for the built-in system");
        }


        (system->x)[0] = 0.970043;
        (system->x)[1] = -0.24308753;
        (system->x)[2] = 0.0;

        (system->x)[3] = -0.970043;
        (system->x)[4] = 0.24308753;
        (system->x)[5] = 0.0;

        (system->x)[6] = 0.0;
        (system->x)[7] = 0.0;
        (system->x)[8] = 0.0;

        (system->v)[0] = 0.466203685;
        (system->v)[1] = 0.43236573;
        (system->v)[2] = 0.0;

        (system->v)[3] = 0.466203685;
        (system->v)[4] = 0.43236573;
        (system->v)[5] = 0.0;

        (system->v)[6] = -0.93240737;
        (system->v)[7] = -0.86473146;
        (system->v)[8] = 0.0;

        (system->m)[0] = 1.0 / system->G;
        (system->m)[1] = 1.0 / system->G;
        (system->m)[2] = 1.0 / system->G;
    }
    else if (strcmp(system_name, "pyth-3-body") == 0) 
    {
        const int system_num_particles = 3;
        if (!is_memory_initialized)
        {
            ErrorStatus error_status = WRAP_TRACEBACK(get_initialized_system(system, system_num_particles));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
        else if (system->num_particles < system_num_particles)
        {
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initialized system is not big enough for the built-in system");
        }

        (system->x)[0] = 1.0;
        (system->x)[1] = 3.0;
        (system->x)[2] = 0.0;

        (system->x)[3] = -2.0;
        (system->x)[4] = -1.0;
        (system->x)[5] = 0.0;

        (system->x)[6] = 1.0;
        (system->x)[7] = -1.0;
        (system->x)[8] = 0.0;

        (system->v)[0] = 0.0;
        (system->v)[1] = 0.0;
        (system->v)[2] = 0.0;

        (system->v)[3] = 0.0;
        (system->v)[4] = 0.0;
        (system->v)[5] = 0.0;

        (system->v)[6] = 0.0;
        (system->v)[7] = 0.0;
        (system->v)[8] = 0.0;

        (system->m)[0] = 3.0 / system->G;
        (system->m)[1] = 4.0 / system->G;
        (system->m)[2] = 5.0 / system->G;
    }
    else if (strcmp(system_name, "solar_system") == 0) 
    {
        const int system_num_particles = 9;
        if (!is_memory_initialized)
        {
            ErrorStatus error_status = WRAP_TRACEBACK(get_initialized_system(system, system_num_particles));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
        else if (system->num_particles < system_num_particles)
        {
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initialized system is not big enough for the built-in system");
        }

        (system->m)[0] = MASS_SUN;
        (system->m)[1] = MASS_MERCURY;
        (system->m)[2] = MASS_VENUS;
        (system->m)[3] = MASS_EARTH;
        (system->m)[4] = MASS_MARS;
        (system->m)[5] = MASS_JUPITER;
        (system->m)[6] = MASS_SATURN;
        (system->m)[7] = MASS_URANUS;
        (system->m)[8] = MASS_NEPTUNE;

        double R_CM[3];
        double V_CM[3];
        const double M = (
            (system->m)[0] + (system->m)[1] + (system->m)[2] + (system->m)[3] + (system->m)[4] 
            + (system->m)[5] + (system->m)[6] + (system->m)[7] + (system->m)[8]
        );

        for (int i = 0; i < 3; i++)
        {
            R_CM[i] = 1 / M * (
                (system->m)[0] * POS_SUN[i] + (system->m)[1] * POS_MERCURY[i] + (system->m)[2] * POS_VENUS[i]
                + (system->m)[3] * POS_EARTH[i] + (system->m)[4] * POS_MARS[i] + (system->m)[5] * POS_JUPITER[i]
                + (system->m)[6] * POS_SATURN[i] + (system->m)[7] * POS_URANUS[i] + (system->m)[8] * POS_NEPTUNE[i]
            );

            V_CM[i] = 1 / M * (
                (system->m)[0] * VEL_SUN[i] + (system->m)[1] * VEL_MERCURY[i] + (system->m)[2] * VEL_VENUS[i]
                + (system->m)[3] * VEL_EARTH[i] + (system->m)[4] * VEL_MARS[i] + (system->m)[5] * VEL_JUPITER[i]
                + (system->m)[6] * VEL_SATURN[i] + (system->m)[7] * VEL_URANUS[i] + (system->m)[8] * VEL_NEPTUNE[i]
            );
        }

        for (int i = 0; i < 3; i++)
        {
            (system->x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
            (system->x)[1 * 3 + i] = POS_MERCURY[i] - R_CM[i];
            (system->x)[2 * 3 + i] = POS_VENUS[i] - R_CM[i];
            (system->x)[3 * 3 + i] = POS_EARTH[i] - R_CM[i];
            (system->x)[4 * 3 + i] = POS_MARS[i] - R_CM[i];
            (system->x)[5 * 3 + i] = POS_JUPITER[i] - R_CM[i];
            (system->x)[6 * 3 + i] = POS_SATURN[i] - R_CM[i];
            (system->x)[7 * 3 + i] = POS_URANUS[i] - R_CM[i];
            (system->x)[8 * 3 + i] = POS_NEPTUNE[i] - R_CM[i];

            (system->v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
            (system->v)[1 * 3 + i] = VEL_MERCURY[i] - V_CM[i];
            (system->v)[2 * 3 + i] = VEL_VENUS[i] - V_CM[i];
            (system->v)[3 * 3 + i] = VEL_EARTH[i] - V_CM[i];
            (system->v)[4 * 3 + i] = VEL_MARS[i] - V_CM[i];
            (system->v)[5 * 3 + i] = VEL_JUPITER[i] - V_CM[i];
            (system->v)[6 * 3 + i] = VEL_SATURN[i] - V_CM[i];
            (system->v)[7 * 3 + i] = VEL_URANUS[i] - V_CM[i];
            (system->v)[8 * 3 + i] = VEL_NEPTUNE[i] - V_CM[i];
        }
    }
    else if (strcmp(system_name, "solar_system_plus") == 0) 
    {
        const int system_num_particles = 12;
        if (!is_memory_initialized)
        {
            ErrorStatus error_status = WRAP_TRACEBACK(get_initialized_system(system, system_num_particles));
            if (error_status.return_code != GRAV_SUCCESS)
            {
                return error_status;
            }
        }
        else if (system->num_particles < system_num_particles)
        {
            return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Initialized system is not big enough for the built-in system");
        }
        
        (system->m)[0] = MASS_SUN;
        (system->m)[1] = MASS_MERCURY;
        (system->m)[2] = MASS_VENUS;
        (system->m)[3] = MASS_EARTH;
        (system->m)[4] = MASS_MARS;
        (system->m)[5] = MASS_JUPITER;
        (system->m)[6] = MASS_SATURN;
        (system->m)[7] = MASS_URANUS;
        (system->m)[8] = MASS_NEPTUNE;
        (system->m)[9] = MASS_PLUTO;
        (system->m)[10] = MASS_CERES;
        (system->m)[11] = MASS_VESTA;

        double R_CM[3];
        double V_CM[3];
        const double M = (
            (system->m)[0] + (system->m)[1] + (system->m)[2] + (system->m)[3] + (system->m)[4] + (system->m)[5] 
            + (system->m)[6] + (system->m)[7] + (system->m)[8] + (system->m)[9] + (system->m)[10] + (system->m)[11]
        );

        for (int i = 0; i < 3; i++)
        {
            R_CM[i] = 1 / M * (
                (system->m)[0] * POS_SUN[i] + (system->m)[1] * POS_MERCURY[i] + (system->m)[2] * POS_VENUS[i]
                + (system->m)[3] * POS_EARTH[i] + (system->m)[4] * POS_MARS[i] + (system->m)[5] * POS_JUPITER[i]
                + (system->m)[6] * POS_SATURN[i] + (system->m)[7] * POS_URANUS[i] + (system->m)[8] * POS_NEPTUNE[i]
                + (system->m)[9] * POS_PLUTO[i] + (system->m)[10] * POS_CERES[i] + (system->m)[11] * POS_VESTA[i]
            );

            V_CM[i] = 1 / M * (
                (system->m)[0] * VEL_SUN[i] + (system->m)[1] * VEL_MERCURY[i] + (system->m)[2] * VEL_VENUS[i]
                + (system->m)[3] * VEL_EARTH[i] + (system->m)[4] * VEL_MARS[i] + (system->m)[5] * VEL_JUPITER[i]
                + (system->m)[6] * VEL_SATURN[i] + (system->m)[7] * VEL_URANUS[i] + (system->m)[8] * VEL_NEPTUNE[i]
                + (system->m)[9] * VEL_PLUTO[i] + (system->m)[10] * VEL_CERES[i] + (system->m)[11] * VEL_VESTA [i]
            );
        }

        for (int i = 0; i < 3; i++)
        {
            (system->x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
            (system->x)[1 * 3 + i] = POS_MERCURY[i] - R_CM[i];
            (system->x)[2 * 3 + i] = POS_VENUS[i] - R_CM[i];
            (system->x)[3 * 3 + i] = POS_EARTH[i] - R_CM[i];
            (system->x)[4 * 3 + i] = POS_MARS[i] - R_CM[i];
            (system->x)[5 * 3 + i] = POS_JUPITER[i] - R_CM[i];
            (system->x)[6 * 3 + i] = POS_SATURN[i] - R_CM[i];
            (system->x)[7 * 3 + i] = POS_URANUS[i] - R_CM[i];
            (system->x)[8 * 3 + i] = POS_NEPTUNE[i] - R_CM[i];
            (system->x)[9 * 3 + i] = POS_PLUTO[i] - R_CM[i];
            (system->x)[10 * 3 + i] = POS_CERES[i] - R_CM[i];
            (system->x)[11 * 3 + i] = POS_VESTA[i] - R_CM[i];

            (system->v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
            (system->v)[1 * 3 + i] = VEL_MERCURY[i] - V_CM[i];
            (system->v)[2 * 3 + i] = VEL_VENUS[i] - V_CM[i];
            (system->v)[3 * 3 + i] = VEL_EARTH[i] - V_CM[i];
            (system->v)[4 * 3 + i] = VEL_MARS[i] - V_CM[i];
            (system->v)[5 * 3 + i] = VEL_JUPITER[i] - V_CM[i];
            (system->v)[6 * 3 + i] = VEL_SATURN[i] - V_CM[i];
            (system->v)[7 * 3 + i] = VEL_URANUS[i] - V_CM[i];
            (system->v)[8 * 3 + i] = VEL_NEPTUNE[i] - V_CM[i];
            (system->v)[9 * 3 + i] = VEL_PLUTO[i] - V_CM[i];
            (system->v)[10 * 3 + i] = VEL_CERES[i] - V_CM[i];
            (system->v)[11 * 3 + i] = VEL_VESTA[i] - V_CM[i];
        }
    }
    else
    {
        return WRAP_RAISE_ERROR_FMT(
            GRAV_VALUE_ERROR,
            "System name not recognized. Got: \"%s\".",
            system_name
        );
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus system_set_center_of_mass_zero(System *__restrict system)
{
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    /* Declare variables */
    const double num_particles = system->num_particles;
    double *__restrict x = system->x;
    double *__restrict m = system->m;
    if (!x || !m)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System member is NULL");
    }
    
    double R_CM[3] = {0.0, 0.0, 0.0};
    double total_mass = 0.0;

    /* Compute the center of mass and total mass */
    for (int i = 0; i < num_particles; i++)
    {
        R_CM[0] += m[i] * x[i * 3 + 0];
        R_CM[1] += m[i] * x[i * 3 + 1];
        R_CM[2] += m[i] * x[i * 3 + 2];
        total_mass += m[i];
    }

    if (total_mass <= 0.0)
    {
        return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Total mass is non-positive");
    }

    R_CM[0] /= total_mass;
    R_CM[1] /= total_mass;
    R_CM[2] /= total_mass;

    if (!isfinite(R_CM[0]) || !isfinite(R_CM[1]) || !isfinite(R_CM[2]))
    {
        return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Invalid value for center of mass");
    }

    /* Set the center of mass to zero */
    for (int i = 0; i < num_particles; i++)
    {
        x[i * 3 + 0] -= R_CM[0];
        x[i * 3 + 1] -= R_CM[1];
        x[i * 3 + 2] -= R_CM[2];
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus system_set_total_momentum_zero(System *__restrict system)
{
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    /* Declare variables */
    const double num_particles = system->num_particles;
    double *__restrict v = system->v;
    double *__restrict m = system->m;
    if (!v || !m)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System member is NULL");
    }

    double V_CM[3] = {0.0, 0.0, 0.0};
    double total_mass = 0.0;

    /* Compute the center of mass and total mass */
    for (int i = 0; i < num_particles; i++)
    {
        V_CM[0] += m[i] * v[i * 3 + 0];
        V_CM[1] += m[i] * v[i * 3 + 1];
        V_CM[2] += m[i] * v[i * 3 + 2];
        total_mass += m[i];
    }

    if (total_mass <= 0.0)
    {
        return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Total mass is non-positive");
    }

    V_CM[0] /= total_mass;
    V_CM[1] /= total_mass;
    V_CM[2] /= total_mass;

    if (!isfinite(V_CM[0]) || !isfinite(V_CM[1]) || !isfinite(V_CM[2]))
    {
        return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Invalid value for V_CM");
    }

    /* Set the center of mass to zero */
    for (int i = 0; i < num_particles; i++)
    {
        v[i * 3 + 0] -= V_CM[0];
        v[i * 3 + 1] -= V_CM[1];
        v[i * 3 + 2] -= V_CM[2];
    }

    return make_success_error_status();
}

IN_FILE int compare_distance(const void *a, const void *b)
{
    const HelperSystemSortByDistanceStruct *d1 = a;
    const HelperSystemSortByDistanceStruct *d2 = b;
    return (d1->distance > d2->distance) - (d1->distance < d2->distance);
}

WIN32DLL_API ErrorStatus system_sort_by_distance(
    System *__restrict system,
    const int primary_particle_id
)
{
    ErrorStatus error_status;

    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System is NULL");
    }

    const int num_particles = system->num_particles;
    int *__restrict particle_ids = system->particle_ids;
    double *__restrict x = system->x;
    double *__restrict v = system->v;
    double *__restrict m = system->m;
    if (!particle_ids || !x || !v || !m)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System member is NULL");
    }

    /* Find the primary particle index */
    int primary_particle_index = -1;
    if (
        primary_particle_id < system->num_particles
        && system->particle_ids[primary_particle_id] == primary_particle_id
    )
    {
        primary_particle_index = primary_particle_id;
    }
    else
    {
        for (int i = 0; i < num_particles; i++)
        {
            if (particle_ids[i] == primary_particle_id)
            {
                primary_particle_index = i;
                break;
            }
        }
    }
    if (primary_particle_index == -1)
    {
        return WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Primary particle ID not found in system");
    }

    HelperSystemSortByDistanceStruct *helper_arr = malloc(num_particles * sizeof(HelperSystemSortByDistanceStruct));
    if (!helper_arr)
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for helper arrays"
        );
        goto err_helper_arr_malloc;
    }

    for (int i = 0; i < num_particles; i++)
    {
        const double diff_vec[3] = {
            x[i * 3 + 0] - x[primary_particle_index * 3 + 0],
            x[i * 3 + 1] - x[primary_particle_index * 3 + 1],
            x[i * 3 + 2] - x[primary_particle_index * 3 + 2]
        };
        helper_arr[i].distance = vec_norm_3d(diff_vec);
        helper_arr[i].index = i;
    }
    helper_arr[primary_particle_index].distance = 0.0;

    /* Sort the helper array by distance using quick sort */
    qsort(helper_arr, num_particles, sizeof(HelperSystemSortByDistanceStruct), compare_distance);

    /* Allocate temporary arrays */
    int *__restrict new_particle_ids = malloc(num_particles * sizeof(int));
    double *__restrict new_x = malloc(num_particles * 3 * sizeof(double));
    double *__restrict new_v = malloc(num_particles * 3 * sizeof(double));
    double *__restrict new_m = malloc(num_particles * sizeof(double));

    if (!new_particle_ids || !new_x || !new_v || !new_m)
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for new arrays"
        );
        goto err_temp_arr_malloc;
    }

    /* Fill the new arrays with sorted data */
    for (int i = 0; i < num_particles; i++)
    {
        const int index = helper_arr[i].index;
        new_particle_ids[i] = particle_ids[index];
        new_x[i * 3 + 0] = x[index * 3 + 0];
        new_x[i * 3 + 1] = x[index * 3 + 1];
        new_x[i * 3 + 2] = x[index * 3 + 2];
        new_v[i * 3 + 0] = v[index * 3 + 0];
        new_v[i * 3 + 1] = v[index * 3 + 1];
        new_v[i * 3 + 2] = v[index * 3 + 2];
        new_m[i] = m[index];
    }

    /* Copy the new arrays back to the system */
    memcpy(system->particle_ids, new_particle_ids, num_particles * sizeof(int));
    memcpy(system->x, new_x, num_particles * 3 * sizeof(double));
    memcpy(system->v, new_v, num_particles * 3 * sizeof(double));
    memcpy(system->m, new_m, num_particles * sizeof(double));
    
    /* Free the temporary arrays */
    free(helper_arr);

    free(new_particle_ids);
    free(new_x);
    free(new_v);
    free(new_m);

    return make_success_error_status();

err_temp_arr_malloc:
    free(new_particle_ids);
    free(new_x);
    free(new_v);
    free(new_m);
err_helper_arr_malloc:
    free(helper_arr);
    return error_status;
}

WIN32DLL_API void set_periodic_boundary_conditions(CosmologicalSystem *__restrict system)
{
    const int num_particles = system->num_particles;
    double *__restrict x = system->x;
    const double *__restrict box_center = system->box_center;
    const double box_width = system->box_width;

    for (int i = 0; i < num_particles; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            const double normalized_position = x[i * 3 + j] - box_center[j];
            if (normalized_position < -box_width)
            {
                if ((int) (normalized_position / box_width) % 2 != 0)
                {
                    x[i * 3 + j] = fmod(normalized_position, box_width) + box_width;
                }
                else
                {
                    x[i * 3 + j] = fmod(normalized_position, box_width);
                }
                x[i * 3 + j] += box_center[j];
            }
            else if (normalized_position > box_width)
            {
                if ((int) (normalized_position / box_width) % 2 != 0)
                {
                    x[i * 3 + j] = fmod(normalized_position, box_width) - box_width;
                }
                else
                {
                    x[i * 3 + j] = fmod(normalized_position, box_width);
                }
                x[i * 3 + j] += box_center[j];
            }
        }
    }
}
