#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "acceleration.h"
#include "error.h"
#include "gravity_sim.h"
#include "integrator.h"
#include "storing.h"

// IN_FILE int initialize_system(
//     System *__restrict system
// );

/**
 * \brief Check if the integrator is a fixed step size integrator
 * 
 * \param integrator Name of the integrator
 * \param is_fixed_step_size_integrator Pointer to the flag that indicates
 *                                      whether the integrator is a fixed 
 *                                      step size integrator
 * 
 * \retval SUCCESS If the integrator is recognized
 * \retval ERROR_UNKNOWN_INTEGRATOR_METHOD If the integrator is not recognized
 */
IN_FILE int check_is_fixed_step_size_integrator(
    const char *__restrict integrator,
    bool *__restrict is_fixed_step_size_integrator
);

/**
 * \brief Launch simulation helper function
 * 
 * \param system Pointer to the gravitational system
 * \param integrator_param Pointer to the integrator parameters
 * \param acceleration_param Pointer to the acceleration parameters
 * \param storing_param Pointer to the storing parameters
 * \param solutions Pointer to the solutions
 * \param simulation_status Pointer to the simulation status
 * \param settings Pointer to the settings
 * \param simulation_param Pointer to the simulation parameters
 * 
 * \retval SUCCESS If the simulation is successful
 * \retval error code If there is any error
 */
IN_FILE int _launch_simulation(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
);

WIN32DLL_API int launch_simulation_python(
    real *x,
    real *v,
    real *m,
    int objects_count,
    real G,
    const char *integrator,
    real dt,
    real tolerance,
    real initial_dt,
    real whfast_kepler_tol,
    int whfast_kepler_max_iter,
    bool whfast_kepler_auto_remove,
    real whfast_kepler_auto_remove_tol,
    const char *acceleration_method,
    real opening_angle,
    real softening_length,
    int order,
    const char *storing_method,
    const char *flush_path,
    int storing_freq,
    double **sol_state,
    double **sol_time,
    double **sol_dt,
    int64 *sol_size_,
    real *t,
    real *simulation_status_last_dt,
    real *run_time_,
    int verbose,
    bool *is_exit,
    real tf
)
{
    System *system = &(System) {
        .x = x,
        .v = v,
        .m = m,
        .objects_count = objects_count,
        .G = G,
    };
    IntegratorParam *integrator_param = &(IntegratorParam) {
        .integrator = integrator,
        .dt = dt,
        .tolerance = tolerance,
        .initial_dt = initial_dt,
        .whfast_kepler_tol = whfast_kepler_tol,
        .whfast_kepler_max_iter = whfast_kepler_max_iter,
        .whfast_kepler_auto_remove = whfast_kepler_auto_remove,
        .whfast_kepler_auto_remove_tol = whfast_kepler_auto_remove_tol
    };
    AccelerationParam *acceleration_param = &(AccelerationParam) {
        .method = acceleration_method,
        .opening_angle = opening_angle,
        .softening_length = softening_length,
        .order = order,
        .acceleration_method_flag_ = 0
    };
    StoringParam *storing_param = &(StoringParam) {
        .method = storing_method,
        .flush_path = flush_path,
        .storing_freq = storing_freq,
        .storing_method_flag_ = 0,
        .flush_file_ = NULL,
        .max_sol_size_ = 0
    };
    Solutions *solutions = &(Solutions) {
        .sol_state = *sol_state,
        .sol_time = *sol_time,
        .sol_dt = *sol_dt,
        .sol_size_ = sol_size_
    };
    SimulationStatus *simulation_status = &(SimulationStatus) {
        .t = t,
        .dt = 0.0,
        .run_time_ = 0.0
    };
    Settings *settings = &(Settings) {
        .verbose = verbose,
        .is_exit = is_exit
    };
    SimulationParam *simulation_param = &(SimulationParam) {
        .tf = tf
    };

    int return_code = launch_simulation(
        system,
        integrator_param,
        acceleration_param,
        storing_param,
        solutions,
        simulation_status,
        settings,
        simulation_param
    );

    *sol_state = solutions->sol_state;
    *sol_time = solutions->sol_time;
    *sol_dt = solutions->sol_dt;
    *simulation_status_last_dt = simulation_status->dt;
    *run_time_ = simulation_status->run_time_;

    return return_code;
}

WIN32DLL_API int launch_simulation(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
)
{
    int return_code;

    /* Check if the system has NULL array */
    if (!system->x || !system->v || !system->m)
    {
        return_code = ERROR_NULL_SYSTEM_POINTER;
        goto error;
    }

    return_code = get_acceleration_method_flag(
        acceleration_param->method,
        &(acceleration_param->acceleration_method_flag_)
    );
    if (return_code != SUCCESS)
    {
        goto error;
    }

    return_code = get_storing_method_flag(
        storing_param->method,
        &(storing_param->storing_method_flag_)
    );
    if (return_code != SUCCESS)
    {
        goto error;
    }

    return_code = _launch_simulation(
        system,
        integrator_param,
        acceleration_param,
        storing_param,
        solutions,
        simulation_status,
        settings,
        simulation_param
    );
    if (return_code != SUCCESS)
    {
        goto error;
    }
    
    return SUCCESS;

error:
    print_error_msg(return_code);
    return ERROR_SIMULATION_FAILURE;
}

IN_FILE int check_is_fixed_step_size_integrator(
    const char *__restrict integrator,
    bool *__restrict is_fixed_step_size_integrator
)
{
    if (
        strcmp(integrator, "euler") == 0
        || strcmp(integrator, "euler_cromer") == 0
        || strcmp(integrator, "rk4") == 0
        || strcmp(integrator, "leapfrog") == 0
        || strcmp(integrator, "whfast") == 0
    )
    {
        *is_fixed_step_size_integrator = true;
        return SUCCESS;
    }
    else if (
        strcmp(integrator, "rkf45") == 0
        || strcmp(integrator, "dopri") == 0
        || strcmp(integrator, "dverk") == 0
        || strcmp(integrator, "rkf78") == 0
        || strcmp(integrator, "ias15") == 0
    )
    {
        *is_fixed_step_size_integrator = false;
        return SUCCESS;
    }
    else
    {
        return ERROR_UNKNOWN_INTEGRATOR_METHOD;
    }
}

IN_FILE int _launch_simulation(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    StoringParam *storing_param,
    Solutions *solutions,
    SimulationStatus *simulation_status,
    Settings *settings,
    SimulationParam *simulation_param
)
{
    int return_code = SUCCESS;

    bool is_fixed_step_size_integrator;
    return_code = check_is_fixed_step_size_integrator(
        integrator_param->integrator,
        &is_fixed_step_size_integrator
    );
    if (return_code != SUCCESS)
    {
        goto error;
    }

    if (is_fixed_step_size_integrator)
    {
        // Number of time steps
        simulation_param->n_steps_ = simulation_param->tf / integrator_param->dt;
        storing_param->max_sol_size_ = simulation_param->n_steps_ / storing_param->storing_freq;
        storing_param->max_sol_size_ += 1;    // +1 for initial state

        // dt
        simulation_status->dt = integrator_param->dt;
    }
    else
    {
        // Initial buffer size
        storing_param->max_sol_size_ = ADAPTIVE_STEP_SIZE_SOL_BUFFER_SIZE;

        if (integrator_param->initial_dt < 0.0)
        {
            return_code = ERROR_INITIAL_DT_NEGATIVE;
            goto error;
        }
        else if (integrator_param->initial_dt > simulation_param->tf)
        {
            integrator_param->initial_dt = simulation_param->tf;
        }
    }

    /* Output storage */
    *(solutions->sol_size_) = 0;
    *(simulation_status->t) = 0.0;
    if (storing_param->storing_method_flag_ == STORING_METHOD_DEFAULT)
    {
        return_code = allocate_solutions_memory(
            solutions,
            storing_param,
            system->objects_count
        );
    }
    else if (storing_param->storing_method_flag_ == STORING_METHOD_FLUSH)
    {
        return_code = open_flush_file(storing_param);
    }
    if (return_code != SUCCESS)
    {
        goto error;
    }
    
    /* Store initial state */
    return_code = store_solution_step(
        storing_param,
        system,
        simulation_status,
        solutions
    );
    if (return_code != SUCCESS)
    {
        goto error;
    }

    // Launch simulation
    int (*integrator)(
        System *system,
        IntegratorParam *integrator_param,
        AccelerationParam *acceleration_param,
        StoringParam *storing_param,
        Solutions *solutions,
        SimulationStatus *simulation_status,
        Settings *settings,
        SimulationParam *simulation_param
    );

    if (strcmp(integrator_param->integrator, "euler") == 0)
    {
        integrator = euler;
    }
    else if (strcmp(integrator_param->integrator, "euler_cromer") == 0)
    {
        integrator = euler_cromer;
    }
    else if (strcmp(integrator_param->integrator, "rk4") == 0)
    {
        integrator = rk4;
    }
    else if (strcmp(integrator_param->integrator, "leapfrog") == 0)
    {
        integrator = leapfrog;
    }
    else if (
        strcmp(integrator_param->integrator, "rkf45") == 0
        || strcmp(integrator_param->integrator, "dopri") == 0
        || strcmp(integrator_param->integrator, "dverk") == 0
        || strcmp(integrator_param->integrator, "rkf78") == 0
    )
    {
        integrator = rk_embedded;
    }
    else if (strcmp(integrator_param->integrator, "ias15") == 0)
    {
        integrator = ias15;
    }
    else if (strcmp(integrator_param->integrator, "whfast") == 0)
    {
        integrator = whfast;
    }
    else
    {
        return_code = ERROR_UNKNOWN_INTEGRATOR_METHOD;
        goto error;
    }

    clock_t start_time;
    clock_t end_time;
    start_time = clock();
    return_code = integrator(
        system,
        integrator_param,
        acceleration_param,
        storing_param,
        solutions,
        simulation_status,
        settings,
        simulation_param
    );
    if (return_code != SUCCESS)
    {
        goto error;
    }
    end_time = clock();
    simulation_status->run_time_ = (real) (end_time - start_time) / CLOCKS_PER_SEC;

    /* Close flush file */
    if (storing_param->storing_method_flag_ == STORING_METHOD_FLUSH)
    {
        return_code = close_flush_file(storing_param);
        if (return_code != SUCCESS)
        {
            goto error;
        }
    }

    return SUCCESS;

error:
    return return_code;
}

/* Commented as initialize_system() is not used */
// IN_FILE int initialize_system(
//     System *__restrict system
// )
// {
//     const char *__restrict system_name = system->system_name_to_be_initialized;
//     if (!system_name)
//     {
//         return ERROR_INITIALIZE_SYSTEM_NAME_NULL;
//     }

//     /* Pre-defined constants */

//     // Conversion factor from km^3 s^-2 to AU^3 d^-2
//     real CONVERSION_FACTOR = ((real) 86400.0L * 86400.0L) / (149597870.7L * 149597870.7L * 149597870.7L);
//     // GM values (km^3 s^-2)
//     // ref: https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf
//     const real GM_SI_SUN = 132712440041.279419L;
//     const real GM_SI_MERCURY = 22031.868551L;
//     const real GM_SI_VENUS = 324858.592000L;
//     const real GM_SI_EARTH = 398600.435507L;
//     const real GM_SI_MARS = 42828.375816L;
//     const real GM_SI_JUPITER = 126712764.100000L;
//     const real GM_SI_SATURN = 37940584.841800L;
//     const real GM_SI_URANUS = 5794556.400000L;
//     const real GM_SI_NEPTUNE = 6836527.100580L;
//     const real GM_SI_MOON = 4902.800118L;
//     const real GM_SI_PLUTO = 975.500000L;
//     const real GM_SI_CERES = 62.62890L;
//     const real GM_SI_VESTA = 17.288245L;

//     // GM values (AU^3 d^-2)
//     const real GM_SUN = 132712440041.279419L * CONVERSION_FACTOR;
//     // const real GM_MERCURY = 22031.868551L * CONVERSION_FACTOR;
//     // const real GM_VENUS = 324858.592000L * CONVERSION_FACTOR;
//     // const real GM_EARTH = 398600.435507L * CONVERSION_FACTOR;
//     // const real GM_MARS = 42828.375816L * CONVERSION_FACTOR;
//     // const real GM_JUPITER = 126712764.100000L * CONVERSION_FACTOR;
//     // const real GM_SATURN = 37940584.841800L * CONVERSION_FACTOR;
//     // const real GM_URANUS = 5794556.400000L * CONVERSION_FACTOR;
//     // const real GM_NEPTUNE = 6836527.100580L * CONVERSION_FACTOR;
//     // const real GM_MOON = 4902.800118L * CONVERSION_FACTOR;
//     // const real GM_PLUTO = 975.500000L * CONVERSION_FACTOR;
//     // const real GM_CERES = 62.62890L * CONVERSION_FACTOR;
//     // const real GM_VESTA = 17.288245L * CONVERSION_FACTOR;

//     // Solar system masses (M_sun^-1)
//     const real MASS_SUN = 1.0;
//     const real MASS_MERCURY = GM_SI_MERCURY / GM_SI_SUN;
//     const real MASS_VENUS = GM_SI_VENUS / GM_SI_SUN;
//     const real MASS_EARTH = GM_SI_EARTH / GM_SI_SUN;
//     const real MASS_MARS = GM_SI_MARS / GM_SI_SUN;
//     const real MASS_JUPITER = GM_SI_JUPITER / GM_SI_SUN;
//     const real MASS_SATURN = GM_SI_SATURN / GM_SI_SUN;
//     const real MASS_URANUS = GM_SI_URANUS / GM_SI_SUN;
//     const real MASS_NEPTUNE = GM_SI_NEPTUNE / GM_SI_SUN;
//     const real MASS_MOON = GM_SI_MOON / GM_SI_SUN;
//     const real MASS_PLUTO = GM_SI_PLUTO / GM_SI_SUN;
//     const real MASS_CERES = GM_SI_CERES / GM_SI_SUN;
//     const real MASS_VESTA = GM_SI_VESTA / GM_SI_SUN;

//     // Gravitational constant (kg^-1 m^3 s^-2):
//     // const real G_SI = 6.67430e-11;
//     // Gravitational constant (M_sun^-1 AU^3 d^-2):
//     system->G = GM_SUN;

//     /*
//     * Solar system position and velocities data
//     * Units: AU-D
//     * Coordinate center: Solar System Barycenter
//     * Data dated on A.D. 2024-Jan-01 00:00:00.0000 TDB
//     * Computational data generated by NASA JPL Horizons System https://ssd.jpl.nasa.gov/horizons/
//     */
//     const real POS_SUN[3] = {-7.967955691533730e-03L, -2.906227441573178e-03L, 2.103054301547123e-04L};
//     const real POS_MERCURY[3] = {-2.825983269538632e-01L, 1.974559795958082e-01L, 4.177433558063677e-02L};
//     const real POS_VENUS[3] = {-7.232103701666379e-01L, -7.948302026312400e-02L, 4.042871428174315e-02L};
//     const real POS_EARTH[3] = {-1.738192017257054e-01L, 9.663245550235138e-01L, 1.553901854897183e-04L};
//     const real POS_MARS[3] = {-3.013262392582653e-01L, -1.454029331393295e00L, -2.300531433991428e-02L};
//     const real POS_JUPITER[3] = {3.485202469657674e00L, 3.552136904413157e00L, -9.271035442798399e-02L};
//     const real POS_SATURN[3] = {8.988104223143450e00L, -3.719064854634689e00L, -2.931937777323593e-01L};
//     const real POS_URANUS[3] = {1.226302417897505e01L, 1.529738792480545e01L, -1.020549026883563e-01L};
//     const real POS_NEPTUNE[3] = {2.983501460984741e01L, -1.793812957956852e00L, -6.506401132254588e-01L};
//     const real POS_MOON[3] = {-1.762788124769829e-01L, 9.674377513177153e-01L, 3.236901585768862e-04L};
//     const real POS_PLUTO[3] = {1.720200478843485e01L, -3.034155683573043e01L, -1.729127607100611e00L};
//     const real POS_CERES[3] = {-1.103880510367569e00L, -2.533340440444230e00L, 1.220283937721780e-01L};
//     const real POS_VESTA[3] = {-8.092549658731499e-02L, 2.558381434460076e00L, -6.695836142398572e-02L};

//     const real VEL_SUN[3] = {4.875094764261564e-06L, -7.057133213976680e-06L, -4.573453713094512e-08L};
//     const real VEL_MERCURY[3] = {-2.232165900189702e-02L, -2.157207103176252e-02L, 2.855193410495743e-04L};
//     const real VEL_VENUS[3] = {2.034068201002341e-03L, -2.020828626592994e-02L, -3.945639843855159e-04L};
//     const real VEL_EARTH[3] = {-1.723001232538228e-02L, -2.967721342618870e-03L, 6.382125383116755e-07L};
//     const real VEL_MARS[3] = {1.424832259345280e-02L, -1.579236181580905e-03L, -3.823722796161561e-04L};
//     const real VEL_JUPITER[3] = {-5.470970658852281e-03L, 5.642487338479145e-03L, 9.896190602066252e-05L};
//     const real VEL_SATURN[3] = {1.822013845554067e-03L, 5.143470425888054e-03L, -1.617235904887937e-04L};
//     const real VEL_URANUS[3] = {-3.097615358317413e-03L, 2.276781932345769e-03L, 4.860433222241686e-05L};
//     const real VEL_NEPTUNE[3] = {1.676536611817232e-04L, 3.152098732861913e-03L, -6.877501095688201e-05L};
//     const real VEL_MOON[3] = {-1.746667306153906e-02L, -3.473438277358121e-03L, -3.359028758606074e-05L};
//     const real VEL_PLUTO[3] = {2.802810313667557e-03L, 8.492056438614633e-04L, -9.060790113327894e-04L};
//     const real VEL_CERES[3] = {8.978653480111301e-03L, -4.873256528198994e-03L, -1.807162046049230e-03L};
//     const real VEL_VESTA[3] = {-1.017876585480054e-02L, -5.452367109338154e-04L, 1.255870551153315e-03L};
    
//     /* Initialize some variables */
//     // Note that G has been initialized in the pre-defined constants
//     int *objects_count = &(system->objects_count);
//     real **x = &(system->x);
//     real **v = &(system->v);
//     real **m = &(system->m);
//     real *G = &(system->G);

//     if (*x || *v || *m)
//     {
//         return ERROR_INITIALIZE_SYSTEM_MEMORY_NOT_NULL;
//     }

//     // Pre-defined systems
//     if (strcmp(system_name, "circular_binary_orbit") == 0) 
//     {
//         *objects_count = 2;
//         *x = malloc(*objects_count * 3 * sizeof(real));
//         *v = malloc(*objects_count * 3 * sizeof(real));
//         *m = malloc(*objects_count * sizeof(real));

//         if (!*x || !*v || !*m)
//         {
//             goto err_memory;
//         }

//         (*x)[0] = 1.0;
//         (*x)[1] = 0.0;
//         (*x)[2] = 0.0;

//         (*x)[3] = -1.0;
//         (*x)[4] = 0.0;
//         (*x)[5] = 0.0; 

//         (*v)[0] = 0.0;
//         (*v)[1] = 0.5;
//         (*v)[2] = 0.0;

//         (*v)[3] = 0.0;
//         (*v)[4] = -0.5;
//         (*v)[5] = 0.0;
        
//         (*m)[0] = 1.0 / *G;
//         (*m)[1] = 1.0 / *G;

//         return SUCCESS;
//     }
//     else if (strcmp(system_name, "eccentric_binary_orbit") == 0) 
//     {
//         *objects_count = 2;
//         *x = malloc(*objects_count * 3 * sizeof(real));
//         *v = malloc(*objects_count * 3 * sizeof(real));
//         *m = malloc(*objects_count * sizeof(real));

//         if (!*x || !*v || !*m)
//         {
//             goto err_memory;
//         }

//         (*x)[0] = 1.0;
//         (*x)[1] = 0.0;
//         (*x)[2] = 0.0;

//         (*x)[3] = -1.25;
//         (*x)[4] = 0.0;
//         (*x)[5] = 0.0; 

//         (*v)[0] = 0.0;
//         (*v)[1] = 0.5;
//         (*v)[2] = 0.0;

//         (*v)[3] = 0.0;
//         (*v)[4] = -0.625;
//         (*v)[5] = 0.0;
        
//         (*m)[0] = 1.0 / *G;
//         (*m)[1] = 0.8 / *G;

//         return SUCCESS;
//     }
//     else if (strcmp(system_name, "3d_helix") == 0) 
//     {
//         *objects_count = 3;
//         *x = malloc(*objects_count * 3 * sizeof(real));
//         *v = malloc(*objects_count * 3 * sizeof(real));
//         *m = malloc(*objects_count * sizeof(real));

//         if (!*x || !*v || !*m)
//         {
//             goto err_memory;
//         }
        
//         (*x)[0] = 0.0;
//         (*x)[1] = 0.0;
//         (*x)[2] = -1.0;

//         (*x)[3] = -sqrt(3.0) / 2.0;
//         (*x)[4] = 0.0;
//         (*x)[5] = 0.5;

//         (*x)[6] = sqrt(3.0) / 2.0;
//         (*x)[7] = 0.0;
//         (*x)[8] = 0.5;

//         real v0 = sqrt(1.0 / sqrt(3));

//         (*v)[0] = -v0;
//         (*v)[1] = 0.5;
//         (*v)[2] = 0.0;

//         (*v)[3] = 0.5 * v0;
//         (*v)[4] = 0.5;
//         (*v)[5] = (sqrt(3.0) / 2.0) * v0;

//         (*v)[6] = 0.5 * v0;
//         (*v)[7] = 0.5;
//         (*v)[8] = -(sqrt(3.0) / 2.0) * v0;
        
//         (*m)[0] = 1.0 / *G;
//         (*m)[1] = 1.0 / *G;
//         (*m)[2] = 1.0 / *G;

//         return SUCCESS;
//     }
//     else if (strcmp(system_name, "sun_earth_moon") == 0) 
//     {
//         *objects_count = 3;
//         *x = malloc(*objects_count * 3 * sizeof(real));
//         *v = malloc(*objects_count * 3 * sizeof(real));
//         *m = malloc(*objects_count * sizeof(real));

//         if (!*x || !*v || !*m)
//         {
//             goto err_memory;
//         }

//         (*m)[0] = MASS_SUN;
//         (*m)[1] = MASS_EARTH;
//         (*m)[2] = MASS_MOON;
        
//         real R_CM[3];
//         real V_CM[3];
//         const real M = (*m)[0] + (*m)[1] + (*m)[2];
//         for (int i = 0; i < 3; i++)
//         {
//             R_CM[i] = 1 / M * ((*m)[0] * POS_SUN[i] + (*m)[1] * POS_EARTH[i] + (*m)[2] * POS_MOON[i]);
//             V_CM[i] = 1 / M * ((*m)[0] * VEL_SUN[i] + (*m)[1] * VEL_EARTH[i] + (*m)[2] * VEL_MOON[i]);
//         }

//         for (int i = 0; i < 3; i++)
//         {
//             (*x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
//             (*x)[1 * 3 + i] = POS_EARTH[i] - R_CM[i];
//             (*x)[2 * 3 + i] = POS_MOON[i] - R_CM[i];

//             (*v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
//             (*v)[1 * 3 + i] = VEL_EARTH[i] - V_CM[i];
//             (*v)[2 * 3 + i] = VEL_MOON[i] - V_CM[i];
//         }

//         return SUCCESS;
//     }
//     else if (strcmp(system_name, "figure-8") == 0) 
//     {
//         *objects_count = 3;
//         *x = malloc(*objects_count * 3 * sizeof(real));
//         *v = malloc(*objects_count * 3 * sizeof(real));
//         *m = malloc(*objects_count * sizeof(real));

//         if (!*x || !*v || !*m)
//         {
//             goto err_memory;
//         }

//         (*x)[0] = 0.970043;
//         (*x)[1] = -0.24308753;
//         (*x)[2] = 0.0;

//         (*x)[3] = -0.970043;
//         (*x)[4] = 0.24308753;
//         (*x)[5] = 0.0;

//         (*x)[6] = 0.0;
//         (*x)[7] = 0.0;
//         (*x)[8] = 0.0;

//         (*v)[0] = 0.466203685;
//         (*v)[1] = 0.43236573;
//         (*v)[2] = 0.0;

//         (*v)[3] = 0.466203685;
//         (*v)[4] = 0.43236573;
//         (*v)[5] = 0.0;

//         (*v)[6] = -0.93240737;
//         (*v)[7] = -0.86473146;
//         (*v)[8] = 0.0;

//         (*m)[0] = 1.0 / *G;
//         (*m)[1] = 1.0 / *G;
//         (*m)[2] = 1.0 / *G;

//         return SUCCESS;
//     }
//     else if (strcmp(system_name, "pyth-3-body") == 0) 
//     {
//         *objects_count = 3;
//         *x = malloc(*objects_count * 3 * sizeof(real));
//         *v = malloc(*objects_count * 3 * sizeof(real));
//         *m = malloc(*objects_count * sizeof(real));

//         if (!*x || !*v || !*m)
//         {
//             goto err_memory;
//         }

//         (*x)[0] = 1.0;
//         (*x)[1] = 3.0;
//         (*x)[2] = 0.0;

//         (*x)[3] = -2.0;
//         (*x)[4] = -1.0;
//         (*x)[5] = 0.0;

//         (*x)[6] = 1.0;
//         (*x)[7] = -1.0;
//         (*x)[8] = 0.0;

//         (*v)[0] = 0.0;
//         (*v)[1] = 0.0;
//         (*v)[2] = 0.0;

//         (*v)[3] = 0.0;
//         (*v)[4] = 0.0;
//         (*v)[5] = 0.0;

//         (*v)[6] = 0.0;
//         (*v)[7] = 0.0;
//         (*v)[8] = 0.0;

//         (*m)[0] = 3.0 / *G;
//         (*m)[1] = 4.0 / *G;
//         (*m)[2] = 5.0 / *G;

//         return SUCCESS;
//     }
//     else if (strcmp(system_name, "solar_system") == 0) 
//     {
//         *objects_count = 9;
//         *x = malloc(*objects_count * 3 * sizeof(real));
//         *v = malloc(*objects_count * 3 * sizeof(real));
//         *m = malloc(*objects_count * sizeof(real));

//         if (!*x || !*v || !*m)
//         {
//             goto err_memory;
//         }

//         (*m)[0] = MASS_SUN;
//         (*m)[1] = MASS_MERCURY;
//         (*m)[2] = MASS_VENUS;
//         (*m)[3] = MASS_EARTH;
//         (*m)[4] = MASS_MARS;
//         (*m)[5] = MASS_JUPITER;
//         (*m)[6] = MASS_SATURN;
//         (*m)[7] = MASS_URANUS;
//         (*m)[8] = MASS_NEPTUNE;

//         real R_CM[3];
//         real V_CM[3];
//         const real M = (
//             (*m)[0] + (*m)[1] + (*m)[2] + (*m)[3] + (*m)[4] 
//             + (*m)[5] + (*m)[6] + (*m)[7] + (*m)[8]
//         );

//         for (int i = 0; i < 3; i++)
//         {
//             R_CM[i] = 1 / M * (
//                 (*m)[0] * POS_SUN[i] + (*m)[1] * POS_MERCURY[i] + (*m)[2] * POS_VENUS[i]
//                 + (*m)[3] * POS_EARTH[i] + (*m)[4] * POS_MARS[i] + (*m)[5] * POS_JUPITER[i]
//                 + (*m)[6] * POS_SATURN[i] + (*m)[7] * POS_URANUS[i] + (*m)[8] * POS_NEPTUNE[i]
//             );

//             V_CM[i] = 1 / M * (
//                 (*m)[0] * VEL_SUN[i] + (*m)[1] * VEL_MERCURY[i] + (*m)[2] * VEL_VENUS[i]
//                 + (*m)[3] * VEL_EARTH[i] + (*m)[4] * VEL_MARS[i] + (*m)[5] * VEL_JUPITER[i]
//                 + (*m)[6] * VEL_SATURN[i] + (*m)[7] * VEL_URANUS[i] + (*m)[8] * VEL_NEPTUNE[i]
//             );
//         }

//         for (int i = 0; i < 3; i++)
//         {
//             (*x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
//             (*x)[1 * 3 + i] = POS_MERCURY[i] - R_CM[i];
//             (*x)[2 * 3 + i] = POS_VENUS[i] - R_CM[i];
//             (*x)[3 * 3 + i] = POS_EARTH[i] - R_CM[i];
//             (*x)[4 * 3 + i] = POS_MARS[i] - R_CM[i];
//             (*x)[5 * 3 + i] = POS_JUPITER[i] - R_CM[i];
//             (*x)[6 * 3 + i] = POS_SATURN[i] - R_CM[i];
//             (*x)[7 * 3 + i] = POS_URANUS[i] - R_CM[i];
//             (*x)[8 * 3 + i] = POS_NEPTUNE[i] - R_CM[i];

//             (*v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
//             (*v)[1 * 3 + i] = VEL_MERCURY[i] - V_CM[i];
//             (*v)[2 * 3 + i] = VEL_VENUS[i] - V_CM[i];
//             (*v)[3 * 3 + i] = VEL_EARTH[i] - V_CM[i];
//             (*v)[4 * 3 + i] = VEL_MARS[i] - V_CM[i];
//             (*v)[5 * 3 + i] = VEL_JUPITER[i] - V_CM[i];
//             (*v)[6 * 3 + i] = VEL_SATURN[i] - V_CM[i];
//             (*v)[7 * 3 + i] = VEL_URANUS[i] - V_CM[i];
//             (*v)[8 * 3 + i] = VEL_NEPTUNE[i] - V_CM[i];
//         }

//         return SUCCESS;
//     }
//     else if (strcmp(system_name, "solar_system_plus") == 0) 
//     {
//         *objects_count = 12;
//         *x = malloc(*objects_count * 3 * sizeof(real));
//         *v = malloc(*objects_count * 3 * sizeof(real));
//         *m = malloc(*objects_count * sizeof(real));

//         if (!*x || !*v || !*m)
//         {
//             goto err_memory;
//         }
        
//         (*m)[0] = MASS_SUN;
//         (*m)[1] = MASS_MERCURY;
//         (*m)[2] = MASS_VENUS;
//         (*m)[3] = MASS_EARTH;
//         (*m)[4] = MASS_MARS;
//         (*m)[5] = MASS_JUPITER;
//         (*m)[6] = MASS_SATURN;
//         (*m)[7] = MASS_URANUS;
//         (*m)[8] = MASS_NEPTUNE;
//         (*m)[9] = MASS_PLUTO;
//         (*m)[10] = MASS_CERES;
//         (*m)[11] = MASS_VESTA;

//         real R_CM[3];
//         real V_CM[3];
//         const real M = (
//             (*m)[0] + (*m)[1] + (*m)[2] + (*m)[3] + (*m)[4] + (*m)[5] 
//             + (*m)[6] + (*m)[7] + (*m)[8] + (*m)[9] + (*m)[10] + (*m)[11]
//         );

//         for (int i = 0; i < 3; i++)
//         {
//             R_CM[i] = 1 / M * (
//                 (*m)[0] * POS_SUN[i] + (*m)[1] * POS_MERCURY[i] + (*m)[2] * POS_VENUS[i]
//                 + (*m)[3] * POS_EARTH[i] + (*m)[4] * POS_MARS[i] + (*m)[5] * POS_JUPITER[i]
//                 + (*m)[6] * POS_SATURN[i] + (*m)[7] * POS_URANUS[i] + (*m)[8] * POS_NEPTUNE[i]
//                 + (*m)[9] * POS_PLUTO[i] + (*m)[10] * POS_CERES[i] + (*m)[11] * POS_VESTA[i]
//             );

//             V_CM[i] = 1 / M * (
//                 (*m)[0] * VEL_SUN[i] + (*m)[1] * VEL_MERCURY[i] + (*m)[2] * VEL_VENUS[i]
//                 + (*m)[3] * VEL_EARTH[i] + (*m)[4] * VEL_MARS[i] + (*m)[5] * VEL_JUPITER[i]
//                 + (*m)[6] * VEL_SATURN[i] + (*m)[7] * VEL_URANUS[i] + (*m)[8] * VEL_NEPTUNE[i]
//                 + (*m)[9] * VEL_PLUTO[i] + (*m)[10] * VEL_CERES[i] + (*m)[11] * VEL_VESTA [i]
//             );
//         }

//         for (int i = 0; i < 3; i++)
//         {
//             (*x)[0 * 3 + i] = POS_SUN[i] - R_CM[i];
//             (*x)[1 * 3 + i] = POS_MERCURY[i] - R_CM[i];
//             (*x)[2 * 3 + i] = POS_VENUS[i] - R_CM[i];
//             (*x)[3 * 3 + i] = POS_EARTH[i] - R_CM[i];
//             (*x)[4 * 3 + i] = POS_MARS[i] - R_CM[i];
//             (*x)[5 * 3 + i] = POS_JUPITER[i] - R_CM[i];
//             (*x)[6 * 3 + i] = POS_SATURN[i] - R_CM[i];
//             (*x)[7 * 3 + i] = POS_URANUS[i] - R_CM[i];
//             (*x)[8 * 3 + i] = POS_NEPTUNE[i] - R_CM[i];
//             (*x)[9 * 3 + i] = POS_PLUTO[i] - R_CM[i];
//             (*x)[10 * 3 + i] = POS_CERES[i] - R_CM[i];
//             (*x)[11 * 3 + i] = POS_VESTA[i] - R_CM[i];

//             (*v)[0 * 3 + i] = VEL_SUN[i] - V_CM[i];
//             (*v)[1 * 3 + i] = VEL_MERCURY[i] - V_CM[i];
//             (*v)[2 * 3 + i] = VEL_VENUS[i] - V_CM[i];
//             (*v)[3 * 3 + i] = VEL_EARTH[i] - V_CM[i];
//             (*v)[4 * 3 + i] = VEL_MARS[i] - V_CM[i];
//             (*v)[5 * 3 + i] = VEL_JUPITER[i] - V_CM[i];
//             (*v)[6 * 3 + i] = VEL_SATURN[i] - V_CM[i];
//             (*v)[7 * 3 + i] = VEL_URANUS[i] - V_CM[i];
//             (*v)[8 * 3 + i] = VEL_NEPTUNE[i] - V_CM[i];
//             (*v)[9 * 3 + i] = VEL_PLUTO[i] - V_CM[i];
//             (*v)[10 * 3 + i] = VEL_CERES[i] - V_CM[i];
//             (*v)[11 * 3 + i] = VEL_VESTA[i] - V_CM[i];
//         }

//         return SUCCESS;
//     }

//     return ERROR_UNKNOWN_INITIALIZE_SYSTEM_NAME;

// err_memory:
//     free(*x);
//     free(*v);
//     free(*m);
//     return ERROR_INITIALIZE_SYSTEM_MEMORY_ALLOC;
// }
