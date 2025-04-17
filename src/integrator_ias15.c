/**
 * \file integrator_ias15.c
 * \brief Implementation for IAS15 integrator
 * 
 * \author Ching-Yin Ng
 * 
 * \ref J. Roa, et al. Moving Planets Around: An Introduction to
 *   N-Body Simulations Applied to Exoplanetary Systems*, MIT
 *   Press, 2020
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration.h"
#include "common.h"
#include "error.h"
#include "math_functions.h"
#include "output.h"
#include "progress_bar.h"
#include "settings.h"
#include "system.h"

#define IAS15_PREDICTOR_CORRECTOR_MAX_ITER 12

/**
 * \brief Initialize the radau spacing nodes for IAS15
 * 
 * \param[out] nodes 1D array of size 8 to be modified
 */
IN_FILE void initialize_radau_spacing(double *restrict nodes);

/**
 * \brief Initialize the auxiliary coefficients aux_c for IAS15
 * 
 * \param[out] aux_c 1D array of length 49 to be modified
 */ 
IN_FILE void initialize_aux_c(double *restrict aux_c);

/**
 * \brief Initialize auxiliary coefficients aux_r for IAS15
 * 
 * \param[out] aux_r 1D array of size 64 to be modified
 */
IN_FILE void initialize_aux_r(double *restrict aux_r);

/**
 * \brief Calculate the initial time step for IAS15 integrator
 * 
 * \param[out] initial_dt Pointer to the initial time step
 * \param[in] power Power of the integrator
 * \param[in] system Pointer to the system
 * \param[in] acceleration_param Pointer to the acceleration parameters
 * \param[in] a Pointer to the acceleration array with already computed values
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if failed to allocate memory for arrays
 * \exception GRAV_VALUE_ERROR if initial_dt is negative
 */
IN_FILE ErrorStatus ias15_initial_dt(
    double *restrict initial_dt,
    const int power,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param,
    const double *restrict a
);

/**
 * \brief Calculate position for the predictor-corrector algorithm
 * 
 * \param[out] x Array of position vectors to be modified
 * \param[in] num_particles Number of particles in the system
 * \param[in] x0 Array of position vectors from the last time step
 * \param[in] v0 Array of velocity vectors from the last time step
 * \param[in] a0 Array of acceleration vectors from the last time step
 * \param[in] node Current node of the predictor-corrector algorithm
 * \param[in] aux_b Auxiliary b array
 * \param[in] dt Current time step of the system
 * \param[in] x_err_comp_sum Array of round off errors of position vectors 
 *                       for compensated summation
 */
IN_FILE void approx_pos_pc(
    double *restrict x,
    const int num_particles,
    const double *restrict x0,
    const double *restrict v0,
    const double *restrict a0,
    const double node,
    const double *restrict aux_b,
    const double dt,
    const double *restrict x_err_comp_sum
);

/**
 * \brief Calculate velocity for the predictor-corrector algorithm
 *
 * \param[out] v Array of velocity vectors to be modified
 * \param[in] num_particles Number of particles in the system
 * \param[in] v0 Array of velocity vectors from the last time step
 * \param[in] a0 Array of acceleration vectors from the last time step
 * \param[in] node Current node of the predictor-corrector algorithm
 * \param[in] aux_b Auxiliary b array
 * \param[in] dt Current time step of the system
 * \param[in] v_err_comp_sum Array of round off errors of velocity vectors 
 *                       for compensated summation
 */
IN_FILE void approx_vel_pc(
    double *restrict v,
    const int num_particles,
    const double *restrict v0,
    const double *restrict a0,
    const double node,
    const double *restrict aux_b,
    const double dt,
    const double *restrict v_err_comp_sum
);

/**
 * \brief Calculate position for the next time step
 * 
 * \param[out] x Array of position vectors to be modified
 * \param[out] temp_x_err_comp_sum Temporary array for compensated summation
 * \param[in] num_particles Number of particles in the system
 * \param[in] x0 Array of position vectors from the last time step
 * \param[in] v0 Array of velocity vectors from the last time step
 * \param[in] a0 Array of acceleration vectors from the last time step
 * \param[in] aux_b Auxiliary coefficients b
 * \param[in] dt Current time step of the system
 */
IN_FILE void approx_pos_step(
    double *restrict x,
    double *restrict temp_x_err_comp_sum,
    const int num_particles,
    const double *restrict x0,
    const double *restrict v0,
    const double *restrict a0,
    const double *restrict aux_b,
    const double dt
);

/**
 * \brief Calculate velocity for the next time step
 * 
 * \param[out] v Array of velocity vectors to be modified
 * \param[out] temp_v_err_comp_sum Temporary array for compensated summation
 * \param[in] num_particles Number of particles in the system
 * \param[in] v0 Array of velocity vectors from the last time step
 * \param[in] a0 Array of acceleration vectors from the last time step
 * \param[in] aux_b Auxiliary coefficients b
 * \param[in] dt Current time step of the system
 */
IN_FILE void approx_vel_step(
    double *restrict v,
    double *restrict temp_v_err_comp_sum,
    const int num_particles,
    const double *restrict v0,
    const double *restrict a0,
    const double *restrict aux_b,
    const double dt
);

/**
 * \brief Calculate the auxiliary coefficients b for IAS15 integrator
 * 
 * \param[out] aux_b Array of auxiliary coefficients b to be modified
 * \param[in] num_particles Number of particles in the system
 * \param[in] dim_nodes_minus_1 Dimension of nodes minus one
 * \param[in] aux_g Array of auxiliary coefficients g
 * \param[in] aux_c Array of auxiliary coefficients c
 * \param[in] i Current iteration of nodes of the predictor-corrector algorithm
 */
IN_FILE void compute_aux_b(
    double *restrict aux_b,
    const int num_particles,
    const int dim_nodes_minus_1,
    const double *restrict aux_g,
    const double *restrict aux_c,
    const int i
);

/**
 * \brief Calculate the auxiliary coefficients g for IAS15 integrator
 * 
 * \param[out] aux_g Array of auxiliary coefficients g to be modified
 * \param[in] num_particles Number of particles in the system
 * \param[in] dim_nodes Dimension of nodes
 * \param[in] aux_r Array of auxiliary coefficients r
 * \param[in] aux_a Array of auxiliary accelerations a
 * \param[in] i Current iteration of nodes of the predictor-corrector algorithm
 * \param[in] F Helper array for this function
 */
IN_FILE void compute_aux_g(
    double *restrict aux_g,
    const int num_particles,
    const int dim_nodes,
    const double *restrict aux_r,
    const double *restrict aux_a,
    const int i,
    double *restrict F
);

/**
 * \brief Refine the auxiliary coefficients b for IAS15 integrator
 * 
 * \param[out] aux_b Array of auxiliary coefficients b to be modified
 * \param[out] aux_e Array of auxiliary coefficients e
 * \param[out] delta_aux_b Helper array for this function
 * \param[in] num_particles Number of particles in the system
 * \param[in] dim_nodes_minus_1 Dimension of nodes minus one
 * \param[in] dt Current time step of the system
 * \param[in] dt_new Next Time Step of the system
 * \param[in] refine_flag Helper flag for this function
 */
IN_FILE void refine_aux_b(
    double *restrict aux_b,
    double *restrict aux_e,
    double *restrict delta_aux_b,
    const int num_particles,
    const int dim_nodes_minus_1,
    const double dt,
    const double dt_new,
    const bool refine_flag
);

WIN32DLL_API ErrorStatus ias15(
    System *system,
    IntegratorParam *integrator_param,
    AccelerationParam *acceleration_param,
    OutputParam *output_param,
    SimulationStatus *simulation_status,
    Settings *settings,
    const double tf
)
{
    ErrorStatus error_status;

    /* Declare variables and allocate memory */
    // tolerance
    double tolerance = integrator_param->tolerance;

    // Safety factors for step-size control
    double safety_fac = 0.25;

    // For fixed step integration, choose exponent = 0
    double exponent = 1.0 / 7.0;

    // Tolerance of predictor-corrector algorithm
    double tolerance_pc = 1e-16;

    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;

    bool is_output = (output_param->method != OUTPUT_METHOD_DISABLED);
    int *restrict output_count_ptr = &(output_param->output_count_);
    const double output_interval = output_param->output_interval;
    double next_output_time = output_interval;

    double *restrict t_ptr = &(simulation_status->t);
    int64 *restrict num_steps_ptr = &(simulation_status->num_steps);

    const bool enable_progress_bar = settings->enable_progress_bar;

    double error;
    double error_b7;
    bool accept_step_flag = false;
    bool refine_flag = false;
    const int dim_nodes = 8;
    const int dim_nodes_minus_1 = 7;
    const int dim_nodes_minus_2 = 6;
    double *restrict nodes = malloc(dim_nodes * sizeof(double));
    double *restrict aux_c = calloc(7 * 7, sizeof(double));
    double *restrict aux_r = calloc(8 * 8, sizeof(double));
    double *restrict aux_b0 = calloc((dim_nodes - 1) * num_particles * 3, sizeof(double));
    double *restrict aux_b = calloc((dim_nodes - 1) * num_particles * 3, sizeof(double));
    double *restrict aux_g = calloc((dim_nodes - 1) * num_particles * 3, sizeof(double));
    double *restrict aux_e = calloc((dim_nodes - 1) * num_particles * 3, sizeof(double));
    if (!nodes || !aux_c || !aux_r || !aux_b0 || !aux_b || !aux_g || !aux_e)
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for auxiliary arrays"
        );
        goto err_aux_memory;
    }
    initialize_radau_spacing(nodes);
    initialize_aux_c(aux_c);
    initialize_aux_r(aux_r);

    // Arrays
    double *restrict a = malloc(num_particles * 3 * sizeof(double));
    double *restrict aux_a = calloc(dim_nodes * num_particles * 3, sizeof(double));
    double *restrict x_1 = calloc(num_particles * 3, sizeof(double));
    double *restrict v_1 = calloc(num_particles * 3, sizeof(double));
    double *restrict a_1 = calloc(num_particles * 3, sizeof(double));
    double *restrict delta_b7 = calloc(num_particles * 3, sizeof(double));

    // Array for compute aux_g
    double *restrict F = calloc(8 * num_particles * 3, sizeof(double));

    // Array for refine aux_b
    double *restrict delta_aux_b = calloc(dim_nodes_minus_1 * num_particles * 3, sizeof(double));

    // Arrays for compensated summation
    double *restrict x_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict v_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict temp_x_err_comp_sum = calloc(num_particles * 3, sizeof(double));
    double *restrict temp_v_err_comp_sum = calloc(num_particles * 3, sizeof(double));

    if (
        !a || 
        !aux_a || 
        !x_1 ||
        !v_1 ||
        !a_1 ||
        !delta_b7 ||
        !F ||
        !delta_aux_b ||
        !x_err_comp_sum ||
        !v_err_comp_sum ||
        !temp_x_err_comp_sum ||
        !temp_v_err_comp_sum
    )
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for arrays"
        );
        goto err_memory;
    }

    System temp_system = {
        .num_particles = num_particles,
        .x = x_1,
        .v = v_1,
        .m = system->m,
        .G = system->G
    };

    /* Get initial dt */
    error_status = WRAP_TRACEBACK(acceleration(
        a,
        system,
        acceleration_param
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto err_acc_error;
    }

    double dt;
    double dt_new;
    if (integrator_param->initial_dt > 0.0)
    {
        dt = integrator_param->initial_dt;
    }
    else
    {
        error_status = WRAP_TRACEBACK(ias15_initial_dt(
            &dt,
            15,
            system,
            acceleration_param,
            a
        ));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_initial_dt;
        }

        if (dt > tf)
        {
            dt = tf;
        }
    }

    /* Initial output */
    if (is_output && output_param->output_initial)
    {
        error_status = WRAP_TRACEBACK(output_snapshot(
            output_param,
            system,
            integrator_param,
            acceleration_param,
            simulation_status,
            settings
        ));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_initial_output;
        }
    }

    /* Main Loop */
    ProgressBarParam progress_bar_param;
    if (enable_progress_bar)
    {
        error_status = WRAP_TRACEBACK(start_progress_bar(&progress_bar_param, tf));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_start_progress_bar;
        }
    }

    *t_ptr = 0.0;
    simulation_status->dt = dt;
    *num_steps_ptr = 0;
    while (*t_ptr < tf)
    {
        /* Loop for predictor-corrector algorithm */
        for (int temp = 0; temp < IAS15_PREDICTOR_CORRECTOR_MAX_ITER; temp++)
        {
            for (int i = 0; i < dim_nodes; i++)
            {
                // Estimate position and velocity with current aux_b and nodes
                approx_pos_pc(x_1, num_particles, x, v, a, nodes[i], aux_b, dt, x_err_comp_sum);
                approx_vel_pc(v_1, num_particles, v, a, nodes[i], aux_b, dt, v_err_comp_sum);

                // Evaluate force function and store result
                error_status = WRAP_TRACEBACK(acceleration(
                    &aux_a[i * num_particles * 3],
                    &temp_system,
                    acceleration_param
                ));
                if (error_status.return_code != GRAV_SUCCESS)
                {
                    goto err_acc_error;
                }

                compute_aux_g(aux_g, num_particles, dim_nodes, aux_r, aux_a, i, F);
                compute_aux_b(aux_b, num_particles, dim_nodes_minus_1, aux_g, aux_c, i);
            }

            // Estimate convergence
            for (int i = 0; i < num_particles; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    delta_b7[i * 3 + j] = aux_b[dim_nodes_minus_2 * num_particles * 3 + i * 3 + j] - aux_b0[dim_nodes_minus_2 * num_particles * 3 + i * 3 + j];
                }
            }
            memcpy(aux_b0, aux_b, dim_nodes_minus_1 * num_particles * 3 * sizeof(double));
            if ((abs_max_vec(delta_b7, num_particles * 3) / abs_max_vec(&aux_a[dim_nodes_minus_1 * num_particles * 3], num_particles * 3)) < tolerance_pc)
            {
                break;
            }
        }
        
        /* Advance step */
        memcpy(temp_x_err_comp_sum, x_err_comp_sum, num_particles * 3 * sizeof(double));
        memcpy(temp_v_err_comp_sum, v_err_comp_sum, num_particles * 3 * sizeof(double));

        approx_pos_step(x_1, temp_x_err_comp_sum, num_particles, x, v, a, aux_b, dt);
        approx_vel_step(v_1, temp_v_err_comp_sum, num_particles, v, a, aux_b, dt);
        error_status = WRAP_TRACEBACK(acceleration(
            a_1,
            &temp_system,
            acceleration_param
        ));
        if (error_status.return_code != GRAV_SUCCESS)
        {
            goto err_acc_error;
        }

        /* Estimate relative error */
        error_b7 = abs_max_vec(&aux_b[dim_nodes_minus_2 * num_particles * 3], num_particles * 3) / abs_max_vec(a_1, num_particles * 3);
        error = pow((error_b7 / tolerance), exponent);

        // Prevent error from being too small
        if (error < 1e-10)
        {
            error = 1e-10;
        }
        dt_new = dt / error;

        /* Check error to accept the step */
        if (error <= 1.0 || dt == tf * 1e-12)
        {
            accept_step_flag = true;
            (*num_steps_ptr)++;
            *t_ptr += dt;

            refine_aux_b(aux_b, aux_e, delta_aux_b, num_particles, dim_nodes_minus_1, dt, dt_new, refine_flag);
            refine_flag = true;

            memcpy(x_err_comp_sum, temp_x_err_comp_sum, num_particles * 3 * sizeof(double));
            memcpy(v_err_comp_sum, temp_v_err_comp_sum, num_particles * 3 * sizeof(double));
        }

        if (accept_step_flag)
        {
            accept_step_flag = false;
            memcpy(x, x_1, num_particles * 3 * sizeof(double));
            memcpy(v, v_1, num_particles * 3 * sizeof(double));
            memcpy(a, a_1, num_particles * 3 * sizeof(double));    
        
            /* Output */
            if (is_output && *t_ptr >= next_output_time)
            {
                error_status = WRAP_TRACEBACK(output_snapshot(
                    output_param,
                    system,
                    integrator_param,
                    acceleration_param,
                    simulation_status,
                    settings
                ));
                if (error_status.return_code != GRAV_SUCCESS)
                {
                    goto err_output;
                }

                next_output_time = (*output_count_ptr) * output_interval;
            }
        }

        if (enable_progress_bar)
        {
            update_progress_bar(&progress_bar_param, *t_ptr, false);
        }

        /* Actual step size for the next step */
        if (dt_new > (dt / safety_fac))
        {
            dt /= safety_fac;
        }
        else if (dt_new < dt * safety_fac)
        {
            dt *= safety_fac;
        }
        else
        {
            dt = dt_new;
        }

        if (dt_new < tf * 1e-12)
        {
            dt = tf * 1e-12;
        }

        // Correct overshooting
        if (((*t_ptr) < tf) && ((*t_ptr) + dt > tf))
        {
            dt = tf - (*t_ptr);
        }
        simulation_status->dt = dt;

        /* Check exit */
        if (*(settings->is_exit_ptr))
        {
            break;
        }
    }

    if (enable_progress_bar)
    {
        update_progress_bar(&progress_bar_param, *t_ptr, true);
    }

    /* Free memory */
    free(nodes);
    free(aux_c);
    free(aux_r);
    free(aux_b0);
    free(aux_b);
    free(aux_g);
    free(aux_e);
    free(a);
    free(aux_a);
    free(x_1);
    free(v_1);
    free(a_1);
    free(delta_b7);
    free(F);
    free(delta_aux_b);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    free(temp_x_err_comp_sum);
    free(temp_v_err_comp_sum);

    return make_success_error_status();

err_output:
err_acc_error:
err_start_progress_bar:
err_initial_output:
err_initial_dt:
err_memory:
    free(x_err_comp_sum);
    free(v_err_comp_sum);
    free(temp_x_err_comp_sum);
    free(temp_v_err_comp_sum);
    free(delta_aux_b);
    free(F);
    free(delta_b7);
    free(a_1);
    free(v_1);
    free(x_1);
    free(aux_a);
    free(a);
err_aux_memory:
    free(aux_e);
    free(aux_g);
    free(aux_b);
    free(aux_b0);
    free(aux_r);
    free(aux_c);
    free(nodes);

    return error_status;
}

IN_FILE void initialize_radau_spacing(double *restrict nodes)
{
    nodes[0] = 0.0L;
    nodes[1] = 0.056262560536922146465652191032L;
    nodes[2] = 0.180240691736892364987579942809L;
    nodes[3] = 0.352624717113169637373907770171L;
    nodes[4] = 0.547153626330555383001448557652L;
    nodes[5] = 0.734210177215410531523210608306L;
    nodes[6] = 0.885320946839095768090359762932L;
    nodes[7] = 0.977520613561287501891174500429L;
}

IN_FILE void initialize_aux_c(double *restrict aux_c)
{
    for (int i = 0; i < 7; i++)
    {
        aux_c[i * 7 + i] = 1.0L;
    }

    aux_c[1 * 7 + 0] = -0.0562625605369221464656522L;

    aux_c[2 * 7 + 0] = 0.01014080283006362998648180399549641417413495311078L;
    aux_c[2 * 7 + 1] = -0.2365032522738145114532321L;

    aux_c[3 * 7 + 0] = -0.0035758977292516175949344589284567187362040464593728L;
    aux_c[3 * 7 + 1] = 0.09353769525946206589574845561035371499343547051116L;
    aux_c[3 * 7 + 2] = -0.5891279693869841488271399L;

    aux_c[4 * 7 + 0] = 0.0019565654099472210769005672379668610648179838140913L;
    aux_c[4 * 7 + 1] = -0.054755386889068686440808430671055022602028382584495L;
    aux_c[4 * 7 + 2] = 0.41588120008230686168862193041156933067050816537030L;
    aux_c[4 * 7 + 3] = -1.1362815957175395318285885L;

    aux_c[5 * 7 + 0] = -0.0014365302363708915424459554194153247134438571962198L;
    aux_c[5 * 7 + 1] = 0.042158527721268707707297347813203202980228135395858L;
    aux_c[5 * 7 + 2] = -0.36009959650205681228976647408968845289781580280782L;
    aux_c[5 * 7 + 3] = 1.2501507118406910258505441186857527694077565516084L;
    aux_c[5 * 7 + 4] = -1.8704917729329500633517991L;

    aux_c[6 * 7 + 0] = 0.0012717903090268677492943117622964220889484666147501L;
    aux_c[6 * 7 + 1] = -0.038760357915906770369904626849901899108502158354383L;
    aux_c[6 * 7 + 2] = 0.36096224345284598322533983078129066420907893718190L;
    aux_c[6 * 7 + 3] = -1.4668842084004269643701553461378480148761655599754L;
    aux_c[6 * 7 + 4] = 2.9061362593084293014237914371173946705384212479246L;
    aux_c[6 * 7 + 5] = -2.7558127197720458314421589L;
}

void initialize_aux_r(double *restrict aux_r)
{
    aux_r[1 * 8 + 0] = 17.773808914078000840752659565672904106978971632681L;
    aux_r[2 * 8 + 0] = 5.5481367185372165056928216140765061758579336941398L;
    aux_r[3 * 8 + 0] = 2.8358760786444386782520104428042437400879003147949L;
    aux_r[4 * 8 + 0] = 1.8276402675175978297946077587371204385651628457154L;
    aux_r[5 * 8 + 0] = 1.3620078160624694969370006292445650994197371928318L;
    aux_r[6 * 8 + 0] = 1.1295338753367899027322861542728593509768148769105L;
    aux_r[7 * 8 + 0] = 1.0229963298234867458386119071939636779024159134103L;

    aux_r[2 * 8 + 1] = 8.0659386483818866885371256689687154412267416180207L;
    aux_r[3 * 8 + 1] = 3.3742499769626352599420358188267460448330087696743L;
    aux_r[4 * 8 + 1] = 2.0371118353585847827949159161566554921841792590404L;
    aux_r[5 * 8 + 1] = 1.4750402175604115479218482480167404024740127431358L;
    aux_r[6 * 8 + 1] = 1.2061876660584456166252036299646227791474203527801L;
    aux_r[7 * 8 + 1] = 1.0854721939386423840467243172568913862030118679827L;

    aux_r[3 * 8 + 2] = 5.8010015592640614823286778893918880155743979164251L;
    aux_r[4 * 8 + 2] = 2.7254422118082262837742722003491334729711450288807L;
    aux_r[5 * 8 + 2] = 1.8051535801402512604391147435448679586574414080693L;
    aux_r[6 * 8 + 2] = 1.4182782637347391537713783674858328433713640692518L;
    aux_r[7 * 8 + 2] = 1.2542646222818777659905422465868249586862369725826L;

    aux_r[4 * 8 + 3] = 5.1406241058109342286363199091504437929335189668304L;
    aux_r[5 * 8 + 3] = 2.6206449263870350811541816031933074696730227729812L;
    aux_r[6 * 8 + 3] = 1.8772424961868100972169920283109658335427446084411L;
    aux_r[7 * 8 + 3] = 1.6002665494908162609916716949161150366323259154408L;

    aux_r[5 * 8 + 4] = 5.3459768998711075141214909632277898045770336660354L;
    aux_r[6 * 8 + 4] = 2.9571160172904557478071040204245556508352776929762L;
    aux_r[7 * 8 + 4] = 2.3235983002196942228325345451091668073608955835034L;

    aux_r[6 * 8 + 5] = 6.6176620137024244874471284891193925737033291491748L;
    aux_r[7 * 8 + 5] = 4.1099757783445590862385761824068782144723082633980L;

    aux_r[7 * 8 + 6] = 10.846026190236844684706431007823415424143683137181L;
}

IN_FILE ErrorStatus ias15_initial_dt(
    double *restrict initial_dt,
    const int power,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param,
    const double *restrict a
)
{
    ErrorStatus error_status;
    *initial_dt = -1.0;

    const int num_particles = system->num_particles;
    double *restrict x = system->x;
    double *restrict v = system->v;

    /* Allocate memory and declare variables */
    double *restrict x_1 = malloc(num_particles * 3 * sizeof(double));
    double *restrict a_1 = malloc(num_particles * 3 * sizeof(double));
    if (!x_1 || !a_1)
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for initial_dt calculation"
        );
        goto error_memory;
    }

    double d_0 = abs_max_vec(x, num_particles * 3);
    double d_1 = abs_max_vec(a, num_particles * 3);
    double d_2;
    double dt_0;
    double dt_1;
    System system_1 = {
        .num_particles = num_particles,
        .x = x_1,
        .v = v,
        .m = system->m,
        .G = system->G,
    };

    if (d_0 < 1e-5 || d_1 < 1e-5)
    {
        dt_0 = 1e-6;
    }
    else
    {
        dt_0 = 0.01 * (d_0 / d_1);
    }
    
    for (int i = 0; i < num_particles; i++)
    {
        x_1[i * 3 + 0] = x[i * 3 + 0] + dt_0 * v[i * 3 + 0];
        x_1[i * 3 + 1] = x[i * 3 + 1] + dt_0 * v[i * 3 + 1];
        x_1[i * 3 + 2] = x[i * 3 + 2] + dt_0 * v[i * 3 + 2];
    }
    
    error_status = WRAP_TRACEBACK(acceleration(
        a_1,
        &system_1,
        acceleration_param
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto error_acc;
    }

    for (int i = 0; i < num_particles; i++)
    {
        a_1[i * 3 + 0] = a_1[i * 3 + 0] - a[i * 3 + 0];
        a_1[i * 3 + 1] = a_1[i * 3 + 1] - a[i * 3 + 1];
        a_1[i * 3 + 2] = a_1[i * 3 + 2] - a[i * 3 + 2];
    }
    d_2 = abs_max_vec(a_1, num_particles * 3) / dt_0;

    if (fmax(d_1, d_2) <= 1e-15)
    {
        dt_1 = fmax(1e-6, dt_0 * 1e-3);
    }
    else
    {
        dt_1 = pow((0.01 / fmax(d_1, d_2)), 1.0 / (1.0 + (double) power));
    }

    *initial_dt = fmin(100.0 * dt_0, dt_1);
    if (*initial_dt <= 0.0)
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_VALUE_ERROR,
            "Initial dt is less than or equal to zero"
        );
        goto error_dt;
    }

    free(x_1);
    free(a_1);

    return make_success_error_status();

error_dt:
error_acc:
error_memory:
    free(x_1);
    free(a_1);
    return error_status;
}

IN_FILE void approx_pos_pc(
    double *restrict x,
    const int num_particles,
    const double *restrict x0,
    const double *restrict v0,
    const double *restrict a0,
    const double node,
    const double *restrict aux_b,
    const double dt,
    const double *restrict x_err_comp_sum
)
{   
    for (int j = 0; j < num_particles; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            /*
            *   Warning: Because x0 is much larger than dx, combining both statements 
            *            would increase floating point error
            *            e.g. x[j][k] = x0[j][k] + x_err_comp_sum[j][k] + ...   (WRONG)
            */

            x[j * 3 + k] = x0[j * 3 + k];     
            x[j * 3 + k] += x_err_comp_sum[j * 3 + k] + dt * node * (
                v0[j * 3 + k]
                + dt
                * node
                * (
                    a0[j * 3 + k]
                    + node 
                    * (
                        aux_b[0 * num_particles * 3 + j * 3 + k] / 3.0
                        + node
                        * (
                            aux_b[1 * num_particles * 3 + j * 3 + k] / 6.0
                            + node
                            * (
                                aux_b[2 * num_particles * 3 + j * 3 + k] / 10.0
                                + node
                                * (
                                    aux_b[3 * num_particles * 3 + j * 3 + k] / 15.0
                                    + node
                                    * (
                                        aux_b[4 * num_particles * 3 + j * 3 + k] / 21.0
                                        + node * (
                                            aux_b[5 * num_particles * 3 + j * 3 + k] / 28.0 
                                            + node * aux_b[6 * num_particles * 3 + j * 3 + k] / 36.0
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
                / 2.0
            );
        }
    }
}

IN_FILE void approx_vel_pc(
    double *restrict v,
    const int num_particles,
    const double *restrict v0,
    const double *restrict a0,
    const double node,
    const double *restrict aux_b,
    const double dt,
    const double *restrict v_err_comp_sum
)
{
    for (int j = 0; j < num_particles; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            /*
            *   Warning: Because v0 is much larger than dv, combining both statements 
            *            would increase floating point error
            *            e.g. v[j][k] = v0[j][k] + v_err_comp_sum[j][k] + ...   (WRONG)
            */

            v[j * 3 + k] = v0[j * 3 + k];
            v[j * 3 + k] += v_err_comp_sum[j * 3 + k] + dt * node * (
                a0[j * 3 + k]
                + node
                * (
                    aux_b[0 * num_particles * 3 + j * 3 + k] / 2.0
                    + node
                    * (
                        aux_b[1 * num_particles * 3 + j * 3 + k] / 3.0
                        + node
                        * (
                            aux_b[2 * num_particles * 3 + j * 3 + k] / 4.0
                            + node
                            * (
                                aux_b[3 * num_particles * 3 + j * 3 + k] / 5.0
                                + node
                                * (
                                    aux_b[4 * num_particles * 3 + j * 3 + k] / 6.0
                                    + node * (
                                        aux_b[5 * num_particles * 3 + j * 3 + k] / 7.0 
                                        + node * aux_b[6 * num_particles * 3 + j * 3 + k] / 8.0
                                    )
                                )
                            )
                        )
                    )
                )
            );
        }
    }
}

IN_FILE void approx_pos_step(
    double *restrict x,
    double *restrict temp_x_err_comp_sum,
    const int num_particles,
    const double *restrict x0,
    const double *restrict v0,
    const double *restrict a0,
    const double *restrict aux_b,
    const double dt
)
{   
    for (int j = 0; j < num_particles; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            temp_x_err_comp_sum[j * 3 + k] += dt * (
                v0[j * 3 + k] + dt * (a0[j * 3 + k]
                    + aux_b[0 * num_particles * 3 + j * 3 + k] / 3.0
                    + aux_b[1 * num_particles * 3 + j * 3 + k] / 6.0
                    + aux_b[2 * num_particles * 3 + j * 3 + k] / 10.0
                    + aux_b[3 * num_particles * 3 + j * 3 + k] / 15.0
                    + aux_b[4 * num_particles * 3 + j * 3 + k] / 21.0
                    + aux_b[5 * num_particles * 3 + j * 3 + k] / 28.0 
                    + aux_b[6 * num_particles * 3 + j * 3 + k] / 36.0
                )
                / 2.0
            );

            x[j * 3 + k] = x0[j * 3 + k] + temp_x_err_comp_sum[j * 3 + k];
            temp_x_err_comp_sum[j * 3 + k] += (x0[j * 3 + k] - x[j * 3 + k]);
        }
    }
}

IN_FILE void approx_vel_step(
    double *restrict v,
    double *restrict temp_v_err_comp_sum,
    const int num_particles,
    const double *restrict v0,
    const double *restrict a0,
    const double *restrict aux_b,
    const double dt
)
{
    for (int j = 0; j < num_particles; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            temp_v_err_comp_sum[j * 3 + k] += dt * (
                a0[j * 3 + k]
                + aux_b[0 * num_particles * 3 + j * 3 + k] / 2.0
                + aux_b[1 * num_particles * 3 + j * 3 + k] / 3.0
                + aux_b[2 * num_particles * 3 + j * 3 + k] / 4.0
                + aux_b[3 * num_particles * 3 + j * 3 + k] / 5.0
                + aux_b[4 * num_particles * 3 + j * 3 + k] / 6.0
                + aux_b[5 * num_particles * 3 + j * 3 + k] / 7.0 
                + aux_b[6 * num_particles * 3 + j * 3 + k] / 8.0
            );
            v[j * 3 + k] = v0[j * 3 + k] + temp_v_err_comp_sum[j * 3 + k];
            temp_v_err_comp_sum[j * 3 + k] += (v0[j * 3 + k] - v[j * 3 + k]);
        }
    }
}

IN_FILE void compute_aux_b(
    double *restrict aux_b,
    const int num_particles,
    const int dim_nodes_minus_1,
    const double *restrict aux_g,
    const double *restrict aux_c,
    const int i
)
{
    for (int j = 0; j < num_particles; j++)
    {
        if (i >= 1) 
        {
            for (int k = 0; k < 3; k++)
            {
                aux_b[0 * num_particles * 3 + j * 3 + k] = (
                    aux_c[0 * dim_nodes_minus_1 + 0] * aux_g[0 * num_particles * 3 + j * 3 + k]
                    + aux_c[1 * dim_nodes_minus_1 + 0] * aux_g[1 * num_particles * 3 + j * 3 + k]
                    + aux_c[2 * dim_nodes_minus_1 + 0] * aux_g[2 * num_particles * 3 + j * 3 + k]
                    + aux_c[3 * dim_nodes_minus_1 + 0] * aux_g[3 * num_particles * 3 + j * 3 + k]
                    + aux_c[4 * dim_nodes_minus_1 + 0] * aux_g[4 * num_particles * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 0] * aux_g[5 * num_particles * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 0] * aux_g[6 * num_particles * 3 + j * 3 + k]
                );
            }
        }
        else
        {
            continue;
        }

        if (i >= 2) 
        {
            for (int k = 0; k < 3; k++)
            {
                aux_b[1 * num_particles * 3 + j * 3 + k] = (
                    aux_c[1 * dim_nodes_minus_1 + 1] * aux_g[1 * num_particles * 3 + j * 3 + k]
                    + aux_c[2 * dim_nodes_minus_1 + 1] * aux_g[2 * num_particles * 3 + j * 3 + k]
                    + aux_c[3 * dim_nodes_minus_1 + 1] * aux_g[3 * num_particles * 3 + j * 3 + k]
                    + aux_c[4 * dim_nodes_minus_1 + 1] * aux_g[4 * num_particles * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 1] * aux_g[5 * num_particles * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 1] * aux_g[6 * num_particles * 3 + j * 3 + k]
                );
            }
        }
        else
        {
            continue;
        }

        if (i >= 3)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_b[2 * num_particles * 3 + j * 3 + k] = (
                    aux_c[2 * dim_nodes_minus_1 + 2] * aux_g[2 * num_particles * 3 + j * 3 + k]
                    + aux_c[3 * dim_nodes_minus_1 + 2] * aux_g[3 * num_particles * 3 + j * 3 + k]
                    + aux_c[4 * dim_nodes_minus_1 + 2] * aux_g[4 * num_particles * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 2] * aux_g[5 * num_particles * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 2] * aux_g[6 * num_particles * 3 + j * 3 + k]
                );
            }
        }
        else
        {
            continue;
        }

        if (i >= 4)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_b[3 * num_particles * 3 + j * 3 + k] = (
                    aux_c[3 * dim_nodes_minus_1 + 3] * aux_g[3 * num_particles * 3 + j * 3 + k]
                    + aux_c[4 * dim_nodes_minus_1 + 3] * aux_g[4 * num_particles * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 3] * aux_g[5 * num_particles * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 3] * aux_g[6 * num_particles * 3 + j * 3 + k]
                );
            }
        }
        else
        {
            continue;
        }

        if (i >= 5)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_b[4 * num_particles * 3 + j * 3 + k] = (
                    aux_c[4 * dim_nodes_minus_1 + 4] * aux_g[4 * num_particles * 3 + j * 3 + k]
                    + aux_c[5 * dim_nodes_minus_1 + 4] * aux_g[5 * num_particles * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 4] * aux_g[6 * num_particles * 3 + j * 3 + k]
                );
            }
        }
        else
        {
            continue;
        }

        if (i >= 6)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_b[5 * num_particles * 3 + j * 3 + k] = (
                    aux_c[5 * dim_nodes_minus_1 + 5] * aux_g[5 * num_particles * 3 + j * 3 + k]
                    + aux_c[6 * dim_nodes_minus_1 + 5] * aux_g[6 * num_particles * 3 + j * 3 + k]
                );
            }
        }
        else
        {
            continue;
        }

        if (i >= 7)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_b[6 * num_particles * 3 + j * 3 + k] = (
                    aux_c[6 * dim_nodes_minus_1 + 6] * aux_g[6 * num_particles * 3 + j * 3 + k]
                );
            }
        }
        else
        {
            continue;
        }
    }
}

IN_FILE void compute_aux_g(
    double *restrict aux_g,
    const int num_particles,
    const int dim_nodes,
    const double *restrict aux_r,
    const double *restrict aux_a,
    const int i,
    double *restrict F
)
{
    // Retrieve required accelerations
    // F is allocated in IAS15
    memcpy(F, aux_a, (i + 1) * num_particles * 3 * sizeof(double));

    // Update aux_g
    for (int j = 0; j < num_particles; j++)
    {
        if (i >= 1)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_g[0 * num_particles * 3 + j * 3 + k] = (
                    (F[1 * num_particles * 3 + j * 3 + k] - F[0 * num_particles * 3 + j * 3 + k]) * aux_r[1 * dim_nodes + 0]
                );
            }
        }
        else
        {
            continue;
        }
                
        if (i >= 2)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_g[1 * num_particles * 3 + j * 3 + k] = (
                    ((F[2 * num_particles * 3 + j * 3 + k] - F[0 * num_particles * 3 + j * 3 + k]) * aux_r[2 * dim_nodes + 0] 
                    - aux_g[0 * num_particles * 3 + j * 3 + k]) * aux_r[2 * dim_nodes + 1]
                );
            }
        }
        else 
        {
            continue;
        }

        if (i >= 3)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_g[2 * num_particles * 3 + j * 3 + k] = (
                    ((F[3 * num_particles * 3 + j * 3 + k] - F[0 * num_particles * 3 + j * 3 + k]) * aux_r[3 * dim_nodes + 0] 
                    - aux_g[0 * num_particles * 3 + j * 3 + k]) * aux_r[3 * dim_nodes + 1] 
                    - aux_g[1 * num_particles * 3 + j * 3 + k]
                ) * aux_r[3 * dim_nodes + 2];
            }
        }
        else
        {
            continue;
        }
                
        if (i >= 4)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_g[3 * num_particles * 3 + j * 3 + k] = (
                    (((F[4 * num_particles * 3 + j * 3 + k] - F[0 * num_particles * 3 + j * 3 + k]) * aux_r[4 * dim_nodes + 0] 
                    - aux_g[0 * num_particles * 3 + j * 3 + k]) * aux_r[4 * dim_nodes + 1] - aux_g[1 * num_particles * 3 + j * 3 + k])
                    * aux_r[4 * dim_nodes + 2]
                    - aux_g[2 * num_particles * 3 + j * 3 + k]
                ) * aux_r[4 * dim_nodes + 3];
            }
        }
        else
        {
            continue;
        }

        if (i >= 5)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_g[4 * num_particles * 3 + j * 3 + k] = (
                    (
                        (((F[5 * num_particles * 3 + j * 3 + k] - F[0 * num_particles * 3 + j * 3 + k]) * aux_r[5 * dim_nodes + 0] 
                        - aux_g[0 * num_particles * 3 + j * 3 + k]) * aux_r[5 * dim_nodes + 1] 
                        - aux_g[1 * num_particles * 3 + j * 3 + k])
                        * aux_r[5 * dim_nodes + 2]
                        - aux_g[2 * num_particles * 3 + j * 3 + k]
                    )
                    * aux_r[5 * dim_nodes + 3]
                    - aux_g[3 * num_particles * 3 + j * 3 + k]
                ) * aux_r[5 * dim_nodes + 4];
            }
        }
        else
        {
            continue;
        }

        if (i >= 6)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_g[5 * num_particles * 3 + j * 3 + k] = (
                    (
                        (
                            (((F[6 * num_particles * 3 + j * 3 + k] - F[0 * num_particles * 3 + j * 3 + k]) * aux_r[6 * dim_nodes + 0] 
                            - aux_g[0 * num_particles * 3 + j * 3 + k]) * aux_r[6 * dim_nodes + 1] 
                            - aux_g[1 * num_particles * 3 + j * 3 + k])
                            * aux_r[6 * dim_nodes + 2]
                            - aux_g[2 * num_particles * 3 + j * 3 + k]
                        )
                        * aux_r[6 * dim_nodes + 3]
                        - aux_g[3 * num_particles * 3 + j * 3 + k]
                    )
                    * aux_r[6 * dim_nodes + 4]
                    - aux_g[4 * num_particles * 3 + j * 3 + k]
                ) * aux_r[6 * dim_nodes + 5];
            }
        }
        else
        {
            continue;
        }

        if (i >= 7)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_g[6 * num_particles * 3 + j * 3 + k] = (
                    (
                        (
                            (
                                (((F[7 * num_particles * 3 + j * 3 + k] - F[0 * num_particles * 3 + j * 3 + k]) * aux_r[7 * dim_nodes + 0] - aux_g[0 * num_particles * 3 + j * 3 + k]) 
                                * aux_r[7 * dim_nodes + 1] 
                                - aux_g[1 * num_particles * 3 + j * 3 + k])
                                * aux_r[7 * dim_nodes + 2]
                                - aux_g[2 * num_particles * 3 + j * 3 + k]
                            )
                            * aux_r[7 * dim_nodes + 3]
                            - aux_g[3 * num_particles * 3 + j * 3 + k]
                        )
                        * aux_r[7 * dim_nodes + 4]
                        - aux_g[4 * num_particles * 3 + j * 3 + k]
                    )
                    * aux_r[7 * dim_nodes + 5]
                    - aux_g[5 * num_particles * 3 + j * 3 + k]
                ) * aux_r[7 * dim_nodes + 6];                
            }
        }
        else
        {
            continue;
        }
    }
}

IN_FILE void refine_aux_b(
    double *restrict aux_b,
    double *restrict aux_e,
    double *restrict delta_aux_b,
    const int num_particles,
    const int dim_nodes_minus_1,
    const double dt,
    const double dt_new,
    const bool refine_flag
)
{
    if (refine_flag)
    {
        for (int i = 0; i < dim_nodes_minus_1; i++)
        {
            for (int j = 0; j < num_particles; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    delta_aux_b[i * num_particles * 3 + j * 3 + k] = (
                        aux_b[i * num_particles * 3 + j * 3 + k] - aux_e[i * num_particles * 3 + j * 3 + k]
                    );
                }
            }
        }
    }
    else
    {
        // Empty delta_aux_b
        for (int i = 0; i < dim_nodes_minus_1 * num_particles; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                delta_aux_b[i * 3 + j] = 0.0;
            }
        }
    }

    double q = dt_new / dt;
    double q2 = q * q, q3 = q2 * q, q4 = q3 * q, q5 = q4 * q, q6 = q5 * q, q7 = q6 * q;
    
    for (int j = 0; j < num_particles; j++) 
    {
        for (int k = 0; k < 3; k++) 
        {
            aux_e[0 * num_particles * 3 + j * 3 + k] = q * (
                aux_b[6 * num_particles * 3 + j * 3 + k] * 7.0
                + aux_b[5 * num_particles * 3 + j * 3 + k] * 6.0
                + aux_b[4 * num_particles * 3 + j * 3 + k] * 5.0
                + aux_b[3 * num_particles * 3 + j * 3 + k] * 4.0
                + aux_b[2 * num_particles * 3 + j * 3 + k] * 3.0
                + aux_b[1 * num_particles * 3 + j * 3 + k] * 2.0
                + aux_b[0 * num_particles * 3 + j * 3 + k]
            );

            aux_e[1 * num_particles * 3 + j * 3 + k] = q2 * (
                aux_b[6 * num_particles * 3 + j * 3 + k] * 21.0
                + aux_b[5 * num_particles * 3 + j * 3 + k] * 15.0
                + aux_b[4 * num_particles * 3 + j * 3 + k] * 10.0
                + aux_b[3 * num_particles * 3 + j * 3 + k] * 6.0
                + aux_b[2 * num_particles * 3 + j * 3 + k] * 3.0
                + aux_b[1 * num_particles * 3 + j * 3 + k]
            );

            aux_e[2 * num_particles * 3 + j * 3 + k] = q3 * (
                aux_b[6 * num_particles * 3 + j * 3 + k] * 35.0
                + aux_b[5 * num_particles * 3 + j * 3 + k] * 20.0
                + aux_b[4 * num_particles * 3 + j * 3 + k] * 10.0
                + aux_b[3 * num_particles * 3 + j * 3 + k] * 4.0
                + aux_b[2 * num_particles * 3 + j * 3 + k]
            );

            aux_e[3 * num_particles * 3 + j * 3 + k] = q4 * (
                aux_b[6 * num_particles * 3 + j * 3 + k] * 35.0
                + aux_b[5 * num_particles * 3 + j * 3 + k] * 15.0
                + aux_b[4 * num_particles * 3 + j * 3 + k] * 5.0
                + aux_b[3 * num_particles * 3 + j * 3 + k]
            );

            aux_e[4 * num_particles * 3 + j * 3 + k] = q5 * (
                aux_b[6 * num_particles * 3 + j * 3 + k] * 21.0
                + aux_b[5 * num_particles * 3 + j * 3 + k] * 6.0
                + aux_b[4 * num_particles * 3 + j * 3 + k]
            );

            aux_e[5 * num_particles * 3 + j * 3 + k] = q6 * (
                aux_b[6 * num_particles * 3 + j * 3 + k] * 7.0
                + aux_b[5 * num_particles * 3 + j * 3 + k]
            );

            aux_e[6 * num_particles * 3 + j * 3 + k] = q7 * aux_b[6 * num_particles * 3 + j * 3 + k];
        }
    }

    for (int i = 0; i < dim_nodes_minus_1; i++)
    {
        for (int j = 0; j < num_particles; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                aux_b[i * num_particles * 3 + j * 3 + k] = (
                    aux_e[i * num_particles * 3 + j * 3 + k] + delta_aux_b[i * num_particles * 3 + j * 3 + k]
                );
            }
        }
    }
}
