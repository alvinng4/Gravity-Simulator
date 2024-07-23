#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"


/**
 * \brief Euler integrator
 * 
 * \param system Name of the system
 * \param dt Time step of the system
 * \param npts Number of time steps to be integrated
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param custom_sys_x Array of the position vectors of customized system 
 * \param custom_sys_v Array of the velocity vectors of customized system 
 * \param custom_sys_m Array of the masses of customized system 
 * \param custom_sys_G Gravitational constant of customized system
 * \param custom_sys_object_count Number of objects in customized system
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int euler(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count,
    Solutions *restrict solution,
    int *restrict is_exit
)
{   
    // Initialize default system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int intitial_sys_flag = initialize_system(system, &x, &v, &m, &objects_count, &G);

    // Custom system
    if (intitial_sys_flag == 2)
    {
        printf("Error: Failed to allocate memory to initialize default system\n");
        goto err_default_sys_memory;
    }
    else if (intitial_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

        if (!x || !v || !m)
        {
            printf("Error: Failed to allocate memory for customized system\n");
            goto err_custom_sys_memory;
        }

        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = custom_sys_x[i * 3 + j];
                v[i * 3 + j] = custom_sys_v[i * 3 + j];
            }
            m[i] = custom_sys_m[i];
        }
        G = custom_sys_G;
    }
    
    // Allocate memory for calculation
    real *temp_x = malloc(objects_count * 3 * sizeof(real));
    real *temp_v = malloc(objects_count * 3 * sizeof(real));
    real *a = malloc(objects_count * 3 * sizeof(real));
    
    // Allocate memory for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (!temp_x || !temp_v || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        printf("Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(store_npts * sizeof(double));
    double *sol_dt = malloc(store_npts * sizeof(double));

    if (!sol_state || !sol_time || !sol_dt)
    {
        printf("Error: Failed to allocate memory for solution output\n");
        goto err_sol_output_memory;
    }

    // Initial value
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sol_state[i * 3 + j] = x[i * 3 + j];
            sol_state[objects_count * 3 + i * 3 + j] = v[i * 3 + j];
        }
    }
    sol_time[0] = 0.0;
    for (int i = 0; i < store_npts; i++)
    {
        sol_dt[i] = dt;
    }

    // Main Loop
    while ((count + 1) <= npts)
    {
        acceleration(objects_count, x, a, m, G);

        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));

        // Calculation
        for (int i = 0; i < objects_count; i++) 
        {
            for (int j = 0; j < 3; j++) 
            {
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;

                x[i * 3 + j] = temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];

                x_err_comp_sum[i * 3 + j] += temp_x[i * 3 + j] - x[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];
            }
        }

        // Store solution
        if ((count + 1) == npts)
        {
            memcpy(&sol_state[(store_npts - 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
            memcpy(&sol_state[(store_npts - 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            sol_time[store_npts - 1] = dt * count;
        }
        else if (((count + 1) % store_every_n == 0))
        {
            memcpy(&sol_state[(*store_count + 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
            memcpy(&sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            sol_time[*store_count + 1] = dt * (count + 1);
            if (!((*store_count + 1) == (store_npts - 1)))
            {
                (*store_count)++;
            }
        }
        count++;

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }
    }

    // Exit after simulation is finished
    free(x);
    free(v);
    free(a);
    free(temp_x);
    free(temp_v);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    solution->sol_state = sol_state;
    solution->sol_time = sol_time;
    solution->sol_dt = sol_dt;
    solution->m = m;
    solution->G = G;
    solution->objects_count = objects_count;

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_sol_output_memory:
    free(sol_state);
    free(sol_time);
    free(sol_dt);
err_calc_memory:
    free(temp_x);
    free(temp_v);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
err_custom_sys_memory:
    free(x);
    free(v);
    free(m);
err_default_sys_memory:
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}

/**
 * \brief Euler-cromer integrator
 * 
 * \param system Name of the system
 * \param dt Time step of the system
 * \param npts Number of time steps to be integrated
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param custom_sys_x Array of the position vectors of customized system 
 * \param custom_sys_v Array of the velocity vectors of customized system 
 * \param custom_sys_m Array of the masses of customized system 
 * \param custom_sys_G Gravitational constant of customized system
 * \param custom_sys_object_count Number of objects in customized system
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int euler_cromer(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count,
    Solutions *restrict solution,
    int *restrict is_exit
)
{   
    // Initialize default system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int intitial_sys_flag = initialize_system(system, &x, &v, &m, &objects_count, &G);

    // Custom system
    if (intitial_sys_flag == 2)
    {
        printf("Error: Failed to allocate memory to initialize default system\n");
        goto err_default_sys_memory;
    }
    else if (intitial_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

        if (!x || !v || !m)
        {
            printf("Error: Failed to allocate memory for customized system\n");
            goto err_custom_sys_memory;
        }

        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = custom_sys_x[i * 3 + j];
                v[i * 3 + j] = custom_sys_v[i * 3 + j];
            }
            m[i] = custom_sys_m[i];
        }
        G = custom_sys_G;
    }

    // Allocate memory for calculation
    real *temp_x = malloc(objects_count * 3 * sizeof(real));
    real *temp_v = malloc(objects_count * 3 * sizeof(real));
    real *a = malloc(objects_count * 3 * sizeof(real));
    
    // Allocate memory for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    
    if (!temp_x || !temp_v || !a || !x_err_comp_sum || !v_err_comp_sum)
    {
        printf("Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(store_npts * sizeof(double));
    double *sol_dt = malloc(store_npts * sizeof(double));

    if (!sol_state || !sol_time || !sol_dt)
    {
        printf("Error: Failed to allocate memory for solution output\n");
        goto err_sol_output_memory;
    }

    // Initial value
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sol_state[i * 3 + j] = x[i * 3 + j];
            sol_state[objects_count * 3 + i * 3 + j] = v[i * 3 + j];
        }
    }
    sol_time[0] = 0.0;
    for (int i = 0; i < store_npts; i++)
    {
        sol_dt[i] = dt;
    }

    // Main Loop
    while ((count + 1) <= npts)
    {
        acceleration(objects_count, x, a, m, G);

        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                // Calculate v
                v_err_comp_sum[i * 3 + j] += a[i * 3 + j] * dt;
                v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];

                // Calculate x
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt;
                x[i * 3 + j] = temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += temp_x[i * 3 + j] - x[i * 3 + j];
            }
        }

        // Store solution
        if ((count + 1) == npts)
        {
            memcpy(&sol_state[(store_npts - 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
            memcpy(&sol_state[(store_npts - 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            sol_time[store_npts - 1] = dt * count;
        }
        else if (((count + 1) % store_every_n == 0))
        {
            memcpy(&sol_state[(*store_count + 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
            memcpy(&sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            sol_time[*store_count + 1] = dt * (count + 1);
            if (!((*store_count + 1) == (store_npts - 1)))
            {
                (*store_count)++;
            }
        }
        count++;

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }
    }

    // Exit after simulation is finished
    free(x);
    free(v);
    free(a);
    free(temp_x);
    free(temp_v);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    solution->sol_state = sol_state;
    solution->sol_time = sol_time;
    solution->sol_dt = sol_dt;
    solution->m = m;
    solution->G = G;
    solution->objects_count = objects_count;

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_sol_output_memory:
    free(sol_state);
    free(sol_time);
    free(sol_dt);
err_calc_memory:
    free(temp_x);
    free(temp_v);
    free(a);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
err_custom_sys_memory:
    free(x);
    free(v);
    free(m);
err_default_sys_memory: 
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}

/**
 * \brief RK4 integrator
 * 
 * \param system Name of the system
 * \param dt Time step of the system
 * \param npts Number of time steps to be integrated
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param custom_sys_x Array of the position vectors of customized system 
 * \param custom_sys_v Array of the velocity vectors of customized system 
 * \param custom_sys_m Array of the masses of customized system 
 * \param custom_sys_G Gravitational constant of customized system
 * \param custom_sys_object_count Number of objects in customized system
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int rk4(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count,
    Solutions *restrict solution,
    int *restrict is_exit
)
{
    // Initialize default system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int intitial_sys_flag = initialize_system(system, &x, &v, &m, &objects_count, &G);

    // Custom system
    if (intitial_sys_flag == 2)
    {
        printf("Error: Failed to allocate memory to initialize default system\n");
        goto err_default_sys_memory;
    }
    else if (intitial_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

        if (!x || !v || !m)
        {
            printf("Error: Failed to allocate memory for customized system\n");
            goto err_custom_sys_memory;
        }

        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = custom_sys_x[i * 3 + j];
                v[i * 3 + j] = custom_sys_v[i * 3 + j];
            }
            m[i] = custom_sys_m[i];
        }
        G = custom_sys_G;
    }

    // Allocate memory for calculation
    real *temp_x = malloc(objects_count * 3 * sizeof(real));
    real *temp_v = malloc(objects_count * 3 * sizeof(real));
    real *a = malloc(objects_count * 3 * sizeof(real));
    real *vk1 = malloc(objects_count * 3 * sizeof(real));
    real *vk2 = malloc(objects_count * 3 * sizeof(real));
    real *vk3 = malloc(objects_count * 3 * sizeof(real));
    real *vk4 = malloc(objects_count * 3 * sizeof(real));
    real *xk1 = malloc(objects_count * 3 * sizeof(real));
    real *xk2 = malloc(objects_count * 3 * sizeof(real));
    real *xk3 = malloc(objects_count * 3 * sizeof(real));
    real *xk4 = malloc(objects_count * 3 * sizeof(real));

    // Allocate memory for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (!temp_x || !temp_v || !a || !vk1 || !vk2 || !vk3 || !vk4 || !xk1 || !xk2 || !xk3 || !xk4 || !x_err_comp_sum || !v_err_comp_sum)
    {
        printf("Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(store_npts * sizeof(double));
    double *sol_dt = malloc(store_npts * sizeof(double));

    if (!sol_state || !sol_time || !sol_dt)
    {
        printf("Error: Failed to allocate memory for solution output\n");
        goto err_sol_output_memory;
    }

    // Initial value
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sol_state[i * 3 + j] = x[i * 3 + j];
            sol_state[objects_count * 3 + i * 3 + j] = v[i * 3 + j];
        }
    }
    sol_time[0] = 0.0;
    for (int i = 0; i < store_npts; i++)
    {
        sol_dt[i] = dt;
    }

    // Main Loop
    while ((count + 1) <= npts)
    {   
        acceleration(objects_count, x, vk1, m, G);
        memcpy(xk1, v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                temp_x[j * 3 + k] = x[j * 3 + k] + 0.5 * xk1[j * 3 + k] * dt;
                temp_v[j * 3 + k] = v[j * 3 + k] + 0.5 * vk1[j * 3 + k] * dt;
            }
        }
        acceleration(objects_count, temp_x, vk2, m, G);
        memcpy(xk2, temp_v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                temp_x[j * 3 + k] = x[j * 3 + k] + 0.5 * xk2[j * 3 + k] * dt;
                temp_v[j * 3 + k] = v[j * 3 + k] + 0.5 * vk2[j * 3 + k] * dt;
            }
        }
        acceleration(objects_count, temp_x, vk3, m, G);
        memcpy(xk3, temp_v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                temp_x[j * 3 + k] = x[j * 3 + k] + xk3[j * 3 + k] * dt;
                temp_v[j * 3 + k] = v[j * 3 + k] + vk3[j * 3 + k] * dt;
            }
        }
        acceleration(objects_count, temp_x, vk4, m, G);
        memcpy(xk4, temp_v, objects_count * 3 * sizeof(real));

        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        memcpy(temp_x, x, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            for (int k = 0; k < 3; k++)
            {
                v_err_comp_sum[j * 3 + k] += (vk1[j * 3 + k] + 2 * vk2[j * 3 + k] + 2 * vk3[j * 3 + k] + vk4[j * 3 + k]) * dt / 6.0;
                x_err_comp_sum[j * 3 + k] += (xk1[j * 3 + k] + 2 * xk2[j * 3 + k] + 2 * xk3[j * 3 + k] + xk4[j * 3 + k]) * dt / 6.0;

                v[j * 3 + k] = temp_v[j * 3 + k] + v_err_comp_sum[j * 3 + k];
                x[j * 3 + k] = temp_x[j * 3 + k] + x_err_comp_sum[j * 3 + k];

                v_err_comp_sum[j * 3 + k] += temp_v[j * 3 + k] - v[j * 3 + k];
                x_err_comp_sum[j * 3 + k] += temp_x[j * 3 + k] - x[j * 3 + k];
            }
        }

        // Store solution
        if ((count + 1) == npts)
        {
            memcpy(&sol_state[(store_npts - 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
            memcpy(&sol_state[(store_npts - 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            sol_time[store_npts - 1] = dt * count;
        }
        else if (((count + 1) % store_every_n == 0))
        {
            memcpy(&sol_state[(*store_count + 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
            memcpy(&sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            sol_time[*store_count + 1] = dt * (count + 1);
            if (!((*store_count + 1) == (store_npts - 1)))
            {
                (*store_count)++;
            }
        }
        count++;

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }
    }

    // Exit after simulation is finished
    free(x);
    free(v);
    free(a);
    free(temp_x);
    free(temp_v);
    free(vk1);
    free(vk2);
    free(vk3);
    free(vk4);
    free(xk1);
    free(xk2);
    free(xk3);
    free(xk4);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    solution->sol_state = sol_state;
    solution->sol_time = sol_time;
    solution->sol_dt = sol_dt;
    solution->m = m;
    solution->G = G;
    solution->objects_count = objects_count;

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_sol_output_memory:
    free(sol_state);
    free(sol_time);
    free(sol_dt);
err_calc_memory:
    free(temp_x);
    free(temp_v);
    free(a);
    free(vk1);
    free(vk2);
    free(vk3);
    free(vk4);
    free(xk1);
    free(xk2);
    free(xk3);
    free(xk4);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
err_custom_sys_memory:
    free(x);
    free(v);
    free(m);
err_default_sys_memory:    
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}

/**
 * \brief LeapFrog integrator
 * 
 * \param system Name of the system
 * \param dt Time step of the system
 * \param npts Number of time steps to be integrated
 * \param store_npts Number of points to be stored
 * \param store_every_n Store every nth point
 * \param store_count Pointer to the store count
 * \param custom_sys_x Array of the position vectors of customized system 
 * \param custom_sys_v Array of the velocity vectors of customized system 
 * \param custom_sys_m Array of the masses of customized system 
 * \param custom_sys_G Gravitational constant of customized system
 * \param custom_sys_object_count Number of objects in customized system
 * \param solution Pointer to a Solution struct, in order to store the solution
 * \param is_exit Pointer to flag that indicates whether user sent 
 *                KeyboardInterrupt in the main thread
 * 
 * \retval 0 If exit successfully
 * \retval 1 If failed to allocate memory
 * \retval 2 If KeyboardInterrupt in the main thread
 */
WIN32DLL_API int leapfrog(
    const char *restrict system,
    double dt,
    int64 npts,
    int store_npts,
    int store_every_n,
    int *restrict store_count,
    const double *restrict custom_sys_x,
    const double *restrict custom_sys_v,
    const double *restrict custom_sys_m,
    double custom_sys_G,
    int custom_sys_objects_count,
    Solutions *restrict solution,
    int *restrict is_exit
)
{   
    // Initialize default system
    real *x = NULL;
    real *v = NULL;
    real *m = NULL;
    int objects_count;
    real G;
    int intitial_sys_flag = initialize_system(system, &x, &v, &m, &objects_count, &G);

    // Custom system
    if (intitial_sys_flag == 2)
    {
        printf("Error: Failed to allocate memory to initialize default system\n");
        goto err_default_sys_memory;
    }
    else if (intitial_sys_flag == 1)
    {
        objects_count = custom_sys_objects_count;
        x = malloc(objects_count * 3 * sizeof(real));
        v = malloc(objects_count * 3 * sizeof(real));
        m = malloc(objects_count * sizeof(real));

        if (!x || !v || !m)
        {
            printf("Error: Failed to allocate memory for customized system\n");
            goto err_custom_sys_memory;
        }

        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x[i * 3 + j] = custom_sys_x[i * 3 + j];
                v[i * 3 + j] = custom_sys_v[i * 3 + j];
            }
            m[i] = custom_sys_m[i];
        }
        G = custom_sys_G;
    }

    // Allocate memory for calculation
    real *temp_x = malloc(objects_count * 3 * sizeof(real));
    real *temp_v = malloc(objects_count * 3 * sizeof(real));
    real *a_0 = malloc(objects_count * 3 * sizeof(real));
    real *a_1 = malloc(objects_count * 3 * sizeof(real));

    // Allocate memory for compensated summation
    real *x_err_comp_sum = calloc(objects_count * 3, sizeof(real));
    real *v_err_comp_sum = calloc(objects_count * 3, sizeof(real));

    if (!temp_x || !temp_v || !a_0 || !a_1 || !x_err_comp_sum || !v_err_comp_sum)
    {
        printf("Error: Failed to allocate memory for calculation\n");
        goto err_calc_memory;
    }

    // Allocate memory for solution output
    int64 count = 0;
    double *sol_state = malloc(store_npts * objects_count * 6 * sizeof(double));
    double *sol_time = malloc(store_npts * sizeof(double));
    double *sol_dt = malloc(store_npts * sizeof(double));

    if (!sol_state || !sol_time || !sol_dt)
    {
        printf("Error: Failed to allocate memory for solution output\n");
        goto err_sol_output_memory;
    }

    // Initial value
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sol_state[i * 3 + j] = x[i * 3 + j];
            sol_state[objects_count * 3 + i * 3 + j] = v[i * 3 + j];
        }
    }
    sol_time[0] = 0.0;
    for (int i = 0; i < store_npts; i++)
    {
        sol_dt[i] = dt;
    }

    // Main Loop
    acceleration(objects_count, x, a_1, m, G);
    while ((count + 1) <= npts)
    {       
        // Use a_1 from last iteration as a_0
        memcpy(a_0, a_1, objects_count * 3 * sizeof(real));

        // Calculate x
        memcpy(temp_x, x, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                x_err_comp_sum[i * 3 + j] += v[i * 3 + j] * dt + 0.5 * a_0[i * 3 + j] * dt * dt;
                x[i * 3 + j] = temp_x[i * 3 + j] + x_err_comp_sum[i * 3 + j];
                x_err_comp_sum[i * 3 + j] += temp_x[i * 3 + j] - x[i * 3 + j];
            }
        }

        // Calculate v
        acceleration(objects_count, x, a_1, m, G);
        memcpy(temp_v, v, objects_count * 3 * sizeof(real));
        for (int i = 0; i < objects_count; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                v_err_comp_sum[i * 3 + j] += 0.5 * (a_0[i * 3 + j] + a_1[i * 3 + j]) * dt;
                v[i * 3 + j] = temp_v[i * 3 + j] + v_err_comp_sum[i * 3 + j];
                v_err_comp_sum[i * 3 + j] += temp_v[i * 3 + j] - v[i * 3 + j];
            }
        }

        // Store solution
        if ((count + 1) == npts)
        {
            memcpy(&sol_state[(store_npts - 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
            memcpy(&sol_state[(store_npts - 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            sol_time[store_npts - 1] = dt * count;
        }
        else if (((count + 1) % store_every_n == 0))
        {
            memcpy(&sol_state[(*store_count + 1) * objects_count * 6], x, objects_count * 6 * sizeof(double));
            memcpy(&sol_state[(*store_count + 1) * objects_count * 6 + objects_count * 3], v, objects_count * 6 * sizeof(double));
            sol_time[*store_count + 1] = dt * (count + 1);
            if (!((*store_count + 1) == (store_npts - 1)))
            {
                (*store_count)++;
            }
        }
        count++;

        // Check if user sends KeyboardInterrupt in main thread
        if (*is_exit)
        {
            goto err_user_exit;
        }
    }

    // Exit after simulation is finished
    free(x);
    free(v);
    free(temp_x);
    free(temp_v);
    free(a_0);
    free(a_1);
    free(x_err_comp_sum);
    free(v_err_comp_sum);

    solution->sol_state = sol_state;
    solution->sol_time = sol_time;
    solution->sol_dt = sol_dt;
    solution->m = m;
    solution->G = G;
    solution->objects_count = objects_count;

    return 0;

err_user_exit: // User sends KeyboardInterrupt in main thread
err_sol_output_memory:
    free(sol_state);
    free(sol_time);
    free(sol_dt);
err_calc_memory:
    free(temp_x);
    free(temp_v);
    free(a_0);
    free(a_1);
    free(x_err_comp_sum);
    free(v_err_comp_sum);
err_custom_sys_memory:
    free(x);
    free(v);
    free(m);
err_default_sys_memory: 
    if (*is_exit)
    {
        return 2;   // User sends KeyboardInterrupt in main thread
    }
    else
    {
        return 1;
    }
}