#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h> // For testing

#ifdef WIN32DLL_EXPORTS
    #define WIN32DLL_API __declspec(dllexport)
#else
    #define WIN32DLL_API 
#endif

#define real double

real abs_max_vec(const real *restrict vec, int vec_length);
real abs_max_vec_array(const real (*restrict arr)[3], int objects_count);
real vec_norm(const real *restrict vec, int vec_length);
void compute_energy(
    int objects_count, 
    int npts,
    int *restrict count, 
    real *restrict energy, 
    const real (*restrict sol_state)[objects_count * 6], 
    const real *restrict m, 
    real G
);
void acceleration(int objects_count, const real (*restrict x)[3], real (*restrict a)[3], const real *restrict m, real G);
void euler(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    real tf, 
    unsigned long npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    unsigned long *restrict count,
    unsigned int *restrict store_count
);
void euler_cromer(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    real tf, 
    unsigned long npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    unsigned long *restrict count,
    unsigned int *restrict store_count
);
void rk4(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3],  
    real dt, 
    real tf, 
    unsigned long npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    unsigned long *restrict count,
    unsigned int *restrict store_count
);
void leapfrog(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    real tf, 
    unsigned long npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    unsigned long *restrict count,
    unsigned int *restrict store_count
);
int rk_embedded(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real *restrict t, 
    real *restrict dt, 
    real tf, 
    int store_every_n,
    unsigned int *restrict store_count,
    unsigned int *restrict count, 
    const real *restrict m, 
    real G, 
    int power,
    int power_test,
    int len_coeff,
    const real (*restrict coeff)[len_coeff],
    int len_weights,
    const real *restrict weights,
    const real *restrict weights_test,
    real abs_tolerance,
    real rel_tolerance,
    real (*restrict sol_state)[6 * objects_count],
    int len_sol_time,
    real *restrict sol_time,
    int len_sol_dt,
    real *restrict sol_dt
);
int ias15(
    int objects_count, 
    int dim_nodes,
    const real *restrict nodes,
    const real *restrict aux_c, 
    const real *restrict aux_r,
    real *restrict aux_b0,
    real *restrict aux_b,
    real *restrict aux_g,
    real *restrict aux_e,
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real (*restrict a)[3],
    const real *restrict m, 
    real G,     
    real *restrict t, 
    real *restrict dt, 
    real tf, 
    int store_every_n,
    unsigned int *restrict store_count,
    unsigned int *restrict count, 
    real tolerance,
    real tolerance_pc,
    real (*restrict sol_state)[6 * objects_count],
    int len_sol_time,
    real *restrict sol_time,
    int len_sol_dt,
    real *restrict sol_dt,
    real safety_fac,
    real exponent,
    int *restrict ias15_refine_flag
);
void ias15_step(
    int objects_count,
    int dim_nodes,
    real (*restrict x)[3],
    real (*restrict v)[3],
    real (*restrict a)[3],
    const real *restrict m,
    real G,
    real *restrict t,
    real *restrict dt,
    real tf,
    const real *restrict nodes,
    real *restrict aux_b0,
    real *restrict aux_b,
    const real *restrict aux_c,
    real *restrict aux_e,
    real *restrict aux_g,
    const real *restrict aux_r,
    real tolerance,
    real tolerance_pc,
    real exponent,
    real safety_fac,
    int *restrict ias15_refine_flag
);
void ias15_approx_pos(
    int objects_count,
    real (*restrict x)[3],
    const real (*restrict v)[3],
    const real (*restrict a)[3],
    real node,
    real *restrict aux_b,
    real dt
);
void ias15_approx_vel(
    int objects_count,
    real (*restrict v)[3],
    const real (*restrict a)[3],
    real node,
    real *restrict aux_b,
    real dt
);
void ias15_compute_aux_b(
    int objects_count,
    int dim_nodes,
    real *restrict aux_b,
    const real *restrict aux_g,
    const real *restrict aux_c,
    int i
);
void ias15_compute_aux_g(
    int objects_count,
    int dim_nodes,
    real *restrict aux_g,
    const real *restrict aux_r,
    const real *restrict aux_a,
    int i
);
void ias15_refine_aux_b(
    int objects_count,
    int dim_nodes,
    real *restrict aux_b,
    real *restrict aux_e,
    real dt,
    real dt_new,
    int ias15_refine_flag
);

// Find the max absolute value in a 1D array
WIN32DLL_API real abs_max_vec(const real *restrict vec, int vec_length)
{
    real max = fabs(vec[0]);
    for (int i = 1; i < vec_length; i++)
    {
        max = fmax(max, fabs(vec[i]));
    }

    return max;
}

// Find the max absolute value in a array with 3D vectors
WIN32DLL_API real abs_max_vec_array(const real (*restrict arr)[3], int objects_count)
{
    real max = 0;
    for (int i = 0; i < objects_count; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            max = fmax(max, fabs(arr[i][j]));
        }
    }

    return max;
}

// Find the norm of 1D vector
WIN32DLL_API real vec_norm(const real *restrict vec, int vec_length)
{   
    real sum = 0.0;
    if (vec_length == 3) 
    {
        sum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    }
    else
    {
        for (int i = 0; i < vec_length; i++) 
        {
            sum += vec[i] * vec[i];
        }
    }
    return sqrt(sum);
}

WIN32DLL_API void compute_energy(
    int objects_count, 
    int npts,
    int *restrict count, 
    real *restrict energy, 
    const real (*restrict sol_state)[objects_count * 6], 
    const real *restrict m, 
    real G
)
{
    int progress_percentage = (int) round((float) *count / (float) npts * 100.0);

    real temp_vec[3];

    while (1)
    {   
        for (int i = 0; i < objects_count; i++)
        {
            // KE
            energy[*count] += (
                0.5 * m[i] 
                * pow(vec_norm(&sol_state[*count][(objects_count + i) * 3], 3), 2)
            );

            // PE
            for (int j = i + 1; j < objects_count; j++)
            {
                temp_vec[0] = (
                    sol_state[*count][i * 3 + 0] 
                    - sol_state[*count][j * 3 + 0]
                );
                temp_vec[1] = (
                    sol_state[*count][i * 3 + 1] 
                    - sol_state[*count][j * 3 + 1]
                );
                temp_vec[2] = (
                    sol_state[*count][i * 3 + 2] 
                    - sol_state[*count][j * 3 + 2]
                );
                energy[*count] -= (
                    G * m[i] * m[j]
                    / vec_norm(temp_vec, 3)
                );
            }
        }
        *count += 1;

        if (*count >= npts)
        {
            break;
        }

        // Exit to update progress bar
        if ((int) round((float) *count / (float) npts * 100.0) > progress_percentage)
        {   
            break;
        }
    }
}

WIN32DLL_API void acceleration(int objects_count, const real (*restrict x)[3], real (*restrict a)[3], const real *restrict m, real G)
{   
    real R_norm;
    real temp_value;
    real temp_vec[3];
    real R[3];

    // Empty the input array
    memset(a, 0, objects_count * 3 * sizeof(real));

    for(int i = 0; i < objects_count; i++)
    {
        for(int j = i + 1; j < objects_count; j++)
        {
            // Calculate \vec{R} and its norm
            R[0] = x[i][0] - x[j][0];
            R[1] = x[i][1] - x[j][1];
            R[2] = x[i][2] - x[j][2];
            R_norm = vec_norm(R, 3);

            // Calculate the acceleration
            temp_value = G / (R_norm * R_norm * R_norm);
            temp_vec[0] = temp_value * R[0];
            temp_vec[1] = temp_value * R[1];
            temp_vec[2] = temp_value * R[2];
            a[i][0] += - temp_vec[0] * m[j];
            a[i][1] += - temp_vec[1] * m[j];
            a[i][2] += - temp_vec[2] * m[j];
            a[j][0] += temp_vec[0] * m[i];
            a[j][1] += temp_vec[1] * m[i];
            a[j][2] += temp_vec[2] * m[i];
        }
    }
}

WIN32DLL_API void euler(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    real tf, 
    unsigned long npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    unsigned long *restrict count,
    unsigned int *restrict store_count
)
{   
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*count * 100 / npts);

    // Main Loop
    while ((*count + 1) <= npts)
    {   
        acceleration(objects_count, x, a, m, G);
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            x[j][0] += v[j][0] * dt;
            x[j][1] += v[j][1] * dt;
            x[j][2] += v[j][2] * dt;
            v[j][0] += a[j][0] * dt;
            v[j][1] += a[j][1] * dt;
            v[j][2] += a[j][2] * dt;

            // Store solution
            if ((*count + 1) % store_every_n == 0)
            {
                sol_state[*store_count + 1][j * 3] = x[j][0];
                sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }

            if ((*count + 2) == npts)
            {
                sol_state[store_npts - 1][j * 3] = x[j][0];
                sol_state[store_npts - 1][j * 3 + 1] = x[j][1];
                sol_state[store_npts - 1][j * 3 + 2] = x[j][2];
                sol_state[store_npts - 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }
        }    
        if ((*count + 1) % store_every_n == 0)
        {
            *store_count += 1;
        }

        *count += 1;

        // Exit to update progress bar
        if ((int) (*count * 100 / npts) > progress_percentage)
        {   
            free(a);
            return;
        }
    }
}

WIN32DLL_API void euler_cromer(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    real tf, 
    unsigned long npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    unsigned long *restrict count,
    unsigned int *restrict store_count
)
{   
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*count * 100 / npts);

    // Main Loop
    while ((*count + 1) <= npts)
    {   
        acceleration(objects_count, x, a, m, G);
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            v[j][0] += a[j][0] * dt;
            v[j][1] += a[j][1] * dt;
            v[j][2] += a[j][2] * dt;
            x[j][0] += v[j][0] * dt;
            x[j][1] += v[j][1] * dt;
            x[j][2] += v[j][2] * dt;

            // Store solution
            if ((*count + 1) % store_every_n == 0)
            {
                sol_state[*store_count + 1][j * 3] = x[j][0];
                sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }

            if ((*count + 2) == npts)
            {
                sol_state[store_npts - 1][j * 3] = x[j][0];
                sol_state[store_npts - 1][j * 3 + 1] = x[j][1];
                sol_state[store_npts - 1][j * 3 + 2] = x[j][2];
                sol_state[store_npts - 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }
        }    
        if ((*count + 1) % store_every_n == 0)
        {
            *store_count += 1;
        }

        *count += 1;

        // Exit to update progress bar
        if ((int) (*count * 100 / npts) > progress_percentage)
        {   
            free(a);
            return;
        }
    }
}

WIN32DLL_API void rk4(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    real tf, 
    unsigned long npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    unsigned long *restrict count,
    unsigned int *restrict store_count
)
{
    real (*temp_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v)[3] = malloc(objects_count * 3 * sizeof(real));

    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    real (*vk1)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*vk2)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*vk3)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*vk4)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*xk1)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*xk2)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*xk3)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*xk4)[3] = malloc(objects_count * 3 * sizeof(real));

    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*count * 100 / npts);

    // Main Loop
    while ((*count + 1) <= npts)
    {   
        acceleration(objects_count, x, a, m, G);
        memcpy(vk1, a, objects_count * 3 * sizeof(real));
        memcpy(xk1, v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            temp_x[j][0] = x[j][0] + 0.5 * xk1[j][0] * dt;
            temp_x[j][1] = x[j][1] + 0.5 * xk1[j][1] * dt;
            temp_x[j][2] = x[j][2] + 0.5 * xk1[j][2] * dt;
            temp_v[j][0] = v[j][0] + 0.5 * vk1[j][0] * dt;
            temp_v[j][1] = v[j][1] + 0.5 * vk1[j][1] * dt;
            temp_v[j][2] = v[j][2] + 0.5 * vk1[j][2] * dt;
        }
        acceleration(objects_count, temp_x, a, m, G);
        memcpy(vk2, a, objects_count * 3 * sizeof(real));
        memcpy(xk2, temp_v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            temp_x[j][0] = x[j][0] + 0.5 * xk2[j][0] * dt;
            temp_x[j][1] = x[j][1] + 0.5 * xk2[j][1] * dt;
            temp_x[j][2] = x[j][2] + 0.5 * xk2[j][2] * dt;
            temp_v[j][0] = v[j][0] + 0.5 * vk2[j][0] * dt;
            temp_v[j][1] = v[j][1] + 0.5 * vk2[j][1] * dt;
            temp_v[j][2] = v[j][2] + 0.5 * vk2[j][2] * dt;
        }
        acceleration(objects_count, temp_x, a, m, G);
        memcpy(vk3, a, objects_count * 3 * sizeof(real));
        memcpy(xk3, temp_v, objects_count * 3 * sizeof(real));

        for (int j = 0; j < objects_count; j++)
        {
            temp_x[j][0] = x[j][0] + xk3[j][0] * dt;
            temp_x[j][1] = x[j][1] + xk3[j][1] * dt;
            temp_x[j][2] = x[j][2] + xk3[j][2] * dt;
            temp_v[j][0] = v[j][0] + vk3[j][0] * dt;
            temp_v[j][1] = v[j][1] + vk3[j][1] * dt;
            temp_v[j][2] = v[j][2] + vk3[j][2] * dt;
        }
        acceleration(objects_count, temp_x, a, m, G);
        memcpy(vk4, a, objects_count * 3 * sizeof(real));
        memcpy(xk4, temp_v, objects_count * 3 * sizeof(real));


        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            v[j][0] += (vk1[j][0] + 2 * vk2[j][0] + 2 * vk3[j][0] + vk4[j][0]) * dt / 6.0;
            v[j][1] += (vk1[j][1] + 2 * vk2[j][1] + 2 * vk3[j][1] + vk4[j][1]) * dt / 6.0;
            v[j][2] += (vk1[j][2] + 2 * vk2[j][2] + 2 * vk3[j][2] + vk4[j][2]) * dt / 6.0;
            x[j][0] += (xk1[j][0] + 2 * xk2[j][0] + 2 * xk3[j][0] + xk4[j][0]) * dt / 6.0;
            x[j][1] += (xk1[j][1] + 2 * xk2[j][1] + 2 * xk3[j][1] + xk4[j][1]) * dt / 6.0;
            x[j][2] += (xk1[j][2] + 2 * xk2[j][2] + 2 * xk3[j][2] + xk4[j][2]) * dt / 6.0;

            // Store solution
            if ((*count + 1) % store_every_n == 0)
            {
                sol_state[*store_count + 1][j * 3] = x[j][0];
                sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }

            if ((*count + 2) == npts)
            {
                sol_state[store_npts - 1][j * 3] = x[j][0];
                sol_state[store_npts - 1][j * 3 + 1] = x[j][1];
                sol_state[store_npts - 1][j * 3 + 2] = x[j][2];
                sol_state[store_npts - 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }
        }    
        if ((*count + 1) % store_every_n == 0)
        {
            *store_count += 1;
        }

        *count += 1;

        // Exit to update progress bar
        if ((int) (*count * 100 / npts) > progress_percentage)
        {   
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
            return;
        }
    }
}

WIN32DLL_API void leapfrog(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real dt, 
    real tf, 
    unsigned long npts, 
    const real *restrict m, 
    real G, 
    real (*restrict sol_state)[6 * objects_count],
    int store_every_n,
    int store_npts,
    unsigned long *restrict count,
    unsigned int *restrict store_count
)
{   
    real (*a_0)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a_1)[3] = malloc(objects_count * 3 * sizeof(real));

    int is_initialize = 1;

    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*count * 100 / npts);

    // Main Loop
    while ((*count + 1) <= npts)
    {       
        if (is_initialize == 1)
        {
            acceleration(objects_count, x, a_0, m, G);
            is_initialize = 0;
        }
        else 
        {
            // Use a_1 from last iteration as a_0
            memcpy(a_0, a_1, objects_count * 3 * sizeof(real));
        }
           
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            x[j][0] += v[j][0] * dt + 0.5 * a_0[j][0] * dt * dt;
            x[j][1] += v[j][1] * dt + 0.5 * a_0[j][1] * dt * dt;
            x[j][2] += v[j][2] * dt + 0.5 * a_0[j][2] * dt * dt;
        }    
        acceleration(objects_count, x, a_1, m, G);
        for (int j = 0; j < objects_count; j++)
        {
            // Calculation
            v[j][0] += 0.5 * (a_0[j][0] + a_1[j][0]) * dt;
            v[j][1] += 0.5 * (a_0[j][1] + a_1[j][1]) * dt;
            v[j][2] += 0.5 * (a_0[j][2] + a_1[j][2]) * dt;

            // Store solution
            if ((*count + 1) % store_every_n == 0)
            {
                sol_state[*store_count + 1][j * 3] = x[j][0];
                sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }

            if ((*count + 2) == npts)
            {
                sol_state[store_npts - 1][j * 3] = x[j][0];
                sol_state[store_npts - 1][j * 3 + 1] = x[j][1];
                sol_state[store_npts - 1][j * 3 + 2] = x[j][2];
                sol_state[store_npts - 1][objects_count * 3 + j * 3] = v[j][0];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[store_npts - 1][objects_count * 3 + j * 3 + 2] = v[j][2];
            }
        }    
        if ((*count + 1) % store_every_n == 0)
        {
            *store_count += 1;
        }

        *count += 1;

        // Exit to update progress bar
        if ((int) (*count * 100 / npts) > progress_percentage)
        {   
            free(a_0);
            free(a_1);
            return;
        }
    }
}

WIN32DLL_API int rk_embedded(
    int objects_count, 
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real *restrict t, 
    real *restrict dt, 
    real tf, 
    int store_every_n,
    unsigned int *restrict store_count,
    unsigned int *restrict count, 
    const real *restrict m, 
    real G, 
    int power,
    int power_test,
    int len_coeff,
    const real (*restrict coeff)[len_coeff],
    int len_weights,
    const real *restrict weights,
    const real *restrict weights_test,
    real abs_tolerance,
    real rel_tolerance,
    real (*restrict sol_state)[6 * objects_count],
    int len_sol_time,
    real *restrict sol_time,
    int len_sol_dt,
    real *restrict sol_dt
)
{
    // Initialization
    int stages = len_weights;
    int min_power = fmin(power, power_test);

    real *error_estimation_delta_weights = malloc(len_weights * sizeof(real));
    for (int stage = 0; stage < stages; stage++)
    {
        error_estimation_delta_weights[stage] = weights[stage] - weights_test[stage];
    }

    // Safety factors for step-size control:
    real safety_fac_max = 6.0;
    real safety_fac_min = 0.33;
    real safety_fac = pow(0.38, (1.0 / (1.0 + (real) min_power)));

    // Initialize arrays and values
    real sum, error, dt_new; 
    real (*v_1)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*x_1)[3] = malloc(objects_count * 3 * sizeof(real));
    real *vk = malloc(stages * objects_count * 3 * sizeof(real));
    real *xk = malloc(stages * objects_count * 3 * sizeof(real));    
    real (*temp_a)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*temp_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*error_estimation_delta_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*error_estimation_delta_x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*tolerance_scale_v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*tolerance_scale_x)[3] = malloc(objects_count * 3 * sizeof(real));

    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*t / tf * 100);

    // Main Loop
    while (1)
    {
        // Calculate xk and vk
        acceleration(objects_count, x, temp_a, m, G);
        memcpy(vk, temp_a, objects_count * 3 * sizeof(real));
        memcpy(xk, v, objects_count * 3 * sizeof(real));

        for (int stage = 1; stage < stages; stage++)
        {
            memset(temp_v, 0, objects_count * 3 * sizeof(real));
            memset(temp_x, 0, objects_count * 3 * sizeof(real));         

            for (int j = 0; j < stage; j++)
            {
                for (int k = 0; k < objects_count; k++)
                {
                    temp_v[k][0] += coeff[stage - 1][j] * vk[j * objects_count * 3 + k * 3];
                    temp_v[k][1] += coeff[stage - 1][j] * vk[j * objects_count * 3 + k * 3 + 1];
                    temp_v[k][2] += coeff[stage - 1][j] * vk[j * objects_count * 3 + k * 3 + 2];
                    temp_x[k][0] += coeff[stage - 1][j] * xk[j * objects_count * 3 + k * 3];
                    temp_x[k][1] += coeff[stage - 1][j] * xk[j * objects_count * 3 + k * 3 + 1];
                    temp_x[k][2] += coeff[stage - 1][j] * xk[j * objects_count * 3 + k * 3 + 2];
                }
            }

            for (int k = 0; k < objects_count; k++)
            {
                temp_x[k][0] = x[k][0] + *dt * temp_x[k][0];
                temp_x[k][1] = x[k][1] + *dt * temp_x[k][1];
                temp_x[k][2] = x[k][2] + *dt * temp_x[k][2];
                temp_v[k][0] = v[k][0] + *dt * temp_v[k][0];
                temp_v[k][1] = v[k][1] + *dt * temp_v[k][1];
                temp_v[k][2] = v[k][2] + *dt * temp_v[k][2];
            }
            acceleration(objects_count, temp_x, temp_a, m, G);
            memcpy(&vk[stage * objects_count * 3], temp_a, objects_count * 3 * sizeof(real));
            memcpy(&xk[stage * objects_count * 3], temp_v, objects_count * 3 * sizeof(real));
        }

        // Calculate x_1, v_1 and also delta x, delta v for error estimation
        memset(temp_v, 0, objects_count * 3 * sizeof(real));
        memset(temp_x, 0, objects_count * 3 * sizeof(real));       
        memset(error_estimation_delta_v, 0, objects_count * 3 * sizeof(real));
        memset(error_estimation_delta_x, 0, objects_count * 3 * sizeof(real));       
        for(int stage = 0; stage < stages; stage++)
        {
            for (int k = 0; k < objects_count; k++)
            {
                temp_v[k][0] += weights[stage] * vk[stage * objects_count * 3 + k * 3];
                temp_v[k][1] += weights[stage] * vk[stage * objects_count * 3 + k * 3 + 1];
                temp_v[k][2] += weights[stage] * vk[stage * objects_count * 3 + k * 3 + 2];
                temp_x[k][0] += weights[stage] * xk[stage * objects_count * 3 + k * 3];
                temp_x[k][1] += weights[stage] * xk[stage * objects_count * 3 + k * 3 + 1];
                temp_x[k][2] += weights[stage] * xk[stage * objects_count * 3 + k * 3 + 2];

                error_estimation_delta_v[k][0] += *dt * error_estimation_delta_weights[stage] * vk[stage * objects_count * 3 + k * 3];
                error_estimation_delta_v[k][1] += *dt * error_estimation_delta_weights[stage] * vk[stage * objects_count * 3 + k * 3 + 1];
                error_estimation_delta_v[k][2] += *dt * error_estimation_delta_weights[stage] * vk[stage * objects_count * 3 + k * 3 + 2];
                error_estimation_delta_x[k][0] += *dt * error_estimation_delta_weights[stage] * xk[stage * objects_count * 3 + k * 3];
                error_estimation_delta_x[k][1] += *dt * error_estimation_delta_weights[stage] * xk[stage * objects_count * 3 + k * 3 + 1];
                error_estimation_delta_x[k][2] += *dt * error_estimation_delta_weights[stage] * xk[stage * objects_count * 3 + k * 3 + 2];
            }
        }

        for (int k = 0; k < objects_count; k++)
        {
            v_1[k][0] = v[k][0] + *dt * temp_v[k][0];
            v_1[k][1] = v[k][1] + *dt * temp_v[k][1];
            v_1[k][2] = v[k][2] + *dt * temp_v[k][2];
            x_1[k][0] = x[k][0] + *dt * temp_x[k][0];
            x_1[k][1] = x[k][1] + *dt * temp_x[k][1];
            x_1[k][2] = x[k][2] + *dt * temp_x[k][2];
        }

        // Error calculation
        for (int k = 0; k < objects_count; k++)
        {
            tolerance_scale_v[k][0] = abs_tolerance + fmax(fabs(v[k][0]), fabs(v_1[k][0])) * rel_tolerance;
            tolerance_scale_v[k][1] = abs_tolerance + fmax(fabs(v[k][1]), fabs(v_1[k][1])) * rel_tolerance;
            tolerance_scale_v[k][2] = abs_tolerance + fmax(fabs(v[k][2]), fabs(v_1[k][2])) * rel_tolerance;
            tolerance_scale_x[k][0] = abs_tolerance + fmax(fabs(x[k][0]), fabs(x_1[k][0])) * rel_tolerance;
            tolerance_scale_x[k][1] = abs_tolerance + fmax(fabs(x[k][1]), fabs(x_1[k][1])) * rel_tolerance;
            tolerance_scale_x[k][2] = abs_tolerance + fmax(fabs(x[k][2]), fabs(x_1[k][2])) * rel_tolerance;
        }

        // Sum up all the elements of x/tol and v/tol, 
        // square and divide by the total number of elements
        sum = 0.0;
        for (int k = 0; k < objects_count; k++)
        {
            sum += (error_estimation_delta_v[k][0] / tolerance_scale_v[k][0]) * (error_estimation_delta_v[k][0] / tolerance_scale_v[k][0]);
            sum += (error_estimation_delta_v[k][1] / tolerance_scale_v[k][1]) * (error_estimation_delta_v[k][1] / tolerance_scale_v[k][1]);
            sum += (error_estimation_delta_v[k][2] / tolerance_scale_v[k][2]) * (error_estimation_delta_v[k][2] / tolerance_scale_v[k][2]); 
            sum += (error_estimation_delta_x[k][0] / tolerance_scale_x[k][0]) * (error_estimation_delta_x[k][0] / tolerance_scale_x[k][0]); 
            sum += (error_estimation_delta_x[k][1] / tolerance_scale_x[k][1]) * (error_estimation_delta_x[k][1] / tolerance_scale_x[k][1]); 
            sum += (error_estimation_delta_x[k][2] / tolerance_scale_x[k][2]) * (error_estimation_delta_x[k][2] / tolerance_scale_x[k][2]); 
        }
        error = sqrt(sum / (objects_count * 3 * 2));

        if (error <= 1 || *dt == tf * 1e-12)
        {
            // Advance step
            *t += *dt; 
            memcpy(x, x_1, objects_count * 3 * sizeof(real));
            memcpy(v, v_1, objects_count * 3 * sizeof(real));
            *count += 1;

            // Store step
            if ((*count + 1) % store_every_n == 0)
            {
                sol_time[*store_count + 1] = *t;
                sol_dt[*store_count + 1] = *dt;
                for (int j = 0; j < objects_count; j++)
                {
                    sol_state[*store_count + 1][j * 3] = x[j][0];
                    sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                    sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
                }  
                *store_count += 1;
            }

            if (*t == tf)
            {
                sol_time[*store_count] = *t;
                sol_dt[*store_count] = *dt;
                for (int j = 0; j < objects_count; j++)
                {
                    sol_state[*store_count][j * 3] = x[j][0];
                    sol_state[*store_count][j * 3 + 1] = x[j][1];
                    sol_state[*store_count][j * 3 + 2] = x[j][2];
                    sol_state[*store_count][objects_count * 3 + j * 3] = v[j][0];
                    sol_state[*store_count][objects_count * 3 + j * 3 + 1] = v[j][1];
                    sol_state[*store_count][objects_count * 3 + j * 3 + 2] = v[j][2];
                }  
            }

            // Check buffer size and quit if full
            if ((*store_count + 1) == len_sol_time)
            {
                free(error_estimation_delta_weights);
                free(v_1);
                free(x_1);
                free(vk);
                free(xk);
                free(temp_a);
                free(temp_v);
                free(temp_x);
                free(error_estimation_delta_v);
                free(error_estimation_delta_x);
                free(tolerance_scale_v);
                free(tolerance_scale_x);
                return 1;
            } 
        }

        // Calculate dt
        if (error != 0.0)   // Prevent division by zero
        {
            dt_new = *dt * safety_fac / pow(error, (1.0 / (1.0 + (real) min_power)));
        }
        else
        {
            dt_new = *dt;
        }
        
        if (dt_new > safety_fac_max * *dt) 
        {
            *dt *= safety_fac_max;
        }
        else if (dt_new < safety_fac_min * *dt)
        {
            *dt *= safety_fac_min;
        }
        else
        {
            *dt = dt_new;
        }

        if (dt_new / tf < 1e-12)
        {
            *dt = tf * 1e-12;
        }

        // Correct overshooting
        if (*t + *dt > tf)
        {
            *dt = tf - *t;
        }
        
        // Exit if finished
        if (*t >= tf)
        {
            free(error_estimation_delta_weights);
            free(v_1);
            free(x_1);
            free(vk);
            free(xk);
            free(temp_a);
            free(temp_v);
            free(temp_x);
            free(error_estimation_delta_v);
            free(error_estimation_delta_x);
            free(tolerance_scale_v);
            free(tolerance_scale_x);
            return 2;
        }

        // Exit to update progress bar
        if ((int) (*t / tf * 100.0) > progress_percentage)
        {
            free(error_estimation_delta_weights);
            free(v_1);
            free(x_1);
            free(vk);
            free(xk);
            free(temp_a);
            free(temp_v);
            free(temp_x);
            free(error_estimation_delta_v);
            free(error_estimation_delta_x);
            free(tolerance_scale_v);
            free(tolerance_scale_x);
            return 0;
        }
    }
}

WIN32DLL_API int ias15(
    int objects_count, 
    int dim_nodes,
    const real *restrict nodes,
    const real *restrict aux_c, 
    const real *restrict aux_r,
    real *restrict aux_b0,
    real *restrict aux_b,
    real *restrict aux_g,
    real *restrict aux_e,
    real (*restrict x)[3], 
    real (*restrict v)[3], 
    real (*restrict a)[3], 
    const real *restrict m, 
    real G,     
    real *restrict t, 
    real *restrict dt, 
    real tf, 
    int store_every_n,
    unsigned int *store_count,
    unsigned int *count, 
    real tolerance,
    real tolerance_pc,
    real (*restrict sol_state)[6 * objects_count],
    int len_sol_time,
    real *restrict sol_time,
    int len_sol_dt,
    real *restrict sol_dt,
    real safety_fac,
    real exponent,
    int *restrict ias15_refine_flag
)
{
    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*t / tf * 100);

    while (1)
    {
        ias15_step(
            objects_count,
            dim_nodes,
            x,
            v,
            a,
            m,
            G,
            t,
            dt,
            tf,
            nodes,
            aux_b0,
            aux_b,
            aux_c,
            aux_e,
            aux_g,
            aux_r,
            tolerance,
            tolerance_pc,
            exponent,
            safety_fac,
            ias15_refine_flag
        );

        // Update count
        *count += 1;

            // Store step
            if ((*count + 1) % store_every_n == 0)
            {
                sol_time[*store_count + 1] = *t;
                sol_dt[*store_count + 1] = *dt;
                for (int j = 0; j < objects_count; j++)
                {
                    sol_state[*store_count + 1][j * 3] = x[j][0];
                    sol_state[*store_count + 1][j * 3 + 1] = x[j][1];
                    sol_state[*store_count + 1][j * 3 + 2] = x[j][2];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3] = v[j][0];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
                    sol_state[*store_count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
                }  
                *store_count += 1;
            }

            if (*t == tf)
            {
                sol_time[*store_count] = *t;
                sol_dt[*store_count] = *dt;
                for (int j = 0; j < objects_count; j++)
                {
                    sol_state[*store_count][j * 3] = x[j][0];
                    sol_state[*store_count][j * 3 + 1] = x[j][1];
                    sol_state[*store_count][j * 3 + 2] = x[j][2];
                    sol_state[*store_count][objects_count * 3 + j * 3] = v[j][0];
                    sol_state[*store_count][objects_count * 3 + j * 3 + 1] = v[j][1];
                    sol_state[*store_count][objects_count * 3 + j * 3 + 2] = v[j][2];
                }  
            }

        // Check buffer size and quit if full
        if ((*store_count + 1) == len_sol_time)
        {
            return 1;
        } 

        // End simulation as t = tf
        if (*t >= tf)
        {
            return 2;
        }

        // Exit to update progress bar
        if ((int) (*t / tf * 100.0) > progress_percentage)
        {
            return 0;
        }
    }
}

// Advance IAS15 for one step
WIN32DLL_API void ias15_step(
    int objects_count,
    int dim_nodes,
    real (*restrict x0)[3],
    real (*restrict v0)[3],
    real (*restrict a0)[3],
    const real *restrict m,
    real G,
    real *restrict t,
    real *restrict dt,
    real tf,
    const real *restrict nodes,
    real *restrict aux_b0,
    real *restrict aux_b,
    const real *restrict aux_c,
    real *restrict aux_e,
    real *restrict aux_g,
    const real *restrict aux_r,
    real tolerance,
    real tolerance_pc,
    real exponent,
    real safety_fac,
    int *restrict ias15_refine_flag
)
{
    real *aux_a = malloc(dim_nodes * objects_count * 3 * sizeof(real));
    real (*temp_a)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*x)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*v)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));
    real *delta_b7 = malloc(objects_count * 3 * sizeof(real));
    real error, error_b7, dt_new;
    // Main Loop
    int ias15_integrate_flag = 0; 
    while (1)
    {   
        // Loop for predictor-corrector algorithm
        // 12 = max iterations
        memcpy(a, a0, objects_count * 3 * sizeof(real));
        for (int temp = 0; temp < 12; temp++)
        {
            for (int i = 0; i < dim_nodes; i++)
            {
                // Estimate position and velocity with current aux_b and nodes
                memcpy(x, x0, objects_count * 3 * sizeof(real));
                memcpy(v, v0, objects_count * 3 * sizeof(real));
                ias15_approx_pos(objects_count, x, v, a, nodes[i], aux_b, *dt);
                ias15_approx_vel(objects_count, v, a, nodes[i], aux_b, *dt);

                // Evaluate force function and store result
                acceleration(objects_count, x, temp_a, m, G);
                memcpy(&aux_a[i * objects_count * 3], temp_a, objects_count * 3 * sizeof(real));
                
                ias15_compute_aux_g(objects_count, dim_nodes, aux_g, aux_r, aux_a, i);
                ias15_compute_aux_b(objects_count, dim_nodes, aux_b, aux_g, aux_c, i);
            }

            // Estimate convergence
            for (int i = 0; i < objects_count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    delta_b7[i * 3 + j] = aux_b[(dim_nodes - 2) * objects_count * 3 + i * 3 + j] - aux_b0[(dim_nodes - 2) * objects_count * 3 + i * 3 + j];
                }
            }
            memcpy(aux_b0, aux_b, (dim_nodes - 1) * objects_count * 3 * sizeof(real));
            if ((abs_max_vec(delta_b7, objects_count * 3) / abs_max_vec(&aux_a[(dim_nodes - 1) * objects_count * 3], objects_count * 3)) < tolerance_pc)
            {
                break;
            }
        }

        // Advance step
        memcpy(x, x0, objects_count * 3 * sizeof(real));
        memcpy(v, v0, objects_count * 3 * sizeof(real));
        memcpy(a, a0, objects_count * 3 * sizeof(real));
        ias15_approx_pos(objects_count, x, v, a, 1.0, aux_b, *dt);
        ias15_approx_vel(objects_count, v, a, 1.0, aux_b, *dt);
        acceleration(objects_count, x, a, m, G);

        // Estimate relative error
        error_b7 = abs_max_vec(&aux_b[(dim_nodes - 2) * objects_count * 3], objects_count * 3) / abs_max_vec_array(a, objects_count);
        error = pow((error_b7 / tolerance), exponent);
        
        // Step-size for the next step
        if (error != 0)
        {
            dt_new = *dt / error;
        }
        else
        {
            dt_new = *dt;
        }

        // Accept the step
        if (error <= 1 || *dt == tf * 1e-12)
        {
            // Report accepted step
            ias15_integrate_flag = 1;
            *t += *dt;

            ias15_refine_aux_b(objects_count, dim_nodes, aux_b, aux_e, *dt, dt_new, *ias15_refine_flag);
            *ias15_refine_flag = 1;

            if (*t >= tf)
            {
                memcpy(x0, x, objects_count * 3 * sizeof(real));
                memcpy(v0, v, objects_count * 3 * sizeof(real));
                memcpy(a0, a, objects_count * 3 * sizeof(real));
                free(aux_a);
                free(temp_a);
                free(x);
                free(v);
                free(a);
                free(delta_b7);      
                break;  
            }
        }

        // Step size for the next iteration
        if (dt_new > (*dt / safety_fac))
        {
            *dt = *dt / safety_fac;
        }
        else if (dt_new < *dt * safety_fac)
        {
            *dt = *dt * safety_fac;
        }
        else
        {
            *dt = dt_new;
        }

        if (dt_new / tf < 1e-12)
        {
            *dt = tf * 1e-12;
        }

        // Correct overshooting
        if (*t + *dt > tf)
        {
            *dt = tf - *t;
        }

        // Exit 
        if (ias15_integrate_flag > 0)
        {   
            memcpy(x0, x, objects_count * 3 * sizeof(real));
            memcpy(v0, v, objects_count * 3 * sizeof(real));
            memcpy(a0, a, objects_count * 3 * sizeof(real));
            free(aux_a);
            free(temp_a);
            free(x);
            free(v);
            free(a);
            free(delta_b7);       
            break;    
        }
    }
}

WIN32DLL_API void ias15_approx_pos(
    int objects_count,
    real (*restrict x)[3],
    const real (*restrict v)[3],
    const real (*restrict a)[3],
    real node,
    real *restrict aux_b,
    real dt
)
{   
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            x[j][k] = x[j][k] + dt * node * (
                v[j][k]
                + dt
                * node
                * (
                    a[j][k]
                    + node 
                    * (
                        aux_b[0 * objects_count * 3 + j * 3 + k] / 3.0
                        + node
                        * (
                            aux_b[1 * objects_count * 3 + j * 3 + k] / 6.0
                            + node
                            * (
                                aux_b[2 * objects_count * 3 + j * 3 + k] / 10.0
                                + node
                                * (
                                    aux_b[3 * objects_count * 3 + j * 3 + k] / 15.0
                                    + node
                                    * (
                                        aux_b[4 * objects_count * 3 + j * 3 + k] / 21.0
                                        + node * (
                                            aux_b[5 * objects_count * 3 + j * 3 + k] / 28.0 
                                            + node * aux_b[6 * objects_count * 3 + j * 3 + k] / 36.0
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

WIN32DLL_API void ias15_approx_vel(
    int objects_count,
    real (*restrict v)[3],
    const real (*restrict a)[3],
    real node,
    real *restrict aux_b,
    real dt
)
{
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            v[j][k] = v[j][k] + dt * node * (
                a[j][k]
                + node
                * (
                    aux_b[0 * objects_count * 3 + j * 3 + k] / 2.0
                    + node
                    * (
                        aux_b[1 * objects_count * 3 + j * 3 + k] / 3.0
                        + node
                        * (
                            aux_b[2 * objects_count * 3 + j * 3 + k] / 4.0
                            + node
                            * (
                                aux_b[3 * objects_count * 3 + j * 3 + k] / 5.0
                                + node
                                * (
                                    aux_b[4 * objects_count * 3 + j * 3 + k] / 6.0
                                    + node * (
                                        aux_b[5 * objects_count * 3 + j * 3 + k] / 7.0 
                                        + node * aux_b[6 * objects_count * 3 + j * 3 + k] / 8.0
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

// Calculate the auxiliary coefficients b for IAS15
WIN32DLL_API void ias15_compute_aux_b(
    int objects_count,
    int dim_nodes,
    real *restrict aux_b,
    const real *restrict aux_g,
    const real *restrict aux_c,
    int i
)
{
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            if (i >= 1) {
                aux_b[0 * objects_count * 3 + j * 3 + k] = (
                    aux_c[0 * (dim_nodes - 1) + 0] * aux_g[0 * objects_count * 3 + j * 3 + k]
                    + aux_c[1 * (dim_nodes - 1) + 0] * aux_g[1 * objects_count * 3 + j * 3 + k]
                    + aux_c[2 * (dim_nodes - 1) + 0] * aux_g[2 * objects_count * 3 + j * 3 + k]
                    + aux_c[3 * (dim_nodes - 1) + 0] * aux_g[3 * objects_count * 3 + j * 3 + k]
                    + aux_c[4 * (dim_nodes - 1) + 0] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * (dim_nodes - 1) + 0] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * (dim_nodes - 1) + 0] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 2) {
                aux_b[1 * objects_count * 3 + j * 3 + k] = (
                    aux_c[1 * (dim_nodes - 1) + 1] * aux_g[1 * objects_count * 3 + j * 3 + k]
                    + aux_c[2 * (dim_nodes - 1) + 1] * aux_g[2 * objects_count * 3 + j * 3 + k]
                    + aux_c[3 * (dim_nodes - 1) + 1] * aux_g[3 * objects_count * 3 + j * 3 + k]
                    + aux_c[4 * (dim_nodes - 1) + 1] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * (dim_nodes - 1) + 1] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * (dim_nodes - 1) + 1] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 3) {
                aux_b[2 * objects_count * 3 + j * 3 + k] = (
                    aux_c[2 * (dim_nodes - 1) + 2] * aux_g[2 * objects_count * 3 + j * 3 + k]
                    + aux_c[3 * (dim_nodes - 1) + 2] * aux_g[3 * objects_count * 3 + j * 3 + k]
                    + aux_c[4 * (dim_nodes - 1) + 2] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * (dim_nodes - 1) + 2] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * (dim_nodes - 1) + 2] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 4) {
                aux_b[3 * objects_count * 3 + j * 3 + k] = (
                    aux_c[3 * (dim_nodes - 1) + 3] * aux_g[3 * objects_count * 3 + j * 3 + k]
                    + aux_c[4 * (dim_nodes - 1) + 3] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * (dim_nodes - 1) + 3] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * (dim_nodes - 1) + 3] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 5)
            {
                aux_b[4 * objects_count * 3 + j * 3 + k] = (
                    aux_c[4 * (dim_nodes - 1) + 4] * aux_g[4 * objects_count * 3 + j * 3 + k]
                    + aux_c[5 * (dim_nodes - 1) + 4] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * (dim_nodes - 1) + 4] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 6)
            {
                aux_b[5 * objects_count * 3 + j * 3 + k] = (
                    aux_c[5 * (dim_nodes - 1) + 5] * aux_g[5 * objects_count * 3 + j * 3 + k]
                    + aux_c[6 * (dim_nodes - 1) + 5] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }

            if (i >= 7)
            {
                aux_b[6 * objects_count * 3 + j * 3 + k] = (
                    aux_c[6 * (dim_nodes - 1) + 6] * aux_g[6 * objects_count * 3 + j * 3 + k]
                );
            }
            else
            {
                continue;
            }
        }
    }
}

WIN32DLL_API void ias15_compute_aux_g(
    int objects_count,
    int dim_nodes,
    real *restrict aux_g,
    const real *restrict aux_r,
    const real *restrict aux_a,
    int i
)
{
    // Retrieve required accelerations
    real *F1 = malloc(objects_count * 3 * sizeof(real));
    real *F2 = malloc(objects_count * 3 * sizeof(real));
    real *F3 = malloc(objects_count * 3 * sizeof(real));
    real *F4 = malloc(objects_count * 3 * sizeof(real));
    real *F5 = malloc(objects_count * 3 * sizeof(real));
    real *F6 = malloc(objects_count * 3 * sizeof(real));
    real *F7 = malloc(objects_count * 3 * sizeof(real));
    real *F8 = malloc(objects_count * 3 * sizeof(real));

    memcpy(F1, &aux_a[0 * objects_count * 3], objects_count * 3 * sizeof(real));
    memcpy(F2, &aux_a[1 * objects_count * 3], objects_count * 3 * sizeof(real));
    memcpy(F3, &aux_a[2 * objects_count * 3], objects_count * 3 * sizeof(real));
    memcpy(F4, &aux_a[3 * objects_count * 3], objects_count * 3 * sizeof(real));
    memcpy(F5, &aux_a[4 * objects_count * 3], objects_count * 3 * sizeof(real));
    memcpy(F6, &aux_a[5 * objects_count * 3], objects_count * 3 * sizeof(real));
    memcpy(F7, &aux_a[6 * objects_count * 3], objects_count * 3 * sizeof(real));
    memcpy(F8, &aux_a[7 * objects_count * 3], objects_count * 3 * sizeof(real)); 

    // Update aux_g
    for (int j = 0; j < objects_count; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            if (i >= 1)
            {
                aux_g[0 * objects_count * 3 + j * 3 + k] = (
                    (F2[j * 3 + k] - F1[j * 3 + k]) * aux_r[1 * dim_nodes + 0]
                );
            }
            else
            {
                continue;
            }
                
            if (i >= 2)
            {
                aux_g[1 * objects_count * 3 + j * 3 + k] = (
                    ((F3[j * 3 + k] - F1[j * 3 + k]) * aux_r[2 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[2 * dim_nodes + 1]
                );
            }
            else 
            {
                continue;
            }

            if (i >= 3)
            {
                aux_g[2 * objects_count * 3 + j * 3 + k] = (
                    ((F4[j * 3 + k] - F1[j * 3 + k]) * aux_r[3 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[3 * dim_nodes + 1] 
                    - aux_g[1 * objects_count * 3 + j * 3 + k]
                ) * aux_r[3 * dim_nodes + 2];
            }
            else
            {
                continue;
            }
                
            if (i >= 4)
            {
                aux_g[3 * objects_count * 3 + j * 3 + k] = (
                    (((F5[j * 3 + k] - F1[j * 3 + k]) * aux_r[4 * dim_nodes + 0] 
                    - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[4 * dim_nodes + 1] - aux_g[1 * objects_count * 3 + j * 3 + k])
                    * aux_r[4 * dim_nodes + 2]
                    - aux_g[2 * objects_count * 3 + j * 3 + k]
                ) * aux_r[4 * dim_nodes + 3];
            }
            else
            {
                continue;
            }

            if (i >= 5)
            {
                aux_g[4 * objects_count * 3 + j * 3 + k] = (
                    (
                        (((F6[j * 3 + k] - F1[j * 3 + k]) * aux_r[5 * dim_nodes + 0] 
                        - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[5 * dim_nodes + 1] 
                        - aux_g[1 * objects_count * 3 + j * 3 + k])
                        * aux_r[5 * dim_nodes + 2]
                        - aux_g[2 * objects_count * 3 + j * 3 + k]
                    )
                    * aux_r[5 * dim_nodes + 3]
                    - aux_g[3 * objects_count * 3 + j * 3 + k]
                ) * aux_r[5 * dim_nodes + 4];
            }
            else
            {
                continue;
            }

            if (i >= 6)
            {
                aux_g[5 * objects_count * 3 + j * 3 + k] = (
                    (
                        (
                            (((F7[j * 3 + k] - F1[j * 3 + k]) * aux_r[6 * dim_nodes + 0] 
                            - aux_g[0 * objects_count * 3 + j * 3 + k]) * aux_r[6 * dim_nodes + 1] 
                            - aux_g[1 * objects_count * 3 + j * 3 + k])
                            * aux_r[6 * dim_nodes + 2]
                            - aux_g[2 * objects_count * 3 + j * 3 + k]
                        )
                        * aux_r[6 * dim_nodes + 3]
                        - aux_g[3 * objects_count * 3 + j * 3 + k]
                    )
                    * aux_r[6 * dim_nodes + 4]
                    - aux_g[4 * objects_count * 3 + j * 3 + k]
                ) * aux_r[6 * dim_nodes + 5];
            }
            else
            {
                continue;
            }

            if (i >= 7)
            {
                aux_g[6 * objects_count * 3 + j * 3 + k] = (
                    (
                        (
                            (
                                (((F8[j * 3 + k] - F1[j * 3 + k]) * aux_r[7 * dim_nodes + 0] - aux_g[0 * objects_count * 3 + j * 3 + k]) 
                                * aux_r[7 * dim_nodes + 1] 
                                - aux_g[1 * objects_count * 3 + j * 3 + k])
                                * aux_r[7 * dim_nodes + 2]
                                - aux_g[2 * objects_count * 3 + j * 3 + k]
                            )
                            * aux_r[7 * dim_nodes + 3]
                            - aux_g[3 * objects_count * 3 + j * 3 + k]
                        )
                        * aux_r[7 * dim_nodes + 4]
                        - aux_g[4 * objects_count * 3 + j * 3 + k]
                    )
                    * aux_r[7 * dim_nodes + 5]
                    - aux_g[5 * objects_count * 3 + j * 3 + k]
                ) * aux_r[7 * dim_nodes + 6];                
            }
            else
            {
                continue;
            }
        }
    }
    free(F1);
    free(F2);
    free(F3);
    free(F4);
    free(F5);
    free(F6);
    free(F7);
    free(F8);
}

WIN32DLL_API void ias15_refine_aux_b(
    int objects_count,
    int dim_nodes,
    real *restrict aux_b,
    real *restrict aux_e,
    real dt,
    real dt_new,
    int ias15_refine_flag
)
{
    real *delta_aux_b = calloc((dim_nodes - 1) * objects_count * 3, sizeof(real));
    if (ias15_refine_flag != 0)
    {
        for (int i = 0; i < (dim_nodes - 1); i++)
        {
            for (int j = 0; j < objects_count; j++)
            {
                delta_aux_b[i * objects_count * 3 + j * 3 + 0] = (
                    aux_b[i * objects_count * 3 + j * 3 + 0] - aux_e[i * objects_count * 3 + j * 3 + 0]
                );
                delta_aux_b[i * objects_count * 3 + j * 3 + 1] = (
                    aux_b[i * objects_count * 3 + j * 3 + 1] - aux_e[i * objects_count * 3 + j * 3 + 1]
                );
                delta_aux_b[i * objects_count * 3 + j * 3 + 2] = (
                    aux_b[i * objects_count * 3 + j * 3 + 2] - aux_e[i * objects_count * 3 + j * 3 + 2]
                );
            }
        }
    }

    real q = dt_new / dt;
    real q2 = q * q, q3 = q2 * q, q4 = q3 * q, q5 = q4 * q, q6 = q5 * q, q7 = q6 * q;
    
    for (int j = 0; j < objects_count; j++) 
    {
        for (int k = 0; k < 3; k++) 
        {
            aux_e[0 * objects_count * 3 + j * 3 + k] = q * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 7.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 6.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] * 5.0
                + aux_b[3 * objects_count * 3 + j * 3 + k] * 4.0
                + aux_b[2 * objects_count * 3 + j * 3 + k] * 3.0
                + aux_b[1 * objects_count * 3 + j * 3 + k] * 2.0
                + aux_b[0 * objects_count * 3 + j * 3 + k]
            );

            aux_e[1 * objects_count * 3 + j * 3 + k] = q2 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 21.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 15.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] * 10.0
                + aux_b[3 * objects_count * 3 + j * 3 + k] * 6.0
                + aux_b[2 * objects_count * 3 + j * 3 + k] * 3.0
                + aux_b[1 * objects_count * 3 + j * 3 + k]
            );

            aux_e[2 * objects_count * 3 + j * 3 + k] = q3 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 35.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 20.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] * 10.0
                + aux_b[3 * objects_count * 3 + j * 3 + k] * 4.0
                + aux_b[2 * objects_count * 3 + j * 3 + k]
            );

            aux_e[3 * objects_count * 3 + j * 3 + k] = q4 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 35.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 15.0
                + aux_b[4 * objects_count * 3 + j * 3 + k] * 5.0
                + aux_b[3 * objects_count * 3 + j * 3 + k]
            );

            aux_e[4 * objects_count * 3 + j * 3 + k] = q5 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 21.0
                + aux_b[5 * objects_count * 3 + j * 3 + k] * 6.0
                + aux_b[4 * objects_count * 3 + j * 3 + k]
            );

            aux_e[5 * objects_count * 3 + j * 3 + k] = q6 * (
                aux_b[6 * objects_count * 3 + j * 3 + k] * 7.0
                + aux_b[5 * objects_count * 3 + j * 3 + k]
            );

            aux_e[6 * objects_count * 3 + j * 3 + k] = q7 * aux_b[6 * objects_count * 3 + j * 3 + k];
        }
    }

    for (int i = 0; i < (dim_nodes - 1); i++)
    {
        for (int j = 0; j < objects_count; j++)
        {
            aux_b[i * objects_count * 3 + j * 3 + 0] = (
                    aux_e[i * objects_count * 3 + j * 3 + 0] + delta_aux_b[i * objects_count * 3 + j * 3 + 0]
                );
            aux_b[i * objects_count * 3 + j * 3 + 1] = (
                aux_e[i * objects_count * 3 + j * 3 + 1] + delta_aux_b[i * objects_count * 3 + j * 3 + 1]
            );
            aux_b[i * objects_count * 3 + j * 3 + 2] = (
                aux_e[i * objects_count * 3 + j * 3 + 2] + delta_aux_b[i * objects_count * 3 + j * 3 + 2]
            );
        }
    }

    free(delta_aux_b);
}