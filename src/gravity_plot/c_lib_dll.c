#include <math.h>
#include <stdlib.h>
#include <string.h>
// #include <stdio.h> // For testing

#define real double

#ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
    #define min(a,b) ((a) < (b) ? (a) : (b))
#endif

real vec_norm(const real *vec, int vec_size);
void acceleration(int objects_count, const real (*x)[3], real (*a)[3], const real *m, real G);
void euler(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real dt, 
    real tf, 
    int npts, 
    const real *m, 
    real G, 
    real (*sol_state)[6 * objects_count]
);
void euler_cromer(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real dt, 
    real tf, 
    int npts, 
    const real *m, 
    real G, 
    real (*sol_state)[6 * objects_count]
);
void rk4(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real dt, 
    real tf, 
    int npts, 
    const real *m, 
    real G, 
    real (*sol_state)[6 * objects_count]
);
void leapfrog(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real dt, 
    real tf, 
    int npts, 
    const real *m, 
    real G, 
    real (*sol_state)[6 * objects_count]
);
int rk_embedded(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real *dt, 
    real tf, 
    int *count, 
    const real *m, 
    real G, 
    int power,
    int power_test,
    int len_coeff,
    const real (*coeff)[len_coeff],
    int len_weights,
    const real *weights,
    const real *weights_test,
    real abs_tolerance,
    real rel_tolerance,
    real (*sol_state)[6 * objects_count],
    int len_sol_time,
    real *sol_time
);


__declspec(dllexport) real vec_norm(const real *vec, int vec_size)
{   
    real sum = 0.0;
    if (vec_size == 3) 
        sum = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
    else
    {
        for (int i = 0; i < vec_size; i++) sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

__declspec(dllexport) void acceleration(int objects_count, const real (*x)[3], real (*a)[3], const real *m, real G)
{   
    real R_norm, temp_value, *temp_vec = (real *) malloc(3 * sizeof(real)), *R = (real *) malloc(3 * sizeof(real));

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

    free(temp_vec);
    free(R);
}

__declspec(dllexport) void euler(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real dt, 
    real tf, 
    int npts, 
    const real *m, 
    real G, 
    real (*sol_state)[6 * objects_count]
)
{   
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*t / tf * 100);

    // Main Loop
    for(int count = (int) round(*t / dt); count < npts; count++)
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
            sol_state[count + 1][j * 3] = x[j][0];
            sol_state[count + 1][j * 3 + 1] = x[j][1];
            sol_state[count + 1][j * 3 + 2] = x[j][2];
            sol_state[count + 1][objects_count * 3 + j * 3] = v[j][0];
            sol_state[count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
            sol_state[count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
        }    
        *t += dt;

        // Exit to update progress bar
        if ((int) (*t / tf * 100) > progress_percentage)
        {   
            free(a);
            return;
        }
    }
}

__declspec(dllexport) void euler_cromer(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real dt, 
    real tf, 
    int npts, 
    const real *m, 
    real G, 
    real (*sol_state)[6 * objects_count]
)
{   
    real (*a)[3] = malloc(objects_count * 3 * sizeof(real));

    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*t / tf * 100);

    // Main Loop
    for(int count = (int) round(*t / dt); count < npts; count++)
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
            sol_state[count + 1][j * 3] = x[j][0];
            sol_state[count + 1][j * 3 + 1] = x[j][1];
            sol_state[count + 1][j * 3 + 2] = x[j][2];
            sol_state[count + 1][objects_count * 3 + j * 3] = v[j][0];
            sol_state[count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
            sol_state[count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
        }    
        *t += dt;

        // Exit to update progress bar
        if ((int) (*t / tf * 100) > progress_percentage)
        {   
            free(a);
            return;
        }
    }
}

__declspec(dllexport) void rk4(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real dt, 
    real tf, 
    int npts, 
    const real *m, 
    real G, 
    real (*sol_state)[6 * objects_count]
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
    int progress_percentage = (int) round(*t / tf * 100);

    // Main Loop
    for(int count = (int) round(*t / dt); count < npts; count++)
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
            sol_state[count + 1][j * 3] = x[j][0];
            sol_state[count + 1][j * 3 + 1] = x[j][1];
            sol_state[count + 1][j * 3 + 2] = x[j][2];
            sol_state[count + 1][objects_count * 3 + j * 3] = v[j][0];
            sol_state[count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
            sol_state[count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
        }    
        *t += dt;

        // Exit to update progress bar
        if ((int) (*t / tf * 100) > progress_percentage)
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

__declspec(dllexport) void leapfrog(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real dt, 
    real tf, 
    int npts, 
    const real *m, 
    real G, 
    real (*sol_state)[6 * objects_count]
)
{   
    real (*a_0)[3] = malloc(objects_count * 3 * sizeof(real));
    real (*a_1)[3] = malloc(objects_count * 3 * sizeof(real));

    int is_initialize = 1;

    // Current progress percentage rounded to int
    int progress_percentage = (int) round(*t / tf * 100);

    // Main Loop
    for(int count = (int) round(*t / dt); count < npts; count++)
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
            sol_state[count + 1][j * 3] = x[j][0];
            sol_state[count + 1][j * 3 + 1] = x[j][1];
            sol_state[count + 1][j * 3 + 2] = x[j][2];
            sol_state[count + 1][objects_count * 3 + j * 3] = v[j][0];
            sol_state[count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
            sol_state[count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
        }    
        *t += dt;

        // Exit to update progress bar
        if ((int) (*t / tf * 100) > progress_percentage)
        {   
            free(a_0);
            free(a_1);
            return;
        }
    }
}

__declspec(dllexport) int rk_embedded(
    int objects_count, 
    real (*x)[3], 
    real (*v)[3], 
    real *t, 
    real *dt, 
    real tf, 
    int *count, 
    const real *m, 
    real G, 
    int power,
    int power_test,
    int len_coeff,
    const real (*coeff)[len_coeff],
    int len_weights,
    const real *weights,
    const real *weights_test,
    real abs_tolerance,
    real rel_tolerance,
    real (*sol_state)[6 * objects_count],
    int len_sol_time,
    real *sol_time
)
{
    // Initialization
    int stages = len_weights;
    int min_power = min(power, power_test);

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
            tolerance_scale_v[k][0] = abs_tolerance + max(fabs(v[k][0]), fabs(v_1[k][0])) * rel_tolerance;
            tolerance_scale_v[k][1] = abs_tolerance + max(fabs(v[k][1]), fabs(v_1[k][1])) * rel_tolerance;
            tolerance_scale_v[k][2] = abs_tolerance + max(fabs(v[k][2]), fabs(v_1[k][2])) * rel_tolerance;
            tolerance_scale_x[k][0] = abs_tolerance + max(fabs(v[k][0]), fabs(v_1[k][0])) * rel_tolerance;
            tolerance_scale_x[k][1] = abs_tolerance + max(fabs(v[k][1]), fabs(v_1[k][1])) * rel_tolerance;
            tolerance_scale_x[k][2] = abs_tolerance + max(fabs(v[k][2]), fabs(v_1[k][2])) * rel_tolerance;
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
            sol_time[*count] = *t;
            for (int j = 0; j < objects_count; j++)
            {
                sol_state[*count][j * 3] = x[j][0];
                sol_state[*count][j * 3 + 1] = x[j][1];
                sol_state[*count][j * 3 + 2] = x[j][2];
                sol_state[*count][objects_count * 3 + j * 3] = v[j][0];
                sol_state[*count][objects_count * 3 + j * 3 + 1] = v[j][1];
                sol_state[*count][objects_count * 3 + j * 3 + 2] = v[j][2];
            }  

            // Check buffer size and quit if full
            if ((*count + 1) == len_sol_time)
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
        }

        // Calculate dt
        dt_new = *dt * safety_fac / pow(error, (1.0 / (1.0 + (real) min_power)));
        if (dt_new > safety_fac_max * *dt) 
        {
            *dt *= safety_fac_max;
        }
        else if (dt_new < safety_fac_min * *dt)
        {
            *dt *= safety_fac_min;
        }
        else if (dt_new / tf < 1e-12)
        {
            *dt = tf * 1e-12;
        }
        else
        {
            *dt = dt_new;
        }

        // Correct overshooting
        if (*t + *dt > tf)
        {
            *dt = tf - *t;
        }
        
        // Exit to update progress bar
        if ((int) (*t / tf * 100) > progress_percentage)
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
            return 0;
        }
    }
}