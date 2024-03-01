#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h> // For testing

#define real double

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
real vec_norm(const real *vec, int vec_size);

void acceleration(int objects_count, const real (*x)[3], real (*a)[3], const real *m, real G)
{   
    real R_norm, temp_value, *temp_vec = (real *) malloc(3 * sizeof(real)), *R = (real *) malloc(3 * sizeof(real));

    memset(a, 0, 3 * objects_count * sizeof(real));

    for(int i = 0; i < objects_count; i++)
    {
        for(int j = i + 1; j < objects_count; j++)
        {
            R[0] = x[i][0] - x[j][0];
            R[1] = x[i][1] - x[j][1];
            R[2] = x[i][2] - x[j][2];
            R_norm = vec_norm(R, 3);
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
)
{   
    real (*a)[3] = malloc(3 * objects_count * sizeof(real));

    int progress_percentage = (int) round(*t / tf * 100);
    for(int count = (int) round(*t / dt); count < npts; count++)
    {   
        acceleration(objects_count, x, a, m, G);
        for (int j = 0; j < objects_count; j++)
        {
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
        if ((int) (*t / tf * 100) > progress_percentage)
        {   
            free(a);
            return;
        }
    }
}

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
)
{   
    real (*a)[3] = malloc(3 * objects_count * sizeof(real));

    int progress_percentage = (int) round(*t / tf * 100);
    for(int count = (int) round(*t / dt); count < npts; count++)
    {   
        acceleration(objects_count, x, a, m, G);
        for (int j = 0; j < objects_count; j++)
        {
            v[j][0] += a[j][0] * dt;
            v[j][1] += a[j][1] * dt;
            v[j][2] += a[j][2] * dt;
            x[j][0] += v[j][0] * dt;
            x[j][1] += v[j][1] * dt;
            x[j][2] += v[j][2] * dt;
            sol_state[count + 1][j * 3] = x[j][0];
            sol_state[count + 1][j * 3 + 1] = x[j][1];
            sol_state[count + 1][j * 3 + 2] = x[j][2];
            sol_state[count + 1][objects_count * 3 + j * 3] = v[j][0];
            sol_state[count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
            sol_state[count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
        }    
            
        *t += dt;
        if ((int) (*t / tf * 100) > progress_percentage)
        {   
            free(a);
            return;
        }
    }
}

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
)
{
    real (*temp_x)[3] = malloc(3 * objects_count * sizeof(real));
    real (*temp_v)[3] = malloc(3 * objects_count * sizeof(real));

    real (*a)[3] = malloc(3 * objects_count * sizeof(real));

    real (*vk1)[3] = malloc(3 * objects_count * sizeof(real));
    real (*vk2)[3] = malloc(3 * objects_count * sizeof(real));
    real (*vk3)[3] = malloc(3 * objects_count * sizeof(real));
    real (*vk4)[3] = malloc(3 * objects_count * sizeof(real));
    real (*xk1)[3] = malloc(3 * objects_count * sizeof(real));
    real (*xk2)[3] = malloc(3 * objects_count * sizeof(real));
    real (*xk3)[3] = malloc(3 * objects_count * sizeof(real));
    real (*xk4)[3] = malloc(3 * objects_count * sizeof(real));

    int progress_percentage = (int) round(*t / tf * 100);
    for(int count = (int) round(*t / dt); count < npts; count++)
    {   
        acceleration(objects_count, x, a, m, G);
        memcpy(vk1, a, 3 * objects_count * sizeof(real));
        memcpy(xk1, v, 3 * objects_count * sizeof(real));

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
        memcpy(vk2, a, 3 * objects_count * sizeof(real));
        memcpy(xk2, temp_v, 3 * objects_count * sizeof(real));

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
        memcpy(vk3, a, 3 * objects_count * sizeof(real));
        memcpy(xk3, temp_v, 3 * objects_count * sizeof(real));

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
        memcpy(vk4, a, 3 * objects_count * sizeof(real));
        memcpy(xk4, temp_v, 3 * objects_count * sizeof(real));


        for (int j = 0; j < objects_count; j++)
        {
            v[j][0] += (vk1[j][0] + 2 * vk2[j][0] + 2 * vk3[j][0] + vk4[j][0]) * dt / 6.0;
            v[j][1] += (vk1[j][1] + 2 * vk2[j][1] + 2 * vk3[j][1] + vk4[j][1]) * dt / 6.0;
            v[j][2] += (vk1[j][2] + 2 * vk2[j][2] + 2 * vk3[j][2] + vk4[j][2]) * dt / 6.0;
            x[j][0] += (xk1[j][0] + 2 * xk2[j][0] + 2 * xk3[j][0] + xk4[j][0]) * dt / 6.0;
            x[j][1] += (xk1[j][1] + 2 * xk2[j][1] + 2 * xk3[j][1] + xk4[j][1]) * dt / 6.0;
            x[j][2] += (xk1[j][2] + 2 * xk2[j][2] + 2 * xk3[j][2] + xk4[j][2]) * dt / 6.0;
            sol_state[count + 1][j * 3] = x[j][0];
            sol_state[count + 1][j * 3 + 1] = x[j][1];
            sol_state[count + 1][j * 3 + 2] = x[j][2];
            sol_state[count + 1][objects_count * 3 + j * 3] = v[j][0];
            sol_state[count + 1][objects_count * 3 + j * 3 + 1] = v[j][1];
            sol_state[count + 1][objects_count * 3 + j * 3 + 2] = v[j][2];
        }    
            
        *t += dt;
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
)
{   
    real (*a_0)[3] = malloc(3 * objects_count * sizeof(real));
    real (*a_1)[3] = malloc(3 * objects_count * sizeof(real));

    int is_initialize = 1;

    int progress_percentage = (int) round(*t / tf * 100);
    for(int count = (int) round(*t / dt); count < npts; count++)
    {       
        if (is_initialize == 1)
        {
            acceleration(objects_count, x, a_0, m, G);
            is_initialize = 0;
        }
        else 
        {
            memcpy(a_0, a_1, 3 * objects_count * sizeof(real));
        }
           

        for (int j = 0; j < objects_count; j++)
        {
            x[j][0] += v[j][0] * dt + 0.5 * a_0[j][0] * dt * dt;
            x[j][1] += v[j][1] * dt + 0.5 * a_0[j][1] * dt * dt;
            x[j][2] += v[j][2] * dt + 0.5 * a_0[j][2] * dt * dt;
        }    
        acceleration(objects_count, x, a_1, m, G);
        for (int j = 0; j < objects_count; j++)
        {
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
        
        if ((int) (*t / tf * 100) > progress_percentage)
        {   
            free(a_0);
            free(a_1);
            return;
        }
    }
}


real vec_norm(const real *vec, int vec_size)
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
