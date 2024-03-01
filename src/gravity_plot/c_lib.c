#include <math.h>
#include <stdlib.h>
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
    real (*sol_state)[6]
);
real vec_norm(const real *vec, int vec_size);

void acceleration(int objects_count, const real (*x)[3], real (*a)[3], const real *m, real G)
{   
    real R_norm, temp_value, *temp_vec = (real *) malloc(3 * sizeof(real)), *R = (real *) malloc(3 * sizeof(real));

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
            a[i][0] = - temp_vec[0] * m[j];
            a[i][1] = - temp_vec[1] * m[j];
            a[i][2] = - temp_vec[2] * m[j];
            a[j][0] = temp_vec[0] * m[i];
            a[j][1] = temp_vec[1] * m[i];
            a[j][2] = temp_vec[2] * m[i];
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
    real (*sol_state)[6]
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
            sol_state[count + 1][j * 6] = x[j][0];
            sol_state[count + 1][j * 6 + 1] = x[j][1];
            sol_state[count + 1][j * 6 + 2] = x[j][2];
            sol_state[count + 1][j * 6 + 3] = v[j][0];
            sol_state[count + 1][j * 6 + 4] = v[j][1];
            sol_state[count + 1][j * 6 + 5] = v[j][2];
        }        
        *t += dt;
        if ((int) (*t / tf * 100) > progress_percentage)
        {
            free(a);
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
