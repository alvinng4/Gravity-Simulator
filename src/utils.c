/**
 * \file utils.c
 * \brief Definitions of utility functions
 * 
 * \author Ching-Yin Ng
 * \date March 2025
 */

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>

#include "common.h"
#include "error.h"
#include "math_functions.h"
#include "pcg_basic.h"
#include "system.h"
#include "utils.h"


WIN32DLL_API double grav_get_current_time(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + (double) ts.tv_nsec / 1.0e9;
}

WIN32DLL_API void free_memory_double(double *__restrict ptr)
{
    free(ptr);
}

WIN32DLL_API double compute_energy(const System *__restrict system)
{
    const double *__restrict x = system->x;
    const double *__restrict v = system->v;
    const double *__restrict m = system->m;
    const double G = system->G;
    const int num_particles = system->num_particles;

    double energy = 0.0;

    for (int i = 0; i < num_particles; i++)
    {
        // KE
        double v_norm = vec_norm_3d(&v[i * 3]);
        energy += (
            0.5 * m[i] 
            * v_norm * v_norm
        );

        // PE
        for (int j = i + 1; j < num_particles; j++)
        {
            double r_ij[3];
            for (int k = 0; k < 3; k++)
            {
                r_ij[k] = (x[i * 3 + k] - x[j * 3 + k]);
            }
            energy -= (G * m[i] * m[j] / vec_norm_3d(r_ij));
        }
    }

    return energy;
}

WIN32DLL_API pcg32_random_t init_pcg_rng(void)
{
    pcg32_random_t rng;
    pcg32_srandom_r(&rng, time(NULL), (intptr_t)&rng);

    return rng;
}

WIN32DLL_API double grav_randrange(
    const double min,
    const double max,
    pcg32_random_t *rng
)
{
    // Generate a random number in the range [0, 1)
    const double random_number = ldexp(pcg32_random_r(rng), -32);

    // Scale and shift to the desired range
    return min + (max - min) * random_number;
}

// WIN32DLL_API void compute_energy_python(
//     const int num_particles,
//     const double *__restrict m,
//     const double G,
//     const int npts,
//     int *__restrict count,
//     double *__restrict energy,
//     const double (*__restrict sol_state)[num_particles * 6],
//     int *__restrict is_exit
// )
// {
//     while (*count < npts)
//     {   
//         for (int i = 0; i < num_particles; i++)
//         {
//             // KE
//             double v_norm = vec_norm_3d(&sol_state[*count][(num_particles + i) * 3]);
//             energy[*count] += (
//                 0.5 * m[i] 
//                 * v_norm * v_norm
//             );

//             // PE
//             for (int j = i + 1; j < num_particles; j++)
//             {
//                 double r_ij[3];
//                 for (int k = 0; k < 3; k++)
//                 {
//                     r_ij[k] = (
//                         sol_state[*count][i * 3 + k]
//                         - sol_state[*count][j * 3 + k]
//                     );
//                 }
//                 energy[*count] -= (
//                     G * m[i] * m[j]
//                     / vec_norm_3d(r_ij)
//                 );
//             }
//         }
//         *count += 1;

//         // Check if user sends KeyboardInterrupt in main thread
//         if (*is_exit)
//         {
//             return;
//         }
//     }
// }

// WIN32DLL_API void compute_linear_momentum_python(
//     const int num_particles,
//     const double *__restrict m,
//     const int npts,
//     int *__restrict count,
//     double *__restrict linear_momentum,
//     const double (*__restrict sol_state)[num_particles * 6],
//     int *__restrict is_exit
// )
// {
//     while (*count < npts)
//     {
//         double linear_momentum_vec_step[3] = {0.0};
//         for (int i = 0; i < num_particles; i++)
//         {
//             for (int j = 0; j < 3; j++)
//             {
//                 linear_momentum_vec_step[j] += m[i] * (sol_state[*count][(num_particles + i) * 3 + j]);
//             }
//         }
//         linear_momentum[*count] = vec_norm_3d(linear_momentum_vec_step);
//         *count += 1;

//         // Check if user sends KeyboardInterrupt in main thread
//         if (*is_exit)
//         {
//             return;
//         }
//     }
// }

// WIN32DLL_API void compute_angular_momentum_python(
//     const int num_particles,
//     const double *__restrict m,
//     const int npts,
//     int *__restrict count,
//     double *__restrict angular_momentum,
//     const double (*__restrict sol_state)[num_particles * 6],
//     int *__restrict is_exit
// )
// {
//     while (*count < npts)
//     {
//         double angular_momentum_vec_step[3] = {0.0};
//         // L = m * r x v
//         for (int i = 0; i < num_particles; i++)
//         {
//             angular_momentum_vec_step[0] += m[i] * (
//                 sol_state[*count][i * 3 + 1] 
//                 * sol_state[*count][(num_particles + i) * 3 + 2]
//                 - sol_state[*count][i * 3 + 2] 
//                 * sol_state[*count][(num_particles + i) * 3 + 1]
//             );
//             angular_momentum_vec_step[1] += m[i] * (
//                 sol_state[*count][i * 3 + 2]
//                 * sol_state[*count][(num_particles + i) * 3 + 0]
//                 - sol_state[*count][i * 3 + 0]
//                 * sol_state[*count][(num_particles + i) * 3 + 2]
//             );
//             angular_momentum_vec_step[2] += m[i] * (
//                 sol_state[*count][i * 3]
//                 * sol_state[*count][(num_particles + i) * 3 + 1]
//                 - sol_state[*count][i * 3 + 1]
//                 * sol_state[*count][(num_particles + i) * 3]
//             );
//         }
//         angular_momentum[*count] = vec_norm_3d(angular_momentum_vec_step);
//         *count += 1;

//         // Check if user sends KeyboardInterrupt in main thread
//         if (*is_exit)
//         {
//             return;
//         }
//     }
// }
