/**
 * \file utils.c
 * \brief Definitions of utility functions
 */

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
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

WIN32DLL_API double compute_energy(const System *restrict system)
{
    const double *restrict x = system->x;
    const double *restrict v = system->v;
    const double *restrict m = system->m;
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

WIN32DLL_API void keplerian_to_cartesian(
    double *restrict x,
    double *restrict v,
    const double semi_major_axis,
    const double eccentricity,
    const double inclination,
    const double argument_of_periapsis,
    const double longitude_of_ascending_node,
    const double true_anomaly,
    const double total_mass,
    const double G
)
{
    const double cos_inc = cos(inclination);
    const double sin_inc = sin(inclination);

    const double cos_arg_periapsis = cos(argument_of_periapsis);
    const double sin_arg_periapsis = sin(argument_of_periapsis);

    const double cos_long_asc_node = cos(longitude_of_ascending_node);
    const double sin_long_asc_node = sin(longitude_of_ascending_node);

    const double cos_true_anomaly = cos(true_anomaly);
    const double sin_true_anomaly = sin(true_anomaly);

    const double ecc_unit_vec[3] = {
        cos_long_asc_node * cos_arg_periapsis
        - sin_long_asc_node * sin_arg_periapsis * cos_inc,
        sin_long_asc_node * cos_arg_periapsis
        + cos_long_asc_node * sin_arg_periapsis * cos_inc,
        sin_arg_periapsis * sin_inc
    };

    const double q_unit_vec[3] = {
        -cos_long_asc_node * sin_arg_periapsis
        - sin_long_asc_node * cos_arg_periapsis * cos_inc,
        -sin_long_asc_node * sin_arg_periapsis
        + cos_long_asc_node * cos_arg_periapsis * cos_inc,
        cos_arg_periapsis * sin_inc
    };

    for (int i = 0; i < 3; i++)
    {
        x[i] = semi_major_axis * (
            (1.0 - eccentricity * eccentricity)
            * (cos_true_anomaly * ecc_unit_vec[i]+ sin_true_anomaly * q_unit_vec[i])
            / (1.0 + eccentricity * cos_true_anomaly)
        );
    }
    for (int i = 0; i < 3; i++)
    {
        v[i] = (
            sqrt(G * total_mass / (semi_major_axis * (1.0 - eccentricity * eccentricity)))
            * (
                -sin_true_anomaly * ecc_unit_vec[i]
                + (eccentricity + cos_true_anomaly) * q_unit_vec[i]
            )
        );
    }
}
