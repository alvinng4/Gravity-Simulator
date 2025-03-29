/**
 * \file math_functions.c
 * 
 * \brief This file contains the implementation of some useful mathematical functions
 * 
 * \author Ching-Yin Ng
 * \date March 2025
 */

#include <math.h>

#include "common.h"

WIN32DLL_API double vec_norm_3d(const double *__restrict vec)
{
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

WIN32DLL_API double vec_norm(const double *__restrict vec, const int vec_length)
{   
    double sum = 0.0;
    for (int i = 0; i < vec_length; i++) 
    {
        sum += vec[i] * vec[i];
    }

    return sqrt(sum);
}

WIN32DLL_API double abs_max_vec(const double *__restrict vec, const int vec_length)
{
    double max = fabs(vec[0]);
    for (int i = 1; i < vec_length; i++)
    {
        max = fmax(max, fabs(vec[i]));
    }

    return max;
}

WIN32DLL_API double vec_dot_3d(
    const double *__restrict vec_1,
    const double *__restrict vec_2
)
{
    return vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2];
}

WIN32DLL_API double vec_dot(
    const double *__restrict vec_1,
    const double *__restrict vec_2,
    const int vec_length
)
{
    double sum = 0.0;
    for (int i = 0; i < vec_length; i++)
    {
        sum += vec_1[i] * vec_2[i];
    }

    return sum;
}
