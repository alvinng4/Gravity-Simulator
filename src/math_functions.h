/**
 * \file math_functions.h
 * 
 * \brief Library for useful mathematical functions
 * 
 * \author Ching-Yin Ng
 * \date March 2025
 */

#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

#include "gravity_sim.h"

/**
 * \brief Calculate the norm of a 3D vector
 * 
 * \param vec 1D array of size 3
 * 
 * \return Norm of the vector
 */
double vec_norm_3d(const double *__restrict vec);

/**
 * \brief Calculate the norm of a vector
 * 
 * \param vec 1D array of size vec_length
 * \param vec_length Length of the vector
 * 
 * \return Norm of the vector
 */
double vec_norm(const double *__restrict vec, const int vec_length);

/**
 * \brief Calculate the absolute maximum value from a vector
 * 
 * \param vec 1D array of size vec_length
 * \param vec_length Length of the vector
 * 
 * \return Absolute maximum value from the vector
 */
double abs_max_vec(const double *__restrict vec, const int vec_length);

/**
 * \brief Calculate the dot product of two 3D vectors
 * 
 * \param vec_1 1D array of size 3
 * \param vec_2 1D array of size 3
 * 
 * \return Dot product of the two vectors
 */
double vec_dot_3d(
    const double *__restrict vec_1,
    const double *__restrict vec_2
);

/**
 * \brief Calculate the dot product of two vectors
 * 
 * \param vec_1 1D array of size vec_length
 * \param vec_2 1D array of size vec_length
 * \param vec_length Length of the vectors
 * 
 * \return Dot product of the two vectors
 */
double vec_dot(
    const double *__restrict vec_1,
    const double *__restrict vec_2,
    const int vec_length
);

#endif
