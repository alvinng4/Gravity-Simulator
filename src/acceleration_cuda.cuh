#ifndef ACCELERATION_CUDA_H
#define ACCELERATION_CUDA_H

/**
 * acceleration_cuda: Acceleration calculation with CUDA
 */

#include <cuda_runtime.h>

#include "gravity_sim.h"

#define BLOCK_SIZE 64

#ifdef __cplusplus
    extern "C" {
#endif
    int acceleration_pairwise_cuda(
        double *__restrict a,
        const System *__restrict system,
        const AccelerationParam *__restrict acceleration_param
    );

    int acceleration_pairwise_cuda_float(
        double *__restrict a,
        const System *__restrict system,
        const AccelerationParam *__restrict acceleration_param
    );
#ifdef __cplusplus
}
#endif


int acceleration_barnes_hut_cuda(
    double *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
);




#endif