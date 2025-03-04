#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "acceleration_cuda.cuh"
#include "error.h"
#include "gravity_sim.h"

__global__ void acceleration_pairwise_kernel(
    double *__restrict a,
    const int objects_count,
    const double *__restrict x,
    const double *__restrict m,
    const double G,
    const double softening_length
)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= objects_count)
    {
        return;
    }

    double3 local_a = make_double3(0.0, 0.0, 0.0);
    const double3 x_i = make_double3(
        x[i * 3 + 0],
        x[i * 3 + 1],
        x[i * 3 + 2]
    );

    for (int j = 0; j < objects_count; j++)
    {
        if (i == j)
        {
            continue;
        }

        double3 dx = make_double3(
            x_i.x - x[j * 3 + 0],
            x_i.y - x[j * 3 + 1],
            x_i.z - x[j * 3 + 2]
        );
        const double r_norm = sqrt(
            dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + softening_length * softening_length
        );
        const double temp_value = G * m[j] / (r_norm * r_norm * r_norm);
        local_a.x -= temp_value * dx.x;
        local_a.y -= temp_value * dx.y;
        local_a.z -= temp_value * dx.z;
    }
    a[i * 3 + 0] = local_a.x;
    a[i * 3 + 1] = local_a.y;
    a[i * 3 + 2] = local_a.z;
}

extern "C"
{
    WIN32DLL_API int acceleration_pairwise_cuda(
        double *__restrict a,
        const System *__restrict system,
        const AccelerationParam *__restrict acceleration_param
    )
    {
        int return_code;

        const int objects_count = system->objects_count;
        const double *x = system->x;
        const double *m = system->m;
        const double G = system->G;
        const double softening_length = acceleration_param->softening_length;

        double *a_device = NULL;
        double *x_device = NULL;
        double *m_device = NULL;
        cudaError_t error;

        error = cudaMalloc((double **) &a_device, objects_count * 3 * sizeof(double));
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMORY_ALLOC;
            goto err_gpu_memory;
        }
        error = cudaMalloc((double **) &x_device, objects_count * 3 * sizeof(double));
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMORY_ALLOC;
            goto err_gpu_memory;
        }
        error = cudaMalloc((double **) &m_device, objects_count * sizeof(double));
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMORY_ALLOC;
            goto err_gpu_memory;
        }

        error = cudaMemcpy(x_device, x, objects_count * 3 * sizeof(double), cudaMemcpyHostToDevice);
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMCPY_CPU_TO_GPU;
            goto err_memcpy_cpu_to_gpu;
        }
        error = cudaMemcpy(m_device, m, objects_count * sizeof(double), cudaMemcpyHostToDevice);
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMCPY_CPU_TO_GPU;
            goto err_memcpy_cpu_to_gpu;
        }

        acceleration_pairwise_kernel <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>>(
            a_device,
            objects_count,
            x_device,
            m_device,
            G,
            softening_length
        );

        error = cudaMemcpy(a, a_device, objects_count * 3 * sizeof(double), cudaMemcpyDeviceToHost);
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMCPY_GPU_TO_CPU;
            goto err_memcpy_gpu_to_cpu;
        }

        cudaFree(a_device);
        cudaFree(x_device);
        cudaFree(m_device);

        return SUCCESS;

    err_memcpy_gpu_to_cpu:
    err_memcpy_cpu_to_gpu:
    err_gpu_memory:
        cudaFree(a_device);
        cudaFree(x_device);
        cudaFree(m_device);
        return return_code;
    }
}

__global__ void memcpy_array_double_to_float(
    const double *__restrict x_double,
    const double *__restrict m_double,
    const int objects_count,
    float *__restrict x,
    float *__restrict m
)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= objects_count)
    {
        return;
    }

    x[i * 3 + 0] = x_double[i * 3 + 0];
    x[i * 3 + 1] = x_double[i * 3 + 1];
    x[i * 3 + 2] = x_double[i * 3 + 2];
    m[i] = m_double[i];

    return;
}

__global__ void acceleration_pairwise_float_kernel(
    double *__restrict a_double,
    const int objects_count,
    const float *__restrict x,
    const float *__restrict m,
    const float G,
    const float softening_length
)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= objects_count)
    {
        return;
    }

    float3 local_a = make_float3(0.0, 0.0, 0.0);
    const float3 x_i = make_float3(
        x[i * 3 + 0],
        x[i * 3 + 1],
        x[i * 3 + 2]
    );

    for (int j = 0; j < objects_count; j++)
    {
        if (i == j)
        {
            continue;
        }

        float3 dx = make_float3(
            x_i.x - x[j * 3 + 0],
            x_i.y - x[j * 3 + 1],
            x_i.z - x[j * 3 + 2]
        );
        const float r_norm = sqrt(
            dx.x * dx.x + dx.y * dx.y + dx.z * dx.z + softening_length * softening_length
        );
        const float temp_value = G * m[j] / (r_norm * r_norm * r_norm);
        local_a.x -= temp_value * dx.x;
        local_a.y -= temp_value * dx.y;
        local_a.z -= temp_value * dx.z;
    }
    a_double[i * 3 + 0] = local_a.x;
    a_double[i * 3 + 1] = local_a.y;
    a_double[i * 3 + 2] = local_a.z;
}

extern "C"
{
    WIN32DLL_API int acceleration_pairwise_cuda_float(
        double *__restrict a,
        const System *__restrict system,
        const AccelerationParam *__restrict acceleration_param
    )
    {
        int return_code;

        const int objects_count = system->objects_count;
        const double *x = system->x;
        const double *m = system->m;
        const double G = system->G;
        const double softening_length = acceleration_param->softening_length;

        double *a_double_device = NULL;
        double *x_double_device = NULL;
        double *m_double_device = NULL;
        float *x_device = NULL;
        float *m_device = NULL;
        cudaError_t error;

        error = cudaMalloc((double **) &a_double_device, objects_count * 3 * sizeof(double));
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMORY_ALLOC;
            goto err_gpu_memory;
        }
        error = cudaMalloc((double **) &x_double_device, objects_count * 3 * sizeof(double));
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMORY_ALLOC;
            goto err_gpu_memory;
        }
        error = cudaMalloc((double **) &m_double_device, objects_count * sizeof(double));
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMORY_ALLOC;
            goto err_gpu_memory;
        }
        error = cudaMalloc((float **) &x_device, objects_count * 3 * sizeof(float));
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMORY_ALLOC;
            goto err_gpu_memory;
        }
        error = cudaMalloc((float **) &m_device, objects_count * sizeof(float));
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMORY_ALLOC;
            goto err_gpu_memory;
        }

        error = cudaMemcpy(x_double_device, x, objects_count * 3 * sizeof(double), cudaMemcpyHostToDevice);
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMCPY_CPU_TO_GPU;
            goto err_memcpy_cpu_to_gpu;
        }
        error = cudaMemcpy(m_double_device, m, objects_count * sizeof(double), cudaMemcpyHostToDevice);
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMCPY_CPU_TO_GPU;
            goto err_memcpy_cpu_to_gpu;
        }

        memcpy_array_double_to_float <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>>(
            x_double_device,
            m_double_device,
            objects_count,
            x_device,
            m_device
        );

        acceleration_pairwise_float_kernel <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>>(
            a_double_device,
            objects_count,
            x_device,
            m_device,
            G,
            softening_length
        );

        error = cudaMemcpy(a, a_double_device, objects_count * 3 * sizeof(double), cudaMemcpyDeviceToHost);
        if (error != cudaSuccess)
        {
            fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
            return_code = ERROR_CUDA_PAIRWISE_MEMCPY_GPU_TO_CPU;
            goto err_memcpy_gpu_to_cpu;
        }

        cudaFree(a_double_device);
        cudaFree(x_double_device);
        cudaFree(m_double_device);
        cudaFree(x_device);
        cudaFree(m_device);

        return SUCCESS;

    err_memcpy_gpu_to_cpu:
    err_memcpy_cpu_to_gpu:
    err_gpu_memory:
        cudaFree(a_double_device);
        cudaFree(x_double_device);
        cudaFree(m_double_device);
        cudaFree(x_device);
        cudaFree(m_device);
        return return_code;
    }
}
