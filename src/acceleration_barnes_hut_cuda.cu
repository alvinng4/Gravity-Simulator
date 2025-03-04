#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration_barnes_hut.h"
#include "acceleration_cuda.cuh"
#include "error.h"
#include "gravity_sim.h"
#include "math_functions.h"


/**
 * \brief Check if two Morton indices are included in the same octant
 * 
 * \param morton_index_i Morton index of the first object at the deepest level
 * \param morton_index_j Morton index of the second object at the deepest level
 * \param level Level of the Morton indices
 */
 inline __device__ bool _check_if_included(
    const int64 morton_index_i,
    const int64 morton_index_j,
    const int level
)
{
    return (morton_index_i >> (3 * (MORTON_MAX_LEVEL - level))) == (morton_index_j >> (3 * (MORTON_MAX_LEVEL - level)));
}

/**
 * \brief Helper kernel function for computing the acceleration of one particle
 * 
 * \param a Array of acceleration vectors
 * \param objects_count Number of objects
 * \param x Array of position vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param softening_length Softening length
 * \param opening_angle Opening angle
 * \param width Width of the bounding box
 * \param leaf_morton_indices_deepest_level Array of Morton indices
 * \param sorted_indices Array of sorted indices
 * \param tree_start_particle_sorted_idx Array of start index of particles in the node
 * \param tree_num_particles Array of number of particles in the node
 * \param tree_num_internal_children Array of number of internal children
 * \param tree_idx_first_internal_child Array of index of the first internal child
 * \param tree_total_mass Array of total mass of the node
 * \param tree_center_of_mass_x Array of x-coordinate of the center of mass
 * \param tree_center_of_mass_y Array of y-coordinate of the center of mass
 * \param tree_center_of_mass_z Array of z-coordinate of the center of mass
 * 
 * \retval SUCCESS if successful
 */
__global__ void _acceleration_helper_kernel(
    double *__restrict a,
    const int objects_count,
    const double *__restrict x,
    const double *__restrict m,
    const double G,
    const double softening_length,
    const double opening_angle,
    const double width,
    const int64 *__restrict leaf_morton_indices_deepest_level,
    const int *__restrict sorted_indices,
    const int *__restrict tree_start_particle_sorted_idx,
    const int *__restrict tree_num_particles,
    const int *__restrict tree_num_internal_children,
    const int *__restrict tree_idx_first_internal_child,
    const double *__restrict tree_total_mass,
    const double *__restrict tree_center_of_mass_x,
    const double *__restrict tree_center_of_mass_y,
    const double *__restrict tree_center_of_mass_z
)
{
    typedef struct Stack
    {
        int node;
        int processed_children;
        struct Stack *last;
    } Stack;

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= objects_count)
    {
        return;
    }

    const int idx_i = sorted_indices[i];
    const int64 morton_index_i = leaf_morton_indices_deepest_level[idx_i];
    const double3 x_i = make_double3(x[idx_i * 3 + 0], x[idx_i * 3 + 1], x[idx_i * 3 + 2]);
    const double softening_squared = softening_length * softening_length;
    
    Stack stack_pool[MORTON_MAX_LEVEL];
    Stack *stack = &(stack_pool[0]);
    stack->processed_children = -1;
    stack->last = NULL;
    stack->node = 0;
    double3 local_a = make_double3(0.0, 0.0, 0.0);

    int level = 1;

    /* Tree walk */
    while (true)
    {
        const int node = stack->node;
        for (int j = (stack->processed_children) + 1; j < tree_num_internal_children[node]; j++)
        {
            const int child_j = tree_idx_first_internal_child[node] + j;
            const int num_children_j = tree_num_internal_children[child_j];
            const int start_idx_j = tree_start_particle_sorted_idx[child_j];

            /* Leaf node */
            if (num_children_j <= 0)
            {
                const int num_particles_j = tree_num_particles[child_j];
                for (int k = 0; k < num_particles_j; k++)
                {
                    const int idx_j = sorted_indices[start_idx_j + k];
                    if (idx_i == idx_j)
                    {
                        continue;
                    }

                    // Calculate \vec{R} and its norm
                    const double3 dr = make_double3(
                        x_i.x - x[idx_j * 3 + 0],
                        x_i.y - x[idx_j * 3 + 1],
                        x_i.z - x[idx_j * 3 + 2]
                    );
                    const double R_norm = sqrt(
                        dr.x * dr.x + dr.y * dr.y + dr.z * dr.z + softening_squared
                    );

                    // Calculate the acceleration
                    const double temp_value = G * m[idx_j] / (R_norm * R_norm * R_norm);
                    local_a.x -= temp_value * dr.x;
                    local_a.y -= temp_value * dr.y;
                    local_a.z -= temp_value * dr.z;
                }

                stack->processed_children = j;
                continue;
            }

            /* Internal node */
            else
            {
                bool criteria_met = false;

                // If object i is included, then we need to traverse deeper
                const bool is_included = _check_if_included(
                    morton_index_i,
                    leaf_morton_indices_deepest_level[sorted_indices[start_idx_j]],
                    level
                );

                // Check Barnes-Hut criteria
                double3 dr;
                double norm_square;
                if (!is_included)
                {
                    dr = make_double3(
                        x_i.x - tree_center_of_mass_x[child_j],
                        x_i.y - tree_center_of_mass_y[child_j],
                        x_i.z - tree_center_of_mass_z[child_j]
                    );
                    const double width_j = width / (2 << level);
                    norm_square = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                    if (width_j / sqrt(norm_square) < opening_angle)
                    {
                        criteria_met = true;
                    }
                }

                // Traverse deeper
                if (!criteria_met)
                {
                    Stack *new_item = &(stack_pool[level + 1]);
                    new_item->node = child_j;
                    new_item->last = stack;
                    new_item->processed_children = -1;

                    stack = new_item;
                    level++;
                    break;
                }

                else
                {
                    const double R_norm = sqrt(
                        norm_square + softening_squared
                    );

                    const double temp_value = G / (R_norm * R_norm * R_norm);
                    local_a.x -= temp_value * dr.x * tree_total_mass[child_j];
                    local_a.y -= temp_value * dr.y * tree_total_mass[child_j];
                    local_a.z -= temp_value * dr.z * tree_total_mass[child_j];

                    stack->processed_children = j;
                    continue;
                }
            }
        }

        if ((stack->processed_children + 1) >= tree_num_internal_children[stack->node])
        {
            Stack *parent = stack->last;
            if (!parent)
            {
                break;
            }            
            stack = parent;
            stack->processed_children += 1;
            level--;
        }
    }

    a[idx_i * 3 + 0] = local_a.x;
    a[idx_i * 3 + 1] = local_a.y;
    a[idx_i * 3 + 2] = local_a.z;

    return;
}

/**
 * \brief Compute the acceleration of the particles
 * 
 * \param a Array of acceleration vectors
 * \param objects_count Number of objects
 * \param x Array of position vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param softening_length Softening length
 * \param opening_angle Opening angle
 * \param width Width of the bounding box
 * \param actual_num_internal_nodes Number of internal nodes
 * \param leaf_morton_indices_deepest_level Array of Morton indices
 * \param sorted_indices Array of sorted indices
 * \param tree_start_particle_sorted_idx Array of start index of particles in the node
 * \param tree_num_particles Array of number of particles in the node
 * \param tree_num_internal_children Array of number of internal children
 * \param tree_idx_first_internal_child Array of index of the first internal child
 * \param tree_total_mass Array of total mass of the node
 * \param tree_center_of_mass_x Array of x-coordinate of the center of mass
 * \param tree_center_of_mass_y Array of y-coordinate of the center of mass
 * \param tree_center_of_mass_z Array of z-coordinate of the center of mass
 * 
 * \retval SUCCESS if successful
 */
IN_FILE int _compute_acceleration(
    double *__restrict a,
    const int objects_count,
    const double *__restrict x,
    const double *__restrict m,
    const double G,
    const double softening_length,
    const double opening_angle,
    const double width,
    const int actual_num_internal_nodes,
    const int64 *__restrict leaf_morton_indices_deepest_level,
    const int *__restrict sorted_indices,
    const int *__restrict tree_start_particle_sorted_idx,
    const int *__restrict tree_num_particles,
    const int *__restrict tree_num_internal_children,
    const int *__restrict tree_idx_first_internal_child,
    const double *__restrict tree_total_mass,
    const double *__restrict tree_center_of_mass_x,
    const double *__restrict tree_center_of_mass_y,
    const double *__restrict tree_center_of_mass_z
)
{
    int return_code;

    double *__restrict a_device = NULL;
    double *__restrict x_device = NULL;
    double *__restrict m_device = NULL;
    int64 *__restrict leaf_morton_indices_deepest_level_device = NULL;
    int *__restrict sorted_indices_device = NULL;
    int *__restrict tree_start_particle_sorted_idx_device = NULL;
    int *__restrict tree_num_particles_device = NULL;
    int *__restrict tree_num_internal_children_device = NULL;
    int *__restrict tree_idx_first_internal_child_device = NULL;
    double *__restrict tree_total_mass_device = NULL;
    double *__restrict tree_center_of_mass_x_device = NULL;
    double *__restrict tree_center_of_mass_y_device = NULL;
    double *__restrict tree_center_of_mass_z_device = NULL;
    cudaError_t error;

    /* Allocate memory on GPU */
    error = cudaMalloc((double **) &a_device, objects_count * 3 * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &x_device, objects_count * 3 * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &m_device, objects_count * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int64 **) &leaf_morton_indices_deepest_level_device, objects_count * sizeof(int64));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &sorted_indices_device, objects_count * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &tree_start_particle_sorted_idx_device, actual_num_internal_nodes * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &tree_num_particles_device, actual_num_internal_nodes * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &tree_num_internal_children_device, actual_num_internal_nodes * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &tree_idx_first_internal_child_device, actual_num_internal_nodes * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &tree_total_mass_device, actual_num_internal_nodes * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &tree_center_of_mass_x_device, actual_num_internal_nodes * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &tree_center_of_mass_y_device, actual_num_internal_nodes * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &tree_center_of_mass_z_device, actual_num_internal_nodes * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }

    error = cudaMemcpy(x_device, x, objects_count * 3 * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(m_device, m, objects_count * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(leaf_morton_indices_deepest_level_device, leaf_morton_indices_deepest_level, objects_count * sizeof(int64), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(sorted_indices_device, sorted_indices, objects_count * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_start_particle_sorted_idx_device, tree_start_particle_sorted_idx, actual_num_internal_nodes * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_num_particles_device, tree_num_particles, actual_num_internal_nodes * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_num_internal_children_device, tree_num_internal_children, actual_num_internal_nodes * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_idx_first_internal_child_device, tree_idx_first_internal_child, actual_num_internal_nodes * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_total_mass_device, tree_total_mass, actual_num_internal_nodes * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_center_of_mass_x_device, tree_center_of_mass_x, actual_num_internal_nodes * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_center_of_mass_y_device, tree_center_of_mass_y, actual_num_internal_nodes * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_center_of_mass_z_device, tree_center_of_mass_z, actual_num_internal_nodes * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    } 

    _acceleration_helper_kernel <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>> (
        a_device,
        objects_count,
        x_device,
        m_device,
        G,
        softening_length,
        opening_angle,
        width,
        leaf_morton_indices_deepest_level_device,
        sorted_indices_device,
        tree_start_particle_sorted_idx_device,
        tree_num_particles_device,
        tree_num_internal_children_device,
        tree_idx_first_internal_child_device,
        tree_total_mass_device,
        tree_center_of_mass_x_device,
        tree_center_of_mass_y_device,
        tree_center_of_mass_z_device
    );

    error = cudaMemcpy(a, a_device, objects_count * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_GPU_TO_CPU;
        goto err_memcpy_gpu_to_cpu;
    }

    cudaFree(a_device);
    cudaFree(x_device);
    cudaFree(m_device);
    cudaFree(leaf_morton_indices_deepest_level_device);
    cudaFree(sorted_indices_device);
    cudaFree(tree_start_particle_sorted_idx_device);
    cudaFree(tree_num_particles_device);
    cudaFree(tree_num_internal_children_device);
    cudaFree(tree_idx_first_internal_child_device);
    cudaFree(tree_total_mass_device);
    cudaFree(tree_center_of_mass_x_device);
    cudaFree(tree_center_of_mass_y_device);
    cudaFree(tree_center_of_mass_z_device);

    return SUCCESS;

err_memcpy_gpu_to_cpu:
err_memcpy_cpu_to_gpu:
err_gpu_memory:
    cudaFree(a_device);
    cudaFree(x_device);
    cudaFree(m_device);
    cudaFree(leaf_morton_indices_deepest_level_device);
    cudaFree(sorted_indices_device);
    cudaFree(tree_start_particle_sorted_idx_device);
    cudaFree(tree_num_particles_device);
    cudaFree(tree_num_internal_children_device);
    cudaFree(tree_idx_first_internal_child_device);
    cudaFree(tree_total_mass_device);
    cudaFree(tree_center_of_mass_x_device);
    cudaFree(tree_center_of_mass_y_device);
    cudaFree(tree_center_of_mass_z_device);
    return return_code;
}

extern "C"
{
    WIN32DLL_API int acceleration_barnes_hut_cuda(
        real *__restrict a,
        const System *__restrict system,
        AccelerationParam *__restrict acceleration_param
    )
    {
        int return_code;
    
        const int objects_count = system->objects_count;
        const real *__restrict x = system->x;
        const real *__restrict m = system->m;
        const real G = system->G;
        const real softening_length = acceleration_param->softening_length; 
        const real opening_angle = acceleration_param->opening_angle;
    
        /* Empty the input array */
        for (int i = 0; i < objects_count; i++)
        {
            a[i * 3 + 0] = 0.0;
            a[i * 3 + 1] = 0.0;
            a[i * 3 + 2] = 0.0;
        }
    
        /* Construct the octree */
        // Allocate memory
        real width;
        int64 *leaf_morton_indices_deepest_level;
        int *sorted_indices;
        int allocated_internal_nodes;
        int actual_num_internal_nodes;
        int *tree_start_particle_sorted_idx;
        int *tree_num_particles;
        int *tree_num_internal_children;
        int *tree_idx_first_internal_child;
        real *tree_total_mass;
        real *tree_center_of_mass_x;
        real *tree_center_of_mass_y;
        real *tree_center_of_mass_z;
    
        return_code = barnes_hut_setup_octree(
            &width,
            &allocated_internal_nodes,
            &actual_num_internal_nodes,
            objects_count,
            x,
            m,
            &leaf_morton_indices_deepest_level,
            &sorted_indices,
            &tree_start_particle_sorted_idx,
            &tree_num_particles,
            &tree_num_internal_children,
            &tree_idx_first_internal_child,
            &tree_total_mass,
            &tree_center_of_mass_x,
            &tree_center_of_mass_y,
            &tree_center_of_mass_z
        );
        if (return_code != SUCCESS)
        {
            goto err_octree;
        }

        /* Compute the acceleration */ 
        return_code = _compute_acceleration(
            a,
            objects_count,
            x,
            m,
            G,
            softening_length,
            opening_angle,
            width,
            actual_num_internal_nodes,
            leaf_morton_indices_deepest_level,
            sorted_indices,
            tree_start_particle_sorted_idx,
            tree_num_particles,
            tree_num_internal_children,
            tree_idx_first_internal_child,
            tree_total_mass,
            tree_center_of_mass_x,
            tree_center_of_mass_y,
            tree_center_of_mass_z
        );
        if (return_code != SUCCESS)
        {
            goto err_acceleration;
        }

        /* Free the memory */
        free(tree_start_particle_sorted_idx);
        free(tree_num_particles);
        free(tree_num_internal_children);
        free(tree_idx_first_internal_child);
        free(tree_total_mass);
        free(tree_center_of_mass_x);
        free(tree_center_of_mass_y);
        free(tree_center_of_mass_z);
        free(leaf_morton_indices_deepest_level);
        free(sorted_indices);
        return SUCCESS;

    err_acceleration:
    err_octree:
        free(tree_start_particle_sorted_idx);
        free(tree_num_particles);
        free(tree_num_internal_children);
        free(tree_idx_first_internal_child);
        free(tree_total_mass);
        free(tree_center_of_mass_x);
        free(tree_center_of_mass_y);
        free(tree_center_of_mass_z);
        free(leaf_morton_indices_deepest_level);
        free(sorted_indices);
        return return_code;
    }
}

__global__ void memcpy_particles_array_double_to_float(
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

__global__ void memcpy_nodes_array_double_to_float(
    const double *__restrict tree_total_mass_double,
    const double *__restrict tree_center_of_mass_x_double,
    const double *__restrict tree_center_of_mass_y_double,
    const double *__restrict tree_center_of_mass_z_double,
    const int actual_num_internal_nodes,
    float *__restrict tree_total_mass,
    float *__restrict tree_center_of_mass_x,
    float *__restrict tree_center_of_mass_y,
    float *__restrict tree_center_of_mass_z
)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= actual_num_internal_nodes)
    {
        return;
    }

    tree_total_mass[i] = tree_total_mass_double[i];
    tree_center_of_mass_x[i] = tree_center_of_mass_x_double[i];
    tree_center_of_mass_y[i] = tree_center_of_mass_y_double[i];
    tree_center_of_mass_z[i] = tree_center_of_mass_z_double[i];

    return;
}

/**
 * \brief Helper kernel function for computing the acceleration of one particle in single precision
 * 
 * \param a Array of acceleration vectors
 * \param objects_count Number of objects
 * \param x Array of position vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param softening_length Softening length
 * \param opening_angle Opening angle
 * \param width Width of the bounding box
 * \param leaf_morton_indices_deepest_level Array of Morton indices
 * \param sorted_indices Array of sorted indices
 * \param tree_start_particle_sorted_idx Array of start index of particles in the node
 * \param tree_num_particles Array of number of particles in the node
 * \param tree_num_internal_children Array of number of internal children
 * \param tree_idx_first_internal_child Array of index of the first internal child
 * \param tree_total_mass Array of total mass of the node
 * \param tree_center_of_mass_x Array of x-coordinate of the center of mass
 * \param tree_center_of_mass_y Array of y-coordinate of the center of mass
 * \param tree_center_of_mass_z Array of z-coordinate of the center of mass
 * 
 * \retval SUCCESS if successful
 */
 __global__ void _acceleration_helper_float_kernel(
    double *__restrict a,
    const int objects_count,
    const float *__restrict x,
    const float *__restrict m,
    const float G,
    const float softening_length,
    const float opening_angle,
    const float width,
    const int64 *__restrict leaf_morton_indices_deepest_level,
    const int *__restrict sorted_indices,
    const int *__restrict tree_start_particle_sorted_idx,
    const int *__restrict tree_num_particles,
    const int *__restrict tree_num_internal_children,
    const int *__restrict tree_idx_first_internal_child,
    const float *__restrict tree_total_mass,
    const float *__restrict tree_center_of_mass_x,
    const float *__restrict tree_center_of_mass_y,
    const float *__restrict tree_center_of_mass_z
)
{
    typedef struct Stack
    {
        int node;
        int processed_children;
        struct Stack *last;
    } Stack;

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= objects_count)
    {
        return;
    }

    const int idx_i = sorted_indices[i];
    const int64 morton_index_i = leaf_morton_indices_deepest_level[idx_i];
    const float3 x_i = make_float3(x[idx_i * 3 + 0], x[idx_i * 3 + 1], x[idx_i * 3 + 2]);
    const float softening_squared = softening_length * softening_length;

    Stack stack_pool[MORTON_MAX_LEVEL];
    Stack *stack = &(stack_pool[0]);
    stack->processed_children = -1;
    stack->last = NULL;
    stack->node = 0;
    float3 local_a = make_float3(0.0, 0.0, 0.0);

    int level = 1;

    /* Tree walk */
    while (true)
    {
        const int node = stack->node;
        for (int j = (stack->processed_children) + 1; j < tree_num_internal_children[node]; j++)
        {
            const int child_j = tree_idx_first_internal_child[node] + j;
            const int num_children_j = tree_num_internal_children[child_j];
            const int start_idx_j = tree_start_particle_sorted_idx[child_j];

            /* Leaf node */
            if (num_children_j <= 0)
            {
                const int num_particles_j = tree_num_particles[child_j];
                for (int k = 0; k < num_particles_j; k++)
                {
                    const int idx_j = sorted_indices[start_idx_j + k];
                    if (idx_i == idx_j)
                    {
                        continue;
                    }

                    // Calculate \vec{R} and its norm
                    const float3 dr = make_float3(
                        x_i.x - x[idx_j * 3 + 0],
                        x_i.y - x[idx_j * 3 + 1],
                        x_i.z - x[idx_j * 3 + 2]
                    );
                    const float R_norm = sqrt(
                        dr.x * dr.x + dr.y * dr.y + dr.z * dr.z + softening_squared
                    );

                    // Calculate the acceleration
                    const float temp_value = G * m[idx_j] / (R_norm * R_norm * R_norm);
                    local_a.x -= temp_value * dr.x;
                    local_a.y -= temp_value * dr.y;
                    local_a.z -= temp_value * dr.z;
                }

                stack->processed_children = j;
                continue;
            }

            /* Internal node */
            else
            {
                bool criteria_met = false;

                // If object i is included, then we need to traverse deeper
                const bool is_included = _check_if_included(
                    morton_index_i,
                    leaf_morton_indices_deepest_level[sorted_indices[start_idx_j]],
                    level
                );

                // Check Barnes-Hut criteria
                float3 dr;
                float norm_square;
                if (!is_included)
                {
                    dr = make_float3(
                        x_i.x - tree_center_of_mass_x[child_j],
                        x_i.y - tree_center_of_mass_y[child_j],
                        x_i.z - tree_center_of_mass_z[child_j]
                    );
                    const float width_j = width / (2 << level);
                    norm_square = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                    if (width_j / sqrt(norm_square) < opening_angle)
                    {
                        criteria_met = true;
                    }
                }

                // Traverse deeper
                if (!criteria_met)
                {
                    Stack *new_item = &(stack_pool[level + 1]);
                    new_item->node = child_j;
                    new_item->last = stack;
                    new_item->processed_children = -1;

                    stack = new_item;
                    level++;
                    break;
                }

                else
                {
                    const float R_norm = sqrt(
                        norm_square + softening_squared
                    );

                    const float temp_value = G / (R_norm * R_norm * R_norm);
                    local_a.x -= temp_value * dr.x * tree_total_mass[child_j];
                    local_a.y -= temp_value * dr.y * tree_total_mass[child_j];
                    local_a.z -= temp_value * dr.z * tree_total_mass[child_j];

                    stack->processed_children = j;
                    continue;
                }
            }
        }

        if ((stack->processed_children + 1) >= tree_num_internal_children[stack->node])
        {
            Stack *parent = stack->last;
            if (!parent)
            {
                break;
            }            
            stack = parent;
            stack->processed_children += 1;
            level--;
        }
    }

    a[idx_i * 3 + 0] = local_a.x;
    a[idx_i * 3 + 1] = local_a.y;
    a[idx_i * 3 + 2] = local_a.z;

    return;
}

/**
 * \brief Compute the acceleration of the particles with single precision
 * 
 * \param a Array of acceleration vectors
 * \param objects_count Number of objects
 * \param x Array of position vectors
 * \param m Array of masses
 * \param G Gravitational constant
 * \param softening_length Softening length
 * \param opening_angle Opening angle
 * \param width Width of the bounding box
 * \param actual_num_internal_nodes Number of internal nodes
 * \param leaf_morton_indices_deepest_level Array of Morton indices
 * \param sorted_indices Array of sorted indices
 * \param tree_start_particle_sorted_idx Array of start index of particles in the node
 * \param tree_num_particles Array of number of particles in the node
 * \param tree_num_internal_children Array of number of internal children
 * \param tree_idx_first_internal_child Array of index of the first internal child
 * \param tree_total_mass Array of total mass of the node
 * \param tree_center_of_mass_x Array of x-coordinate of the center of mass
 * \param tree_center_of_mass_y Array of y-coordinate of the center of mass
 * \param tree_center_of_mass_z Array of z-coordinate of the center of mass
 * 
 * \retval SUCCESS if successful
 */
 IN_FILE int _compute_acceleration_float(
    double *__restrict a,
    const int objects_count,
    const double *__restrict x,
    const double *__restrict m,
    const double G,
    const double softening_length,
    const double opening_angle,
    const double width,
    const int actual_num_internal_nodes,
    const int64 *__restrict leaf_morton_indices_deepest_level,
    const int *__restrict sorted_indices,
    const int *__restrict tree_start_particle_sorted_idx,
    const int *__restrict tree_num_particles,
    const int *__restrict tree_num_internal_children,
    const int *__restrict tree_idx_first_internal_child,
    const double *__restrict tree_total_mass,
    const double *__restrict tree_center_of_mass_x,
    const double *__restrict tree_center_of_mass_y,
    const double *__restrict tree_center_of_mass_z
)
{
    int return_code;

    double *__restrict a_double_device = NULL;
    double *__restrict x_double_device = NULL;
    double *__restrict m_double_device = NULL;
    double *__restrict tree_total_mass_double_device = NULL;
    double *__restrict tree_center_of_mass_x_double_device = NULL;
    double *__restrict tree_center_of_mass_y_double_device = NULL;
    double *__restrict tree_center_of_mass_z_double_device = NULL;

    float *__restrict x_device = NULL;
    float *__restrict m_device = NULL;
    int64 *__restrict leaf_morton_indices_deepest_level_device = NULL;
    int *__restrict sorted_indices_device = NULL;
    int *__restrict tree_start_particle_sorted_idx_device = NULL;
    int *__restrict tree_num_particles_device = NULL;
    int *__restrict tree_num_internal_children_device = NULL;
    int *__restrict tree_idx_first_internal_child_device = NULL;
    float *__restrict tree_total_mass_device = NULL;
    float *__restrict tree_center_of_mass_x_device = NULL;
    float *__restrict tree_center_of_mass_y_device = NULL;
    float *__restrict tree_center_of_mass_z_device = NULL;
    cudaError_t error;

    /* Allocate memory on GPU */
    error = cudaMalloc((double **) &a_double_device, objects_count * 3 * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &x_double_device, objects_count * 3 * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &m_double_device, objects_count * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &tree_total_mass_double_device, actual_num_internal_nodes * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &tree_center_of_mass_x_double_device, actual_num_internal_nodes * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &tree_center_of_mass_y_double_device, actual_num_internal_nodes * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((double **) &tree_center_of_mass_z_double_device, actual_num_internal_nodes * sizeof(double));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }

    error = cudaMalloc((float **) &x_device, objects_count * 3 * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &m_device, objects_count * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int64 **) &leaf_morton_indices_deepest_level_device, objects_count * sizeof(int64));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &sorted_indices_device, objects_count * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &tree_start_particle_sorted_idx_device, actual_num_internal_nodes * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &tree_num_particles_device, actual_num_internal_nodes * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &tree_num_internal_children_device, actual_num_internal_nodes * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((int **) &tree_idx_first_internal_child_device, actual_num_internal_nodes * sizeof(int));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &tree_total_mass_device, actual_num_internal_nodes * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &tree_center_of_mass_x_device, actual_num_internal_nodes * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &tree_center_of_mass_y_device, actual_num_internal_nodes * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }
    error = cudaMalloc((float **) &tree_center_of_mass_z_device, actual_num_internal_nodes * sizeof(float));
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMORY_ALLOC;
        goto err_gpu_memory;
    }

    error = cudaMemcpy(x_double_device, x, objects_count * 3 * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(m_double_device, m, objects_count * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_total_mass_double_device, tree_total_mass, actual_num_internal_nodes * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_center_of_mass_x_double_device, tree_center_of_mass_x, actual_num_internal_nodes * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_center_of_mass_y_double_device, tree_center_of_mass_y, actual_num_internal_nodes * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_center_of_mass_z_double_device, tree_center_of_mass_z, actual_num_internal_nodes * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    } 

    error = cudaMemcpy(leaf_morton_indices_deepest_level_device, leaf_morton_indices_deepest_level, objects_count * sizeof(int64), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(sorted_indices_device, sorted_indices, objects_count * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_start_particle_sorted_idx_device, tree_start_particle_sorted_idx, actual_num_internal_nodes * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_num_particles_device, tree_num_particles, actual_num_internal_nodes * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_num_internal_children_device, tree_num_internal_children, actual_num_internal_nodes * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }
    error = cudaMemcpy(tree_idx_first_internal_child_device, tree_idx_first_internal_child, actual_num_internal_nodes * sizeof(int), cudaMemcpyHostToDevice);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_CPU_TO_GPU;
        goto err_memcpy_cpu_to_gpu;
    }

    memcpy_particles_array_double_to_float <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>>(
        x_double_device,
        m_double_device,
        objects_count,
        x_device,
        m_device
    );

    memcpy_nodes_array_double_to_float <<< (actual_num_internal_nodes + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>>(
        tree_total_mass_double_device,
        tree_center_of_mass_x_double_device,
        tree_center_of_mass_y_double_device,
        tree_center_of_mass_z_double_device,
        actual_num_internal_nodes,
        tree_total_mass_device,
        tree_center_of_mass_x_device,
        tree_center_of_mass_y_device,
        tree_center_of_mass_z_device
    );

    cudaFree(x_double_device);
    cudaFree(m_double_device);
    cudaFree(tree_total_mass_double_device);
    cudaFree(tree_center_of_mass_x_double_device);
    cudaFree(tree_center_of_mass_y_double_device);
    cudaFree(tree_center_of_mass_z_double_device);

    _acceleration_helper_float_kernel <<< (objects_count + BLOCK_SIZE - 1) / BLOCK_SIZE, BLOCK_SIZE >>> (
        a_double_device,
        objects_count,
        x_device,
        m_device,
        G,
        softening_length,
        opening_angle,
        width,
        leaf_morton_indices_deepest_level_device,
        sorted_indices_device,
        tree_start_particle_sorted_idx_device,
        tree_num_particles_device,
        tree_num_internal_children_device,
        tree_idx_first_internal_child_device,
        tree_total_mass_device,
        tree_center_of_mass_x_device,
        tree_center_of_mass_y_device,
        tree_center_of_mass_z_device
    );

    error = cudaMemcpy(a, a_double_device, objects_count * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    if (error != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
        return_code = ERROR_CUDA_BARNES_HUT_MEMCPY_GPU_TO_CPU;
        goto err_memcpy_gpu_to_cpu;
    }

    cudaFree(a_double_device);
    cudaFree(x_device);
    cudaFree(m_device);
    cudaFree(leaf_morton_indices_deepest_level_device);
    cudaFree(sorted_indices_device);
    cudaFree(tree_start_particle_sorted_idx_device);
    cudaFree(tree_num_particles_device);
    cudaFree(tree_num_internal_children_device);
    cudaFree(tree_idx_first_internal_child_device);
    cudaFree(tree_total_mass_device);
    cudaFree(tree_center_of_mass_x_device);
    cudaFree(tree_center_of_mass_y_device);
    cudaFree(tree_center_of_mass_z_device);

    return SUCCESS;

err_memcpy_cpu_to_gpu:
err_gpu_memory:
    cudaFree(x_double_device);
    cudaFree(m_double_device);
    cudaFree(tree_total_mass_double_device);
    cudaFree(tree_center_of_mass_x_double_device);
    cudaFree(tree_center_of_mass_y_double_device);
    cudaFree(tree_center_of_mass_z_double_device);
err_memcpy_gpu_to_cpu:
    cudaFree(a_double_device);
    cudaFree(x_device);
    cudaFree(m_device);
    cudaFree(leaf_morton_indices_deepest_level_device);
    cudaFree(sorted_indices_device);
    cudaFree(tree_start_particle_sorted_idx_device);
    cudaFree(tree_num_particles_device);
    cudaFree(tree_num_internal_children_device);
    cudaFree(tree_idx_first_internal_child_device);
    cudaFree(tree_total_mass_device);
    cudaFree(tree_center_of_mass_x_device);
    cudaFree(tree_center_of_mass_y_device);
    cudaFree(tree_center_of_mass_z_device);
    return return_code;
}


extern "C"
{
    WIN32DLL_API int acceleration_barnes_hut_cuda_float(
        real *__restrict a,
        const System *__restrict system,
        AccelerationParam *__restrict acceleration_param
    )
    {
        int return_code;
    
        const int objects_count = system->objects_count;
        const real *__restrict x = system->x;
        const real *__restrict m = system->m;
        const real G = system->G;
        const real softening_length = acceleration_param->softening_length; 
        const real opening_angle = acceleration_param->opening_angle;
    
        /* Empty the input array */
        for (int i = 0; i < objects_count; i++)
        {
            a[i * 3 + 0] = 0.0;
            a[i * 3 + 1] = 0.0;
            a[i * 3 + 2] = 0.0;
        }
    
        /* Construct the octree */
        // Allocate memory
        real width;
        int64 *leaf_morton_indices_deepest_level;
        int *sorted_indices;
        int allocated_internal_nodes;
        int actual_num_internal_nodes;
        int *tree_start_particle_sorted_idx;
        int *tree_num_particles;
        int *tree_num_internal_children;
        int *tree_idx_first_internal_child;
        real *tree_total_mass;
        real *tree_center_of_mass_x;
        real *tree_center_of_mass_y;
        real *tree_center_of_mass_z;
    
        return_code = barnes_hut_setup_octree(
            &width,
            &allocated_internal_nodes,
            &actual_num_internal_nodes,
            objects_count,
            x,
            m,
            &leaf_morton_indices_deepest_level,
            &sorted_indices,
            &tree_start_particle_sorted_idx,
            &tree_num_particles,
            &tree_num_internal_children,
            &tree_idx_first_internal_child,
            &tree_total_mass,
            &tree_center_of_mass_x,
            &tree_center_of_mass_y,
            &tree_center_of_mass_z
        );
        if (return_code != SUCCESS)
        {
            goto err_octree;
        }

        /* Compute the acceleration */ 
        return_code = _compute_acceleration_float(
            a,
            objects_count,
            x,
            m,
            G,
            softening_length,
            opening_angle,
            width,
            actual_num_internal_nodes,
            leaf_morton_indices_deepest_level,
            sorted_indices,
            tree_start_particle_sorted_idx,
            tree_num_particles,
            tree_num_internal_children,
            tree_idx_first_internal_child,
            tree_total_mass,
            tree_center_of_mass_x,
            tree_center_of_mass_y,
            tree_center_of_mass_z
        );
        if (return_code != SUCCESS)
        {
            goto err_acceleration;
        }

        /* Free the memory */
        free(tree_start_particle_sorted_idx);
        free(tree_num_particles);
        free(tree_num_internal_children);
        free(tree_idx_first_internal_child);
        free(tree_total_mass);
        free(tree_center_of_mass_x);
        free(tree_center_of_mass_y);
        free(tree_center_of_mass_z);
        free(leaf_morton_indices_deepest_level);
        free(sorted_indices);
        return SUCCESS;

    err_acceleration:
    err_octree:
        free(tree_start_particle_sorted_idx);
        free(tree_num_particles);
        free(tree_num_internal_children);
        free(tree_idx_first_internal_child);
        free(tree_total_mass);
        free(tree_center_of_mass_x);
        free(tree_center_of_mass_y);
        free(tree_center_of_mass_z);
        free(leaf_morton_indices_deepest_level);
        free(sorted_indices);
        return return_code;
    }
}
