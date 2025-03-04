#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acceleration_cuda.cuh"
#include "error.h"
#include "gravity_sim.h"
#include "math_functions.h"

#define MAX_NUM_PARTICLES_PER_LEAF 1 // Note: Potential optimization, fix this to 1 we can remove some loops and we may not need tree_num_particles array
#define MORTON_MAX_LEVEL 21 // Maximum level for 64-bit Morton index, don't change

/**
 * \brief Calculate the bounding box of the system
 * 
 * \param objects_count Number of objects
 * \param x Array of position vectors
 * \param center 3D vector of the center of the bounding box
 * \param width Width of the bounding box
 */
IN_FILE void _calculate_bounding_box(
    const int objects_count,
    const double *__restrict x,
    double *__restrict center,
    double *__restrict width
)
{
    /* Find the width of the bounding box */
    double min_x = x[0];
    double max_x = x[0];
    double min_y = x[1];
    double max_y = x[1];
    double min_z = x[2];
    double max_z = x[2];

    for (int i = 1; i < objects_count; i++)
    {
        min_x = fmin(min_x, x[i * 3 + 0]);
        max_x = fmax(max_x, x[i * 3 + 0]);
        min_y = fmin(min_y, x[i * 3 + 1]);
        max_y = fmax(max_y, x[i * 3 + 1]);
        min_z = fmin(min_z, x[i * 3 + 2]);
        max_z = fmax(max_z, x[i * 3 + 2]);
    }

    center[0] = (max_x + min_x) / 2.0;
    center[1] = (max_y + min_y) / 2.0;
    center[2] = (max_z + min_z) / 2.0;

    double width_x = max_x - min_x;
    double width_y = max_y - min_y;
    double width_z = max_z - min_z;
    *width = fmax(fmax(width_x, width_y), width_z);
}

/**
 * \brief Compute the 3D Morton indices at level 21 using magic number
 * 
 * \param morton_indices Array of Morton indices
 * \param object_count Number of objects
 * \param x Array of position vectors
 * \param center 3D vector of the center of the bounding box
 * \param width Width of the bounding box
 * 
 * \ref https://stackoverflow.com/a/18528775, Stack Overflow
 */
IN_FILE void _compute_3d_morton_indices_level_21(
    int64 *__restrict morton_indices,
    const int object_count,
    const double *__restrict x,
    const double *__restrict center,
    const double width
)
{
    for (int i = 0; i < object_count; i++)
    {
        /* Normalize the position */
        const double x_i = (x[i * 3 + 0] - center[0]) / width + 0.5;
        const double y_i = (x[i * 3 + 1] - center[1]) / width + 0.5;
        const double z_i = (x[i * 3 + 2] - center[2]) / width + 0.5;

        /* Compute the morton indices */
        int64 n_x = x_i * (1 << 21);
        int64 n_y = y_i * (1 << 21);
        int64 n_z = z_i * (1 << 21);

        n_x &= 0x1fffff;
        n_x = (n_x | n_x << 32) & 0x1f00000000ffff;
        n_x = (n_x | n_x << 16) & 0x1f0000ff0000ff;
        n_x = (n_x | n_x << 8)  & 0x100f00f00f00f00f;
        n_x = (n_x | n_x << 4)  & 0x10c30c30c30c30c3;
        n_x = (n_x | n_x << 2)  & 0x1249249249249249;
        
        n_y &= 0x1fffff;
        n_y = (n_y | n_y << 32) & 0x1f00000000ffff;
        n_y = (n_y | n_y << 16) & 0x1f0000ff0000ff;
        n_y = (n_y | n_y << 8)  & 0x100f00f00f00f00f;
        n_y = (n_y | n_y << 4)  & 0x10c30c30c30c30c3;
        n_y = (n_y | n_y << 2)  & 0x1249249249249249;

        n_z &= 0x1fffff;
        n_z = (n_z | n_z << 32) & 0x1f00000000ffff;
        n_z = (n_z | n_z << 16) & 0x1f0000ff0000ff;
        n_z = (n_z | n_z << 8)  & 0x100f00f00f00f00f;
        n_z = (n_z | n_z << 4)  & 0x10c30c30c30c30c3;
        n_z = (n_z | n_z << 2)  & 0x1249249249249249;

        morton_indices[i] = n_x | (n_y << 1) | (n_z << 2);
    }
}

/**
 * \brief Perform radix sort on the particles based on their Morton indices
 * 
 * \param object_count Number of objects
 * \param morton_indices Array of Morton indices
 * \param indices Array of indices
 * \param level Level of the Morton indices
 * 
 * \retval SUCCESS if successful
 * \retval ERROR_BARNES_HUT_RADIX_SORT_MEMORY_ALLOC if memory allocation fails
 */
IN_FILE int _radix_sort_particles_morton_index(
    const int object_count,
    int64 *__restrict morton_indices,
    int *__restrict indices,
    const int level
)
{
    int return_code;

    /* Calculate constnats */
    const int RADIX_BITS = 9;
    const int RADIX_SIZE = 1 << RADIX_BITS;
    const int RADIX_MASK = RADIX_SIZE - 1;
    
    const int num_significant_bits = 3 * level;
    const int num_passes = (num_significant_bits + RADIX_BITS - 1) / RADIX_BITS;

    // Flag to indicate whether the sorted array is in temp arrays
    // This can reduce the number of memcpy to O(1) instead of O(num_passes)
    bool is_temp = false; 

    /* Allocate memory */
    int64 *__restrict temp_morton_indices = (int64*) malloc(object_count * sizeof(int64));
    int *__restrict temp_indices = (int*) malloc(object_count * sizeof(int));
    int *__restrict count = (int*) malloc(RADIX_SIZE * sizeof(int));
    if (!temp_morton_indices || !temp_indices || !count)
    {
        return_code = ERROR_BARNES_HUT_RADIX_SORT_MEMORY_ALLOC;
        goto err_memory;
    }
    
    /* Perform LSB radix sort */    
    for (int i = 0; i < num_passes; i++) 
    {
        // Empty count array
        for (int j = 0; j < RADIX_SIZE; j++)
        {
            count[j] = 0;
        }

        // Calculate shift for this pass (start from least significant bits)
        const int shift = i * RADIX_BITS;
        
        // Count occurrences of each radix value
        if (is_temp)
        {
            for (int j = 0; j < object_count; j++) 
            {
                count[(temp_morton_indices[j] >> shift) & RADIX_MASK]++;
            }
        }
        else
        {
            for (int j = 0; j < object_count; j++) 
            {
                count[(morton_indices[j] >> shift) & RADIX_MASK]++;
            }
        }

        // Get cumulative count
        int total = 0;
        for (int j = 0; j < RADIX_SIZE; j++) 
        {
            int old_count = count[j];
            count[j] = total;
            total += old_count;
        }
        
        // Sort elements into temporary arrays
        if (is_temp)
        {
            for (int j = 0; j < object_count; j++) 
            {
                const int dest = count[(temp_morton_indices[j] >> shift) & RADIX_MASK]++;
                
                morton_indices[dest] = temp_morton_indices[j];
                indices[dest] = temp_indices[j];
            }
        }
        else
        {
            for (int j = 0; j < object_count; j++) 
            {
                const int dest = count[(morton_indices[j] >> shift) & RADIX_MASK]++;
                
                temp_morton_indices[dest] = morton_indices[j];
                temp_indices[dest] = indices[j];
            }
        }
        
        is_temp = !is_temp;
    }

    // Copy the sorted array to the original array
    if (is_temp)
    {
        memcpy(morton_indices, temp_morton_indices, object_count * sizeof(int64));
        memcpy(indices, temp_indices, object_count * sizeof(int));
    }
    
    free(count);
    free(temp_morton_indices);
    free(temp_indices);

    return SUCCESS;

err_memory:
    free(count);
    free(temp_morton_indices);
    free(temp_indices);

    return return_code;
}

/**
 * \brief Perform binary search to find the number of particles in each octant
 * 
 * \param leaf_morton_indices_deepest_level Array of Morton indices
 * \param node_morton_index_level Morton index of the node
 * \param start_idx Start index of the particles in the node
 * \param end_idx End index of the particles in the node
 * \param leaf_level Level of the leaf nodes
 * \param num_particles_per_octant Array to store the number of particles in each octant
 */
IN_FILE void _binary_search_num_particles_per_octant(
    const int64 *__restrict leaf_morton_indices_deepest_level,
    const int64 node_morton_index_level,
    const int start_idx,
    const int end_idx,
    const int leaf_level,
    int *__restrict num_particles_per_octant
)
{
    const int64 prefix = node_morton_index_level * 8;
    const int level_shift = 3 * (MORTON_MAX_LEVEL - leaf_level);

    int cumulative_count = 0;

    for (int i = 0; i < 8; i++)
    {
        // Binary search for the index of last i
        int left = start_idx + cumulative_count;
        int right = end_idx;
        while (left <= right)
        {
            const int mid = left + (right - left) / 2;
            const int mid_octant = ((leaf_morton_indices_deepest_level[mid] >> level_shift) - prefix);

            if (mid_octant < 0 || mid_octant > 7)
            {
                printf("Warning: mid_octant out of range: %d\n", mid_octant);
            }

            if (mid_octant == i && (mid == end_idx || (((leaf_morton_indices_deepest_level[mid + 1] >> level_shift) - prefix)) > i))
            {
                num_particles_per_octant[i] = mid - (start_idx + cumulative_count) + 1;
                cumulative_count += num_particles_per_octant[i];
                break;
            }
            else if (mid_octant <= i)
            {
                left = mid + 1;
            }
            else
            {
                right = mid - 1;
            }
        }
    }
}

/**
 * \brief Set up a new internal node
 * 
 * \param allocated_internal_nodes Pointer to the number of allocated internal nodes
 * \param internal_node_count Pointer to the number of internal nodes
 * \param level Node level
 * \param width Node width
 * \param node Node index
 * \param node_morton_index_level Morton index of the node at the current level
 * \param leaf_morton_indices_deepest_level Pointer to array of Morton indices
 * \param tree_num_internal_children Pointer to array of number of internal children
 * \param tree_idx_first_internal_child Pointer to array of index of the first internal child
 * \param tree_start_particle_sorted_idx Pointer to array of start index of particles in the node
 * \param tree_num_particles Pointer to array of number of particles in the node
 * \param tree_total_mass Pointer to array of total mass of the node
 * \param tree_center_of_mass_x Pointer to array of x-coordinate of the center of mass
 * \param tree_center_of_mass_y Pointer to array of y-coordinate of the center of mass
 * \param tree_center_of_mass_z Pointer to array of z-coordinate of the center of mass
 * 
 * \retval SUCCESS if successful
 * \retval error_code if there is an error
 */
IN_FILE int _setup_node(
    int *__restrict allocated_internal_nodes,
    int *__restrict internal_node_count,
    const int level,
    const double width,
    const int node,
    const int64 node_morton_index_level,
    const int64 *__restrict leaf_morton_indices_deepest_level,
    int **tree_num_internal_children,
    int **tree_idx_first_internal_child,
    int **tree_start_particle_sorted_idx,
    int **tree_num_particles,
    double **tree_total_mass,
    double **tree_center_of_mass_x,
    double **tree_center_of_mass_y,
    double **tree_center_of_mass_z
)
{
    int return_code;

    int num_particles_per_octant[8] = {0};
    const int start_idx = (*tree_start_particle_sorted_idx)[node];
    const int end_idx = start_idx + (*tree_num_particles)[node] - 1;
    const int child_level = level + 1;
    _binary_search_num_particles_per_octant(
        leaf_morton_indices_deepest_level,
        node_morton_index_level,
        start_idx,
        end_idx,
        child_level,
        num_particles_per_octant
    );

    bool first_child_found = false;
    int cumulative_count = 0;
    for (int i = 0; i < 8; i++)
    {
        if (num_particles_per_octant[i] == 0)
        {
            continue;
        }

        const int child = *internal_node_count;

        // Reallocate memory if necessary
        if (child >= *allocated_internal_nodes)
        {
            *allocated_internal_nodes *= 2;
            int *tmp_tree_num_internal_children = (int*) realloc(*tree_num_internal_children, *allocated_internal_nodes * sizeof(int));
            int *tmp_tree_idx_first_internal_child = (int*) realloc(*tree_idx_first_internal_child, *allocated_internal_nodes * sizeof(int));
            int *tmp_tree_start_particle_sorted_idx = (int*) realloc(*tree_start_particle_sorted_idx, *allocated_internal_nodes * sizeof(int));
            int *tmp_tree_num_particles = (int*) realloc(*tree_num_particles, *allocated_internal_nodes * sizeof(int));
            double *tmp_tree_total_mass = (double*) realloc(*tree_total_mass, *allocated_internal_nodes * sizeof(double));
            double *tmp_tree_center_of_mass_x = (double*) realloc(*tree_center_of_mass_x, *allocated_internal_nodes * sizeof(double));
            double *tmp_tree_center_of_mass_y = (double*) realloc(*tree_center_of_mass_y, *allocated_internal_nodes * sizeof(double));
            double *tmp_tree_center_of_mass_z = (double*) realloc(*tree_center_of_mass_z, *allocated_internal_nodes * sizeof(double));

            if (
                !tmp_tree_num_internal_children ||
                !tmp_tree_idx_first_internal_child ||
                !tmp_tree_start_particle_sorted_idx ||
                !tmp_tree_num_particles ||
                !tmp_tree_total_mass ||
                !tmp_tree_center_of_mass_x ||
                !tmp_tree_center_of_mass_y ||
                !tmp_tree_center_of_mass_z
            )
            {
                return_code = ERROR_BARNES_HUT_SETUP_NODE_MEMORY_REALLOC;
                goto err_memory_realloc;
            }

            *tree_num_internal_children = tmp_tree_num_internal_children;
            *tree_idx_first_internal_child = tmp_tree_idx_first_internal_child;
            *tree_start_particle_sorted_idx = tmp_tree_start_particle_sorted_idx;
            *tree_num_particles = tmp_tree_num_particles;
            *tree_total_mass = tmp_tree_total_mass;
            *tree_center_of_mass_x = tmp_tree_center_of_mass_x;
            *tree_center_of_mass_y = tmp_tree_center_of_mass_y;
            *tree_center_of_mass_z = tmp_tree_center_of_mass_z;
        }

        if (!first_child_found)
        {
            first_child_found = true;
            (*tree_idx_first_internal_child)[node] = child;
            (*tree_num_internal_children)[node] = 0;
        }

        // Create a new internal node
        (*tree_num_internal_children)[node] += 1;
        (*tree_start_particle_sorted_idx)[child] = start_idx + cumulative_count;
        (*tree_num_particles)[child] = num_particles_per_octant[i];

        (*tree_center_of_mass_x)[child] = (*tree_center_of_mass_x)[node];
        (*tree_center_of_mass_y)[child] = (*tree_center_of_mass_y)[node];
        (*tree_center_of_mass_z)[child] = (*tree_center_of_mass_z)[node];

        const double child_half_width = width / (2 << child_level);

        if (i & 1)
        {
            (*tree_center_of_mass_x)[child] += child_half_width;
        }
        else
        {
            (*tree_center_of_mass_x)[child] -= child_half_width;
        }

        if (i & 2)
        {
            (*tree_center_of_mass_y)[child] += child_half_width;
        }
        else
        {
            (*tree_center_of_mass_y)[child] -= child_half_width;
        }

        if (i & 4)
        {
            (*tree_center_of_mass_z)[child] += child_half_width;
        }
        else
        {
            (*tree_center_of_mass_z)[child] -= child_half_width;
        }
        
        // Update counters
        (*internal_node_count) += 1;
        cumulative_count += num_particles_per_octant[i];
    }

    return SUCCESS;

err_memory_realloc:
    return return_code;
}

/**
 * \brief Construct the octree
 * 
 * \param allocated_internal_nodes Pointer to the number of allocated internal nodes
 * \param actual_num_internal_nodes Pointer to the count of internal nodes
 * \param x Array of position vectors
 * \param m Array of masses
 * \param width Width of the bounding box
 * \param sorted_indices Array of sorted indices
 * \param leaf_morton_indices_deepest_level Array of Morton indices
 * \param morton_max_level Maximum level of the Morton indices
 * \param tree_start_particle_sorted_idx Pointer to array of start index of particles in the node
 * \param tree_num_particles Pointer to array of number of particles in the node
 * \param tree_num_internal_children Pointer to array of number of internal children
 * \param tree_idx_first_internal_child Pointer to array of index of the first internal child
 * \param tree_total_mass Pointer to array of total mass of the node
 * \param tree_center_of_mass_x Pointer to array of x-coordinate of the center of mass
 * \param tree_center_of_mass_y Pointer to array of y-coordinate of the center of mass
 * \param tree_center_of_mass_z Pointer to array of z-coordinate of the center of mass
 * 
 * \retval SUCCESS if successful
 * \retval error_code if there is an error
 */
IN_FILE int _construct_octree(
    int *__restrict allocated_internal_nodes,
    int *__restrict actual_num_internal_nodes,
    const double *__restrict x,
    const double *__restrict m,
    const double width,
    const int *__restrict sorted_indices,
    const int64 *__restrict leaf_morton_indices_deepest_level,
    const int morton_max_level,
    int **tree_start_particle_sorted_idx,
    int **tree_num_particles,
    int **tree_num_internal_children,
    int **tree_idx_first_internal_child,
    double **tree_total_mass,
    double **tree_center_of_mass_x,
    double **tree_center_of_mass_y,
    double **tree_center_of_mass_z
)
{
    typedef struct Stack
    {
        int node;
        int processed_children;
        struct Stack *last;
        double total_mass;
        double mass_times_distance[3];
    } Stack;

    int return_code;

    /* Create a stack pool */
    Stack stack_pool[MORTON_MAX_LEVEL];
    Stack *stack = &(stack_pool[0]);

    stack->node = 0;
    stack->processed_children = -1;
    stack->last = NULL;
    stack->total_mass = 0.0;
    stack->mass_times_distance[0] = 0.0;
    stack->mass_times_distance[1] = 0.0;
    stack->mass_times_distance[2] = 0.0;

    int level = 0;
    int internal_node_count = 1;

    /* Set up the root node */
    return_code = _setup_node(
        allocated_internal_nodes,
        &internal_node_count,
        level,
        width,
        stack->node,
        0,
        leaf_morton_indices_deepest_level,
        tree_num_internal_children,
        tree_idx_first_internal_child,
        tree_start_particle_sorted_idx,
        tree_num_particles,
        tree_total_mass,
        tree_center_of_mass_x,
        tree_center_of_mass_y,
        tree_center_of_mass_z
    );
    if (return_code != SUCCESS)
    {
        goto err_setup_node;
    }

    level++;

    while (true)
    {
        const int node = stack->node;
        for (int i = stack->processed_children + 1; i < (*tree_num_internal_children)[node]; i++)
        {
            const int child = (*tree_idx_first_internal_child)[node] + i;
            const int start_idx = (*tree_start_particle_sorted_idx)[child];
            const int num_particles = (*tree_num_particles)[child];

            const int64 child_morton_index_level = (leaf_morton_indices_deepest_level[start_idx] >> (3 * (MORTON_MAX_LEVEL - level)));

            if (num_particles <= MAX_NUM_PARTICLES_PER_LEAF || level >= morton_max_level)
            {
                // Leaf node
                (*tree_num_internal_children)[child] = 0;

                // Update the stack
                stack->processed_children = i;

                for (int j = 0; j < num_particles; j++)
                {
                    const int particle_idx = sorted_indices[start_idx + j];
                    stack->total_mass += m[particle_idx];
                    stack->mass_times_distance[0] += m[particle_idx] * x[particle_idx * 3 + 0];
                    stack->mass_times_distance[1] += m[particle_idx] * x[particle_idx * 3 + 1];
                    stack->mass_times_distance[2] += m[particle_idx] * x[particle_idx * 3 + 2];
                }

                continue;
            }
            else
            {
                return_code = _setup_node(
                    allocated_internal_nodes,
                    &internal_node_count,
                    level,
                    width,
                    child,
                    child_morton_index_level,
                    leaf_morton_indices_deepest_level,
                    tree_num_internal_children,
                    tree_idx_first_internal_child,
                    tree_start_particle_sorted_idx,
                    tree_num_particles,
                    tree_total_mass,
                    tree_center_of_mass_x,
                    tree_center_of_mass_y,
                    tree_center_of_mass_z
                );
                if (return_code != SUCCESS)
                {
                    goto err_setup_node;
                }

                Stack *new_item = &(stack_pool[level + 1]);
                new_item->node = child;
                new_item->last = stack;
                new_item->processed_children = -1;
                new_item->total_mass = 0.0;
                new_item->mass_times_distance[0] = 0.0;
                new_item->mass_times_distance[1] = 0.0;
                new_item->mass_times_distance[2] = 0.0;

                stack = new_item;
                level++;

                break;
            }
        }

        if ((stack->processed_children + 1) >= (*tree_num_internal_children)[stack->node])
        {
            /* Update center of mass */
            (*tree_total_mass)[node] = stack->total_mass;
            (*tree_center_of_mass_x)[node] = stack->mass_times_distance[0] / stack->total_mass;
            (*tree_center_of_mass_y)[node] = stack->mass_times_distance[1] / stack->total_mass;
            (*tree_center_of_mass_z)[node] = stack->mass_times_distance[2] / stack->total_mass;

            Stack *parent = stack->last;
            if (!parent)
            {
                break;
            }

            parent->total_mass += stack->total_mass;
            parent->mass_times_distance[0] += stack->mass_times_distance[0];
            parent->mass_times_distance[1] += stack->mass_times_distance[1];
            parent->mass_times_distance[2] += stack->mass_times_distance[2];
            
            stack = parent;
            stack->processed_children += 1;
            level--;
        }
    }

    *actual_num_internal_nodes = internal_node_count;

    return SUCCESS;

err_setup_node:
    return return_code;
}

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
        double *__restrict a,
        const System *__restrict system,
        AccelerationParam *__restrict acceleration_param
    )
    {
        int return_code;

        const int objects_count = system->objects_count;
        const double *__restrict x = system->x;
        const double *__restrict m = system->m;
        const double G = system->G;
        const double softening_length = acceleration_param->softening_length; 
        const double opening_angle = acceleration_param->opening_angle;

        /* Empty the input array */
        for (int i = 0; i < objects_count; i++)
        {
            a[i * 3 + 0] = 0.0;
            a[i * 3 + 1] = 0.0;
            a[i * 3 + 2] = 0.0;
        }

        /* Find the width and center of the bounding box */
        double center[3];
        double width;
        _calculate_bounding_box(objects_count, x, center, &width);

        /* Construct the octree */
        // Variables to be used in the octree
        int factor;
        int allocated_internal_nodes;
        int actual_num_internal_nodes;
        int *tree_start_particle_sorted_idx;
        int *tree_num_particles;
        int *tree_num_internal_children;
        int *tree_idx_first_internal_child;
        double *tree_total_mass;
        double *tree_center_of_mass_x;
        double *tree_center_of_mass_y;
        double *tree_center_of_mass_z;

        // Allocate memory
        int64 *leaf_morton_indices_deepest_level = (int64*) malloc(objects_count * sizeof(int64));
        int *sorted_indices = (int*) malloc(objects_count * sizeof(int));
        if (!leaf_morton_indices_deepest_level || !sorted_indices)
        {
            return_code = ERROR_BARNES_HUT_MORTON_INDICES_MEMORY_ALLOC;
            goto err_morton_indices_memory_alloc;
        }

        for (int i = 0; i < objects_count; i++)
        {
            sorted_indices[i] = i;
        }
        _compute_3d_morton_indices_level_21(
            leaf_morton_indices_deepest_level,
            objects_count,
            x,
            center,
            width
        );
        return_code = _radix_sort_particles_morton_index(
            objects_count,
            leaf_morton_indices_deepest_level,
            sorted_indices,
            MORTON_MAX_LEVEL
        );
        if (return_code != SUCCESS)
        {
            goto err_radix_sort;
        }

        // Allocate memory for the octree
        factor = 1;
        if (MAX_NUM_PARTICLES_PER_LEAF <= 2)
        {
            factor = 2;
        }
        allocated_internal_nodes = factor * objects_count;

        // Start index of the particles in the node
        tree_start_particle_sorted_idx = (int*) malloc(allocated_internal_nodes * sizeof(int));

        // Number of particles in the node
        tree_num_particles = (int*) malloc(allocated_internal_nodes * sizeof(int));

        // Number of internal children of the node (i.e. not leaf)
        tree_num_internal_children = (int*) malloc(allocated_internal_nodes * sizeof(int));

        // Index to the first internal child of the node
        tree_idx_first_internal_child = (int*) malloc(allocated_internal_nodes * sizeof(int));

        // Total mass of the node
        tree_total_mass = (double*) malloc(allocated_internal_nodes * sizeof(double));

        // Center of mass of the node
        tree_center_of_mass_x = (double*) malloc(allocated_internal_nodes * sizeof(double));
        tree_center_of_mass_y = (double*) malloc(allocated_internal_nodes * sizeof(double));
        tree_center_of_mass_z = (double*) malloc(allocated_internal_nodes * sizeof(double));

        if (
            !tree_start_particle_sorted_idx ||
            !tree_num_particles ||
            !tree_num_internal_children ||
            !tree_idx_first_internal_child ||
            !tree_total_mass ||
            !tree_center_of_mass_x ||
            !tree_center_of_mass_y ||
            !tree_center_of_mass_z
        )
        {
            return_code = ERROR_BARNES_HUT_OCTREE_MEMORY_ALLOC;
            goto err_octree_memory_alloc;
        }

        tree_start_particle_sorted_idx[0] = 0;
        tree_num_particles[0] = objects_count;
        tree_num_internal_children[0] = 0;
        tree_center_of_mass_x[0] = center[0];
        tree_center_of_mass_y[0] = center[1];
        tree_center_of_mass_z[0] = center[2];

        return_code = _construct_octree(
            &allocated_internal_nodes,
            &actual_num_internal_nodes,
            x,
            m,
            width,
            sorted_indices,
            leaf_morton_indices_deepest_level,
            MORTON_MAX_LEVEL,
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
    err_octree_memory_alloc:
        free(tree_start_particle_sorted_idx);
        free(tree_num_particles);
        free(tree_num_internal_children);
        free(tree_idx_first_internal_child);
        free(tree_total_mass);
        free(tree_center_of_mass_x);
        free(tree_center_of_mass_y);
        free(tree_center_of_mass_z);
    err_radix_sort:
    err_morton_indices_memory_alloc:
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
        double *__restrict a,
        const System *__restrict system,
        AccelerationParam *__restrict acceleration_param
    )
    {
        int return_code;

        const int objects_count = system->objects_count;
        const double *__restrict x = system->x;
        const double *__restrict m = system->m;
        const double G = system->G;
        const double softening_length = acceleration_param->softening_length; 
        const double opening_angle = acceleration_param->opening_angle;

        /* Empty the input array */
        for (int i = 0; i < objects_count; i++)
        {
            a[i * 3 + 0] = 0.0;
            a[i * 3 + 1] = 0.0;
            a[i * 3 + 2] = 0.0;
        }

        /* Find the width and center of the bounding box */
        double center[3];
        double width;
        _calculate_bounding_box(objects_count, x, center, &width);

        /* Construct the octree */
        // Variables to be used in the octree
        int factor;
        int allocated_internal_nodes;
        int actual_num_internal_nodes;
        int *tree_start_particle_sorted_idx;
        int *tree_num_particles;
        int *tree_num_internal_children;
        int *tree_idx_first_internal_child;
        double *tree_total_mass;
        double *tree_center_of_mass_x;
        double *tree_center_of_mass_y;
        double *tree_center_of_mass_z;

        // Allocate memory
        int64 *leaf_morton_indices_deepest_level = (int64*) malloc(objects_count * sizeof(int64));
        int *sorted_indices = (int*) malloc(objects_count * sizeof(int));
        if (!leaf_morton_indices_deepest_level || !sorted_indices)
        {
            return_code = ERROR_BARNES_HUT_MORTON_INDICES_MEMORY_ALLOC;
            goto err_morton_indices_memory_alloc;
        }

        for (int i = 0; i < objects_count; i++)
        {
            sorted_indices[i] = i;
        }
        _compute_3d_morton_indices_level_21(
            leaf_morton_indices_deepest_level,
            objects_count,
            x,
            center,
            width
        );
        return_code = _radix_sort_particles_morton_index(
            objects_count,
            leaf_morton_indices_deepest_level,
            sorted_indices,
            MORTON_MAX_LEVEL
        );
        if (return_code != SUCCESS)
        {
            goto err_radix_sort;
        }

        // Allocate memory for the octree
        factor = 1;
        if (MAX_NUM_PARTICLES_PER_LEAF <= 2)
        {
            factor = 2;
        }
        allocated_internal_nodes = factor * objects_count;

        // Start index of the particles in the node
        tree_start_particle_sorted_idx = (int*) malloc(allocated_internal_nodes * sizeof(int));

        // Number of particles in the node
        tree_num_particles = (int*) malloc(allocated_internal_nodes * sizeof(int));

        // Number of internal children of the node (i.e. not leaf)
        tree_num_internal_children = (int*) malloc(allocated_internal_nodes * sizeof(int));

        // Index to the first internal child of the node
        tree_idx_first_internal_child = (int*) malloc(allocated_internal_nodes * sizeof(int));

        // Total mass of the node
        tree_total_mass = (double*) malloc(allocated_internal_nodes * sizeof(double));

        // Center of mass of the node
        tree_center_of_mass_x = (double*) malloc(allocated_internal_nodes * sizeof(double));
        tree_center_of_mass_y = (double*) malloc(allocated_internal_nodes * sizeof(double));
        tree_center_of_mass_z = (double*) malloc(allocated_internal_nodes * sizeof(double));

        if (
            !tree_start_particle_sorted_idx ||
            !tree_num_particles ||
            !tree_num_internal_children ||
            !tree_idx_first_internal_child ||
            !tree_total_mass ||
            !tree_center_of_mass_x ||
            !tree_center_of_mass_y ||
            !tree_center_of_mass_z
        )
        {
            return_code = ERROR_BARNES_HUT_OCTREE_MEMORY_ALLOC;
            goto err_octree_memory_alloc;
        }

        tree_start_particle_sorted_idx[0] = 0;
        tree_num_particles[0] = objects_count;
        tree_num_internal_children[0] = 0;
        tree_center_of_mass_x[0] = center[0];
        tree_center_of_mass_y[0] = center[1];
        tree_center_of_mass_z[0] = center[2];

        return_code = _construct_octree(
            &allocated_internal_nodes,
            &actual_num_internal_nodes,
            x,
            m,
            width,
            sorted_indices,
            leaf_morton_indices_deepest_level,
            MORTON_MAX_LEVEL,
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
    err_octree_memory_alloc:
        free(tree_start_particle_sorted_idx);
        free(tree_num_particles);
        free(tree_num_internal_children);
        free(tree_idx_first_internal_child);
        free(tree_total_mass);
        free(tree_center_of_mass_x);
        free(tree_center_of_mass_y);
        free(tree_center_of_mass_z);
    err_radix_sort:
    err_morton_indices_memory_alloc:
        free(leaf_morton_indices_deepest_level);
        free(sorted_indices);
        return return_code;
    }
}
