#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "gravity_sim.h"
#include "math_functions.h"

#define MAX_NUM_PARTICLES_PER_LEAF 1
#define MORTON_MAX_LEVEL 21 // Maximum level for 64-bit Morton index, don't change


// // For debug only
// void print_octree_node(
//     real *restrict x,
//     real *restrict m,
//     const int *restrict sorted_indices,
//     int node_idx,
//     int *restrict tree_start_particle_sorted_idx,
//     int *restrict tree_num_particles,
//     int *restrict tree_num_internal_children,
//     int *restrict tree_idx_first_internal_child,
//     real *restrict tree_total_mass,
//     real *restrict tree_center_of_mass_x,
//     real *restrict tree_center_of_mass_y,
//     real *restrict tree_center_of_mass_z,
//     int indent
// )
// {
//     // Print indent spaces
//     for (int i = 0; i < indent; ++i)
//     printf("  ");

//     // Print summary info about the node
//     printf("Node %d:\n", node_idx);

//     for (int i = 0; i < indent; ++i) printf("  ");
//     printf("  Num Particles: %d, Num children: %d\n",
//         tree_num_particles[node_idx],
//         tree_num_internal_children[node_idx]
//     );

//     if (tree_num_internal_children[node_idx] > 0)
//     {
//         for (int i = 0; i < indent; ++i) printf("  ");
//         printf("  Center of Mass: (%.4g, %.4g, %.4g), Total Mass: %.7g\n",
//             tree_center_of_mass_x[node_idx],
//             tree_center_of_mass_y[node_idx],
//             tree_center_of_mass_z[node_idx],
//             tree_total_mass[node_idx]
//         );
//     }
//     else
//     {
//         for (int i = 0; i < tree_num_particles[node_idx]; i++)
//         {
//             int particle_idx = sorted_indices[tree_start_particle_sorted_idx[node_idx] + i];
//             for (int j = 0; j < indent; ++j) printf("  ");
//             printf("  Particle %d: (%.4g, %.4g, %.4g), m = %.4g\n",
//                 particle_idx,
//                 x[particle_idx * 3 + 0],
//                 x[particle_idx * 3 + 1],
//                 x[particle_idx * 3 + 2],
//                 m[particle_idx]
//             );
//         }
//     }

//     // Recurse on internal children (if any)
//     int num_children = tree_num_internal_children[node_idx];
//     if (num_children > 0) 
//     {
//         int first_child = tree_idx_first_internal_child[node_idx];
//         for (int i = 0; i < num_children; ++i)
//         {
//             int child_idx = first_child + i;
//             print_octree_node(
//                 x,
//                 m,
//                 sorted_indices,
//                 child_idx,
//                 tree_start_particle_sorted_idx,
//                 tree_num_particles,
//                 tree_num_internal_children,
//                 tree_idx_first_internal_child,
//                 tree_total_mass,
//                 tree_center_of_mass_x,
//                 tree_center_of_mass_y,
//                 tree_center_of_mass_z,
//                 indent + 1
//             );
//         }
//     }
// }

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
    const real *restrict x,
    real *restrict center,
    real *restrict width
)
{
    /* Find the width of the bounding box */
    real min_x = x[0];
    real max_x = x[0];
    real min_y = x[1];
    real max_y = x[1];
    real min_z = x[2];
    real max_z = x[2];

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

    real width_x = max_x - min_x;
    real width_y = max_y - min_y;
    real width_z = max_z - min_z;
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
    int64 *restrict morton_indices,
    const int object_count,
    const real *restrict x,
    const real *restrict center,
    const real width
)
{
    for (int i = 0; i < object_count; i++)
    {
        /* Normalize the position */
        const real x_i = (x[i * 3 + 0] - center[0]) / width + 0.5;
        const real y_i = (x[i * 3 + 1] - center[1]) / width + 0.5;
        const real z_i = (x[i * 3 + 2] - center[2]) / width + 0.5;

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
    int64 *restrict morton_indices,
    int *restrict indices,
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

    /* Allocate memory */
    int64 *restrict temp_morton_indices = malloc(object_count * sizeof(int64));
    int *restrict temp_indices = malloc(object_count * sizeof(int));
    int *restrict count = malloc(RADIX_SIZE * sizeof(int));
    if (!temp_morton_indices || !temp_indices || !count)
    {
        return_code = ERROR_BARNES_HUT_RADIX_SORT_MEMORY_ALLOC;
        goto err_memory;
    }
    
    /* Perform LSB radix sort */

    // Flag to indicate whether the sorted array is in temp arrays
    // This can reduce the number of memcpy to O(1) instead of O(num_passes)
    bool is_temp = false; 
    
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
    const int64 *restrict leaf_morton_indices_deepest_level,
    const int64 node_morton_index_level,
    const int start_idx,
    const int end_idx,
    const int leaf_level,
    int *restrict num_particles_per_octant
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
    int *restrict allocated_internal_nodes,
    int *restrict internal_node_count,
    const int level,
    const real width,
    const int node,
    const int64 node_morton_index_level,
    const int64 *restrict leaf_morton_indices_deepest_level,
    int **tree_num_internal_children,
    int **tree_idx_first_internal_child,
    int **tree_start_particle_sorted_idx,
    int **tree_num_particles,
    real **tree_total_mass,
    real **tree_center_of_mass_x,
    real **tree_center_of_mass_y,
    real **tree_center_of_mass_z
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
            int *tmp_tree_num_internal_children = realloc(*tree_num_internal_children, *allocated_internal_nodes * sizeof(int));
            int *tmp_tree_idx_first_internal_child = realloc(*tree_idx_first_internal_child, *allocated_internal_nodes * sizeof(int));
            int *tmp_tree_start_particle_sorted_idx = realloc(*tree_start_particle_sorted_idx, *allocated_internal_nodes * sizeof(int));
            int *tmp_tree_num_particles = realloc(*tree_num_particles, *allocated_internal_nodes * sizeof(int));
            real *tmp_tree_total_mass = realloc(*tree_total_mass, *allocated_internal_nodes * sizeof(real));
            real *tmp_tree_center_of_mass_x = realloc(*tree_center_of_mass_x, *allocated_internal_nodes * sizeof(real));
            real *tmp_tree_center_of_mass_y = realloc(*tree_center_of_mass_y, *allocated_internal_nodes * sizeof(real));
            real *tmp_tree_center_of_mass_z = realloc(*tree_center_of_mass_z, *allocated_internal_nodes * sizeof(real));

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

        const real child_half_width = width / (2 << child_level);

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
    int *restrict allocated_internal_nodes,
    const real *restrict x,
    const real *restrict m,
    const real width,
    const int *restrict sorted_indices,
    const int64 *restrict leaf_morton_indices_deepest_level,
    const int morton_max_level,
    int **tree_start_particle_sorted_idx,
    int **tree_num_particles,
    int **tree_num_internal_children,
    int **tree_idx_first_internal_child,
    real **tree_total_mass,
    real **tree_center_of_mass_x,
    real **tree_center_of_mass_y,
    real **tree_center_of_mass_z
)
{
    typedef struct Stack
    {
        int node;
        int processed_children;
        struct Stack *last;
        real total_mass;
        real mass_times_distance[3];
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
IN_FILE bool _check_if_included(
    const int64 morton_index_i,
    const int64 morton_index_j,
    const int level
)
{
    return (morton_index_i >> (3 * (MORTON_MAX_LEVEL - level))) == (morton_index_j >> (3 * (MORTON_MAX_LEVEL - level)));
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
    real *restrict a,
    const int objects_count,
    const real *restrict x,
    const real *restrict m,
    const real G,
    const real softening_length,
    const real opening_angle,
    const real width,
    const int64 *restrict leaf_morton_indices_deepest_level,
    const int *restrict sorted_indices,
    const int *restrict tree_start_particle_sorted_idx,
    const int *restrict tree_num_particles,
    const int *restrict tree_num_internal_children,
    const int *restrict tree_idx_first_internal_child,
    const real *restrict tree_total_mass,
    const real *restrict tree_center_of_mass_x,
    const real *restrict tree_center_of_mass_y,
    const real *restrict tree_center_of_mass_z
)
{
    typedef struct Stack
    {
        real acceleration[3];
        int node;
        int processed_children;
        struct Stack *last;
    } Stack;

    for (int i = 0; i < objects_count; i++)
    {
        const int idx_i = sorted_indices[i];    // Actually not necessary, we can use i directly
        const int64 morton_index_i = leaf_morton_indices_deepest_level[idx_i];
        const real x_i[3] = {x[idx_i * 3 + 0], x[idx_i * 3 + 1], x[idx_i * 3 + 2]};
        
        Stack stack_pool[MORTON_MAX_LEVEL];
        Stack *stack = &(stack_pool[0]);
        stack->processed_children = -1;
        stack->last = NULL;
        stack->node = 0;
        stack->acceleration[0] = 0.0;
        stack->acceleration[1] = 0.0;
        stack->acceleration[2] = 0.0;

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

                        const real m_j = m[idx_j];
                        real R[3];

                        // Calculate \vec{R} and its norm
                        R[0] = x_i[0] - x[idx_j * 3 + 0];
                        R[1] = x_i[1] - x[idx_j * 3 + 1];
                        R[2] = x_i[2] - x[idx_j * 3 + 2];
                        const real R_norm = sqrt(
                            R[0] * R[0] + 
                            R[1] * R[1] + 
                            R[2] * R[2] +
                            softening_length * softening_length
                        );

                        // Calculate the acceleration
                        const real temp_value = G / (R_norm * R_norm * R_norm);
                        stack->acceleration[0] -= temp_value * R[0] * m_j;
                        stack->acceleration[1] -= temp_value * R[1] * m_j;
                        stack->acceleration[2] -= temp_value * R[2] * m_j;
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
                    real R[3];
                    if (!is_included)
                    {
                        R[0] = x_i[0] - tree_center_of_mass_x[child_j];
                        R[1] = x_i[1] - tree_center_of_mass_y[child_j];
                        R[2] = x_i[2] - tree_center_of_mass_z[child_j];
                        const real width_j = width / (2 << level);

                        if (width_j / sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]) < opening_angle)
                        {
                            criteria_met = true;
                        }
                    }

                    // Traverse deeper
                    if (!criteria_met)
                    {
                        Stack *new_item = &(stack_pool[level + 1]);
                        new_item->acceleration[0] = 0.0;
                        new_item->acceleration[1] = 0.0;
                        new_item->acceleration[2] = 0.0;
                        new_item->node = child_j;
                        new_item->last = stack;
                        new_item->processed_children = -1;

                        stack = new_item;
                        level++;
                        break;
                    }

                    else
                    {
                        const real R_norm = sqrt(
                            R[0] * R[0] + 
                            R[1] * R[1] + 
                            R[2] * R[2] +
                            softening_length * softening_length
                        );

                        const real temp_value = G / (R_norm * R_norm * R_norm);
                        stack->acceleration[0] -= temp_value * R[0] * tree_total_mass[child_j];
                        stack->acceleration[1] -= temp_value * R[1] * tree_total_mass[child_j];
                        stack->acceleration[2] -= temp_value * R[2] * tree_total_mass[child_j];

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

                parent->acceleration[0] += stack->acceleration[0];
                parent->acceleration[1] += stack->acceleration[1];
                parent->acceleration[2] += stack->acceleration[2];
                
                stack = parent;
                stack->processed_children += 1;
                level--;
            }
        }

        a[idx_i * 3 + 0] = stack->acceleration[0];
        a[idx_i * 3 + 1] = stack->acceleration[1];
        a[idx_i * 3 + 2] = stack->acceleration[2];
    }

    return SUCCESS;
}

WIN32DLL_API int acceleration_barnes_hut(
    real *restrict a,
    const System *restrict system,
    AccelerationParam *restrict acceleration_param
)
{
    int return_code;

    const int objects_count = system->objects_count;
    const real *restrict x = system->x;
    const real *restrict m = system->m;
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

    /* Find the width and center of the bounding box */
    real center[3];
    real width;
    _calculate_bounding_box(objects_count, x, center, &width);

    /* Construct the octree */
    // Allocate memory
    int64 *leaf_morton_indices_deepest_level = malloc(objects_count * sizeof(int64));
    int *sorted_indices = malloc(objects_count * sizeof(int));
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
    int factor = 1;
    if (MAX_NUM_PARTICLES_PER_LEAF <= 2)
    {
        factor = 2;
    }
    int allocated_internal_nodes = factor * objects_count;

    // Start index of the particles in the node
    int *tree_start_particle_sorted_idx = malloc(allocated_internal_nodes * sizeof(int));

    // Number of particles in the node
    int *tree_num_particles = malloc(allocated_internal_nodes * sizeof(int));

    // Number of internal children of the node (i.e. not leaf)
    int *tree_num_internal_children = malloc(allocated_internal_nodes * sizeof(int));

    // Index to the first internal child of the node
    int *tree_idx_first_internal_child = malloc(allocated_internal_nodes * sizeof(int));

    // Total mass of the node
    real *tree_total_mass = malloc(allocated_internal_nodes * sizeof(real));

    // Center of mass of the node
    real *tree_center_of_mass_x = malloc(allocated_internal_nodes * sizeof(real));
    real *tree_center_of_mass_y = malloc(allocated_internal_nodes * sizeof(real));
    real *tree_center_of_mass_z = malloc(allocated_internal_nodes * sizeof(real));

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

    return_code = _compute_acceleration(
        a,
        objects_count,
        x,
        m,
        G,
        softening_length,
        opening_angle,
        width,
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

    // For debug only
    // print_octree_node(
    //     x,
    //     m,
    //     sorted_indices,
    //     0,
    //     tree_start_particle_sorted_idx,
    //     tree_num_particles,
    //     tree_num_internal_children,
    //     tree_idx_first_internal_child,
    //     tree_total_mass,
    //     tree_center_of_mass_x,
    //     tree_center_of_mass_y,
    //     tree_center_of_mass_z,
    //     0
    // );

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
