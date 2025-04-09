/**
 * \file linear_octree.c
 * \brief Implementation of linear octree for Barnes-Hut algorithm
 * 
 * \author Ching-Yin Ng
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_OPENMP
    #include <omp.h>
#endif

#include "acceleration.h"
#include "common.h"
#include "error.h"
#include "linear_octree.h"


// // For debug only
// IN_FILE void print_octree_nodes(
//     const LinearOctree *__restrict octree,
//     const double *__restrict x,
//     const double *__restrict m,
//     const int node_idx,
//     const int indent
// )
// {
//     // Print indent spaces
//     for (int i = 0; i < indent; ++i)
//     printf("  ");

//     // Print summary info about the node
//     printf("Node %d:\n", node_idx);

//     for (int i = 0; i < indent; ++i) printf("  ");
//     printf("  Num Particles: %d, Num children: %d\n",
//         octree->tree_num_particles[node_idx],
//         octree->tree_num_internal_children[node_idx]
//     );

//     if (octree->tree_num_internal_children[node_idx] > 0)
//     {
//         for (int i = 0; i < indent; ++i) printf("  ");
//         printf("  Center of Mass: (%.16g, %.16g, %.16g), Total Mass: %.16g\n",
//             octree->tree_center_of_mass_x[node_idx],
//             octree->tree_center_of_mass_y[node_idx],
//             octree->tree_center_of_mass_z[node_idx],
//             octree->tree_mass[node_idx]
//         );
//     }
//     else
//     {
//         for (int i = 0; i < octree->tree_num_particles[node_idx]; i++)
//         {
//             int particle_idx = octree->sorted_indices[octree->tree_first_particle_sorted_idx[node_idx] + i];
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
//     int num_children = octree->tree_num_internal_children[node_idx];
//     if (num_children > 0) 
//     {
//         int first_child = octree->tree_first_internal_children_idx[node_idx];
//         for (int i = 0; i < num_children; ++i)
//         {
//             int child_idx = first_child + i;
//             print_octree_nodes(
//                 octree,
//                 x,
//                 m,
//                 child_idx,
//                 indent + 1
//             );
//         }
//     }
// }

LinearOctree get_new_linear_octree(void)
{
    LinearOctree linear_octree;
    linear_octree.particle_morton_indices_deepest_level = NULL;
    linear_octree.sorted_indices = NULL;
    linear_octree.tree_num_particles = NULL;
    linear_octree.tree_num_internal_children = NULL;
    linear_octree.tree_first_internal_children_idx = NULL;
    linear_octree.tree_mass = NULL;
    linear_octree.tree_center_of_mass_x = NULL;
    linear_octree.tree_center_of_mass_y = NULL;
    linear_octree.tree_center_of_mass_z = NULL;
    return linear_octree;
}

/**
 * \brief Calculate the bounding box of the system
 * 
 * \param[out] center 3D vector of the center of the bounding box
 * \param[out] width Width of the bounding box
 * \param[in] num_particles Number of particles
 * \param[in] x Array of position vectors
 */
IN_FILE void calculate_bounding_box(
    double *__restrict center,
    double *__restrict width,
    const int num_particles,
    const double *__restrict x
)
{
    /* Find the width of the bounding box */
    double min_x = x[0];
    double max_x = x[0];
    double min_y = x[1];
    double max_y = x[1];
    double min_z = x[2];
    double max_z = x[2];

    for (int i = 1; i < num_particles; i++)
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

    const double width_x = max_x - min_x;
    const double width_y = max_y - min_y;
    const double width_z = max_z - min_z;
    *width = fmax(fmax(width_x, width_y), width_z);
}

/**
 * \brief Compute the 3D Morton indices at level 21 using magic number
 * 
 * \param[out] morton_indices Array of Morton indices
 * \param[in] object_count Number of particles
 * \param[in] x Array of position vectors
 * \param[in] center 3D vector of the center of the bounding box
 * \param[in] width Width of the bounding box
 * 
 * \ref https://stackoverflow.com/a/18528775, Stack Overflow
 */
IN_FILE void compute_3d_particle_morton_indices_deepest_level(
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
 * \param morton_indices Array of Morton indices
 * \param indices Array of indices
 * \param object_count Number of particles
 * \param level Level of the Morton indices
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_MEMORY_ERROR if memory allocation for temporary arrays failed
 */
IN_FILE ErrorStatus radix_sort_particles_morton_index(
    int64 *__restrict morton_indices,
    int *__restrict indices,
    const int object_count,
    const int level
)
{
    /* Calculate constnats */
    const int RADIX_BITS = 9;
    const int RADIX_SIZE = 1 << RADIX_BITS;
    const int RADIX_MASK = RADIX_SIZE - 1;
    
    const int num_significant_bits = 3 * level;
    const int num_passes = (num_significant_bits + RADIX_BITS - 1) / RADIX_BITS;

    /* Allocate memory */
    int64 *__restrict temp_morton_indices = malloc(object_count * sizeof(int64));
    int *__restrict temp_indices = malloc(object_count * sizeof(int));
    int *__restrict count = malloc(RADIX_SIZE * sizeof(int));
    if (!temp_morton_indices || !temp_indices || !count)
    {
        free(count);
        free(temp_morton_indices);
        free(temp_indices);

        return WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for temporary arrays"
        );
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

    return make_success_error_status();
}

/**
 * \brief Perform binary search to find the number of particles in each octant
 * 
 * \param[out] num_particles_per_octant Array to store the number of particles in each octant
 * \param[in] particle_morton_indices_deepest_level Array of Morton indices at the deepest level
 * \param[in] node_morton_index_level Morton index of the node
 * \param[in] start_idx Start index of the particles in the node
 * \param[in] end_idx End index of the particles in the node
 * \param[in] leaf_level Level of the leaf nodes
 * 
 * \return ErrorStatus
 * 
 * \exception GRAV_VALUE_ERROR if the Morton index is out of range
 */
IN_FILE ErrorStatus binary_search_num_particles_per_octant(
    int *__restrict num_particles_per_octant,
    const int64 *__restrict particle_morton_indices_deepest_level,
    const int64 node_morton_index_level,
    const int start_idx,
    const int end_idx,
    const int leaf_level
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
            const int mid_octant = ((particle_morton_indices_deepest_level[mid] >> level_shift) - prefix);

            if (mid_octant > 7 || mid_octant < 0)
            {
                return WRAP_RAISE_ERROR_FMT(
                    GRAV_VALUE_ERROR,
                    "Morton index %d is out of range [0, 7]",
                    mid_octant
                );
            }

            if (mid_octant == i && (mid == end_idx || (((particle_morton_indices_deepest_level[mid + 1] >> level_shift) - prefix)) > i))
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

    return make_success_error_status();
}

/**
 * \brief Set up a new internal node
 * 
 * \param[out] octree Pointer to the linear octree
 * \param[out] allocated_internal_nodes_ptr Pointer to the number of allocated internal nodes
 * \param[in] level Node level
 * \param[in] node Node index
 * \param[in] node_morton_index_level Morton index of the node at the current level
 * 
 * \return ErrorStatus
 */
IN_FILE ErrorStatus setup_node(
    LinearOctree *__restrict octree,
    int *__restrict allocated_internal_nodes_ptr,
    const int level,
    const int node,
    const int64 node_morton_index_level
)
{
    ErrorStatus error_status;

    /* Declare variables */
    int *__restrict num_internal_nodes_ptr = &octree->num_internal_nodes;
    int *__restrict tree_num_particles = octree->tree_num_particles;
    int *__restrict tree_num_internal_children = octree->tree_num_internal_children;
    int *__restrict tree_first_particle_sorted_idx = octree->tree_first_particle_sorted_idx;
    int *__restrict tree_first_internal_children_idx = octree->tree_first_internal_children_idx;

    double *__restrict tree_mass = octree->tree_mass;
    double *__restrict tree_center_of_mass_x = octree->tree_center_of_mass_x;
    double *__restrict tree_center_of_mass_y = octree->tree_center_of_mass_y;
    double *__restrict tree_center_of_mass_z = octree->tree_center_of_mass_z;

    int num_particles_per_octant[8] = {0};

    /* Find the number of particles in each octant */
    const int start_idx = tree_first_particle_sorted_idx[node];
    const int end_idx = start_idx + tree_num_particles[node] - 1;
    const int child_level = level + 1;
    error_status = WRAP_TRACEBACK(binary_search_num_particles_per_octant(
        num_particles_per_octant,
        octree->particle_morton_indices_deepest_level,
        node_morton_index_level,
        start_idx,
        end_idx,
        child_level
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Set up child nodes */
    bool first_child_found = false;
    int cumulative_count = 0;
    for (int i = 0; i < 8; i++)
    {
        if (num_particles_per_octant[i] <= 0)
        {
            continue;
        }

        const int child = *num_internal_nodes_ptr;

        // Reallocate memory if necessary
        if (child >= *allocated_internal_nodes_ptr)
        {
            *allocated_internal_nodes_ptr *= 2;
            int *tmp_tree_num_particles = realloc(tree_num_particles, *allocated_internal_nodes_ptr * sizeof(int));
            if (!tmp_tree_num_particles)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR,
                    "Failed to reallocate memory for tree_num_particles"
                );
            }
            tree_num_particles = tmp_tree_num_particles;
            octree->tree_num_particles = tree_num_particles;

            int *tmp_tree_num_internal_children = realloc(tree_num_internal_children, *allocated_internal_nodes_ptr * sizeof(int));
            if (!tmp_tree_num_internal_children)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR,
                    "Failed to reallocate memory for tree_num_internal_children"
                );
            }
            tree_num_internal_children = tmp_tree_num_internal_children;
            octree->tree_num_internal_children = tree_num_internal_children;

            int *tmp_tree_first_particle_sorted_idx = realloc(tree_first_particle_sorted_idx, *allocated_internal_nodes_ptr * sizeof(int));
            if (!tmp_tree_first_particle_sorted_idx)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR,
                    "Failed to reallocate memory for tree_first_particle_sorted_idx"
                );
            }
            tree_first_particle_sorted_idx = tmp_tree_first_particle_sorted_idx;
            octree->tree_first_particle_sorted_idx = tree_first_particle_sorted_idx;

            int *tmp_tree_first_internal_children_idx = realloc(tree_first_internal_children_idx, *allocated_internal_nodes_ptr * sizeof(int));
            if (!tmp_tree_first_internal_children_idx)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR,
                    "Failed to reallocate memory for tree_first_internal_children_idx"
                );
            }
            tree_first_internal_children_idx = tmp_tree_first_internal_children_idx;
            octree->tree_first_internal_children_idx = tree_first_internal_children_idx;

            double *tmp_tree_mass = realloc(octree->tree_mass, *allocated_internal_nodes_ptr * sizeof(double));
            if (!tmp_tree_mass)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR,
                    "Failed to reallocate memory for tree_mass"
                );
            }
            tree_mass = tmp_tree_mass;
            octree->tree_mass = tmp_tree_mass;
            
            double *tmp_tree_center_of_mass_x = realloc(octree->tree_center_of_mass_x, *allocated_internal_nodes_ptr * sizeof(double));
            if (!tmp_tree_center_of_mass_x)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR,
                    "Failed to reallocate memory for tree_center_of_mass_x"
                );
            }
            octree->tree_center_of_mass_x = tmp_tree_center_of_mass_x;
            tree_center_of_mass_x = tmp_tree_center_of_mass_x;

            double *tmp_tree_center_of_mass_y = realloc(octree->tree_center_of_mass_y, *allocated_internal_nodes_ptr * sizeof(double));
            if (!tmp_tree_center_of_mass_y)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR,
                    "Failed to reallocate memory for tree_center_of_mass_y"
                );
            }
            octree->tree_center_of_mass_y = tmp_tree_center_of_mass_y;
            tree_center_of_mass_y = tmp_tree_center_of_mass_y;

            double *tmp_tree_center_of_mass_z = realloc(octree->tree_center_of_mass_z, *allocated_internal_nodes_ptr * sizeof(double));
            if (!tmp_tree_center_of_mass_z)
            {
                return WRAP_RAISE_ERROR(
                    GRAV_MEMORY_ERROR,
                    "Failed to reallocate memory for tree_center_of_mass_z"
                );
            }
            octree->tree_center_of_mass_z = tmp_tree_center_of_mass_z;
            tree_center_of_mass_z = tmp_tree_center_of_mass_z;
        }

        if (!first_child_found)
        {
            first_child_found = true;
            tree_first_internal_children_idx[node] = child;
            tree_num_internal_children[node] = 0;
        }

        // Create a new internal node
        (*num_internal_nodes_ptr)++;
        (tree_num_internal_children[node])++;

        tree_num_internal_children[child] = 0;
        tree_num_particles[child] = num_particles_per_octant[i];
        tree_first_particle_sorted_idx[child] = start_idx + cumulative_count;

        tree_mass[child] = 0.0;
        tree_center_of_mass_x[child] = 0.0;
        tree_center_of_mass_y[child] = 0.0;
        tree_center_of_mass_z[child] = 0.0;

        cumulative_count += num_particles_per_octant[i];
    }

    return make_success_error_status();
}

/**
 * \brief Helper function to construct the octree
 * 
 * \param[out] octree Pointer to the linear octree
 * \param[in] allocated_internal_nodes Number of allocated internal nodes
 * \param[in] max_num_particles_per_leaf Maximum number of particles per leaf
 * \param[in] num_particles Number of particles
 * \param[in] x Array of position vectors
 * \param[in] m Array of masses
 * 
 * \return ErrorStatus
 */
IN_FILE ErrorStatus helper_construct_octree(
    LinearOctree *__restrict octree,
    int allocated_internal_nodes,
    const int max_num_particles_per_leaf,
    const int num_particles,
    const double *__restrict x,
    const double *__restrict m
)
{
    typedef struct Stack
    {
        int node;
        int processed_children;
        double total_mass;
        double mass_times_distance[3];
        struct Stack *parent;
    } Stack;

    ErrorStatus error_status;

    /* Create a stack */
    Stack stack[MORTON_MAX_LEVEL + 1];
    Stack *__restrict current_stack = &(stack[0]);

    current_stack->node = 0;
    current_stack->processed_children = -1;
    current_stack->total_mass = 0.0;
    current_stack->mass_times_distance[0] = 0.0;
    current_stack->mass_times_distance[1] = 0.0;
    current_stack->mass_times_distance[2] = 0.0;
    current_stack->parent = NULL;

    /* Declare variables */
    int *__restrict num_internal_nodes_ptr = &(octree->num_internal_nodes);
    const int64 *__restrict particle_morton_indices_deepest_level = octree->particle_morton_indices_deepest_level;
    const int *__restrict sorted_indices = octree->sorted_indices;

    /* Set up the root node */
    int level = 0;
    *num_internal_nodes_ptr = 1;

    octree->tree_num_particles[0] = num_particles;
    octree->tree_num_internal_children[0] = 0;
    octree->tree_first_particle_sorted_idx[0] = 0;
    octree->tree_mass[0] = 0.0;
    octree->tree_center_of_mass_x[0] = 0.0;
    octree->tree_center_of_mass_y[0] = 0.0;
    octree->tree_center_of_mass_z[0] = 0.0;

    error_status = WRAP_TRACEBACK(setup_node(
        octree,
        &allocated_internal_nodes,
        level,
        current_stack->node,
        0
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }
    level++;

    while (true)
    {
        const int current_node = current_stack->node;
        for (int i = current_stack->processed_children + 1; i < octree->tree_num_internal_children[current_node]; i++)
        {
            const int child = octree->tree_first_internal_children_idx[current_node] + i;
            const int start_idx = octree->tree_first_particle_sorted_idx[child];
            const int num_particles = octree->tree_num_particles[child];

            /* Leaf node */
            if (num_particles <= max_num_particles_per_leaf || level >= MORTON_MAX_LEVEL)
            {
                octree->tree_num_internal_children[child] = 0;

                // Update the stack
                for (int j = 0; j < num_particles; j++)
                {
                    const int particle_idx = sorted_indices[start_idx + j];
                    current_stack->total_mass += m[particle_idx];
                    current_stack->mass_times_distance[0] += m[particle_idx] * x[particle_idx * 3 + 0];
                    current_stack->mass_times_distance[1] += m[particle_idx] * x[particle_idx * 3 + 1];
                    current_stack->mass_times_distance[2] += m[particle_idx] * x[particle_idx * 3 + 2];
                }
                current_stack->processed_children = i;

                continue;
            }

            /* Internal node */
            else
            {
                const int64 child_morton_index_level = (particle_morton_indices_deepest_level[start_idx] >> (3 * (MORTON_MAX_LEVEL - level)));
                error_status = WRAP_TRACEBACK(setup_node(
                    octree,
                    &allocated_internal_nodes,
                    level,
                    child,
                    child_morton_index_level
                ));
                if (error_status.return_code != GRAV_SUCCESS)
                {
                    return error_status;
                }

                Stack *__restrict new_item = &(stack[level + 1]);
                new_item->node = child;
                new_item->processed_children = -1;
                new_item->total_mass = 0.0;
                new_item->mass_times_distance[0] = 0.0;
                new_item->mass_times_distance[1] = 0.0;
                new_item->mass_times_distance[2] = 0.0;
                new_item->parent = current_stack;

                current_stack = new_item;
                level++;

                break;
            }
        }

        /* Processed all children */
        if ((current_stack->processed_children + 1) >= octree->tree_num_internal_children[current_stack->node])
        {
            /* Update center of mass */
            octree->tree_mass[current_node] = current_stack->total_mass;
            octree->tree_center_of_mass_x[current_node] = current_stack->mass_times_distance[0] / current_stack->total_mass;
            octree->tree_center_of_mass_y[current_node] = current_stack->mass_times_distance[1] / current_stack->total_mass;
            octree->tree_center_of_mass_z[current_node] = current_stack->mass_times_distance[2] / current_stack->total_mass;

            Stack *parent = current_stack->parent;
            if (!parent)
            {
                break;
            }

            parent->total_mass += current_stack->total_mass;
            parent->mass_times_distance[0] += current_stack->mass_times_distance[0];
            parent->mass_times_distance[1] += current_stack->mass_times_distance[1];
            parent->mass_times_distance[2] += current_stack->mass_times_distance[2];
            
            current_stack = parent;
            (current_stack->processed_children)++;
            level--;
        }
    }

    /* Release unused memory */
    if (allocated_internal_nodes > (*num_internal_nodes_ptr))
    {
        int *__restrict tmp_tree_num_particles = realloc(octree->tree_num_particles, *num_internal_nodes_ptr * sizeof(int));
        if (!tmp_tree_num_particles)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to reallocate memory for tree_num_particles"
            );
        }
        octree->tree_num_particles = tmp_tree_num_particles;

        int *__restrict tmp_tree_num_internal_children = realloc(octree->tree_num_internal_children, *num_internal_nodes_ptr * sizeof(int));
        if (!tmp_tree_num_internal_children)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to reallocate memory for tree_num_internal_children"
            );
        }
        octree->tree_num_internal_children = tmp_tree_num_internal_children;

        int *__restrict tmp_tree_first_particle_sorted_idx = realloc(octree->tree_first_particle_sorted_idx, *num_internal_nodes_ptr * sizeof(int));
        if (!tmp_tree_first_particle_sorted_idx)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to reallocate memory for tree_first_particle_sorted_idx"
            );
        }
        octree->tree_first_particle_sorted_idx = tmp_tree_first_particle_sorted_idx;

        int *__restrict tmp_tree_first_internal_children_idx = realloc(octree->tree_first_internal_children_idx, *num_internal_nodes_ptr * sizeof(int));
        if (!tmp_tree_first_internal_children_idx)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to reallocate memory for tree_first_internal_children_idx"
            );
        }
        octree->tree_first_internal_children_idx = tmp_tree_first_internal_children_idx;

        double *__restrict tmp_tree_mass = realloc(octree->tree_mass, *num_internal_nodes_ptr * sizeof(double));
        if (!tmp_tree_mass)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to reallocate memory for tree_mass"
            );
        }
        octree->tree_mass = tmp_tree_mass;

        double *__restrict tmp_tree_center_of_mass_x = realloc(octree->tree_center_of_mass_x, *num_internal_nodes_ptr * sizeof(double));
        if (!tmp_tree_center_of_mass_x)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to reallocate memory for tree_center_of_mass_x"
            );
        }
        octree->tree_center_of_mass_x = tmp_tree_center_of_mass_x;

        double *__restrict tmp_tree_center_of_mass_y = realloc(octree->tree_center_of_mass_y, *num_internal_nodes_ptr * sizeof(double));
        if (!tmp_tree_center_of_mass_y)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to reallocate memory for tree_center_of_mass_y"
            );
        }
        octree->tree_center_of_mass_y = tmp_tree_center_of_mass_y;

        double *__restrict tmp_tree_center_of_mass_z = realloc(octree->tree_center_of_mass_z, *num_internal_nodes_ptr * sizeof(double));
        if (!tmp_tree_center_of_mass_z)
        {
            return WRAP_RAISE_ERROR(
                GRAV_MEMORY_ERROR,
                "Failed to reallocate memory for tree_center_of_mass_z"
            );
        }
        octree->tree_center_of_mass_z = tmp_tree_center_of_mass_z;
    }

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus construct_octree(
    LinearOctree *__restrict octree,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param,
    const double *__restrict box_center,
    const double box_width
)
{
    ErrorStatus error_status;

    /* Check for pointers */
    if (!octree)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Octree pointer is NULL");
    }
    if (!system)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "System pointer is NULL");
    }
    if (!acceleration_param)
    {
        return WRAP_RAISE_ERROR(GRAV_POINTER_ERROR, "Acceleration parameter pointer is NULL");
    }

    const int num_particles = system->num_particles;
    const double *__restrict x = system->x;
    const double *__restrict m = system->m;
    const int max_num_particles_per_leaf = acceleration_param->max_num_particles_per_leaf;

    /* Find the width and center of the bounding box */
    double center[3];
    if (!box_center || box_width <= 0.0)
    {
        calculate_bounding_box(center, &(octree->box_width), num_particles, x);
        box_center = center;
    }
    else
    {
        octree->box_width = box_width;
        center[0] = box_center[0];
        center[1] = box_center[1];
        center[2] = box_center[2];
    }

    /* Allocate memory */
    // Indices
    octree->particle_morton_indices_deepest_level = malloc(num_particles * sizeof(int64));
    octree->sorted_indices = malloc(num_particles * sizeof(int));
    if (!octree->particle_morton_indices_deepest_level || !octree->sorted_indices)
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for Morton indices and sorted indices"
        );
        goto err_indices_memory_alloc;
    }

    // Internal nodes
    // int allocated_internal_nodes = num_particles * 2 / max_num_particles_per_leaf;
    int allocated_internal_nodes = num_particles;

    octree->tree_num_particles = malloc(allocated_internal_nodes * sizeof(int));
    octree->tree_num_internal_children = malloc(allocated_internal_nodes * sizeof(int));
    octree->tree_first_internal_children_idx = malloc(allocated_internal_nodes * sizeof(int));
    octree->tree_first_particle_sorted_idx = malloc(allocated_internal_nodes * sizeof(int));
    octree->tree_mass = malloc(allocated_internal_nodes * sizeof(double));
    octree->tree_center_of_mass_x = malloc(allocated_internal_nodes * sizeof(double));
    octree->tree_center_of_mass_y = malloc(allocated_internal_nodes * sizeof(double));
    octree->tree_center_of_mass_z = malloc(allocated_internal_nodes * sizeof(double));
    if (
        !octree->tree_num_particles ||
        !octree->tree_num_internal_children ||
        !octree->tree_first_internal_children_idx ||
        !octree->tree_first_particle_sorted_idx ||
        !octree->tree_mass ||
        !octree->tree_center_of_mass_x ||
        !octree->tree_center_of_mass_y ||
        !octree->tree_center_of_mass_z
    )
    {
        error_status = WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for internal nodes"
        );
        goto err_internal_nodes_memory_alloc;
    }

    /* Initialize the sorted indices */
    for (int i = 0; i < num_particles; i++)
    {
        octree->sorted_indices[i] = i;
    }

    /* Compute the 3D Morton indices at level 21 */
    compute_3d_particle_morton_indices_deepest_level(
        octree->particle_morton_indices_deepest_level,
        num_particles,
        x,
        center,
        octree->box_width
    );

    /* Sort the particles based on their Morton indices */
    error_status = WRAP_TRACEBACK(radix_sort_particles_morton_index(
        octree->particle_morton_indices_deepest_level,
        octree->sorted_indices,
        num_particles,
        MORTON_MAX_LEVEL
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto err_radix_sort;
    }

    /* Construct the octree */
    error_status = WRAP_TRACEBACK(helper_construct_octree(
        octree,
        allocated_internal_nodes,
        max_num_particles_per_leaf,
        num_particles,
        x,
        m
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        goto err_construct_octree;
    }

    return make_success_error_status();

err_construct_octree:
err_radix_sort:
err_internal_nodes_memory_alloc:
err_indices_memory_alloc:
    free_linear_octree(octree);

    return error_status;
}

WIN32DLL_API void free_linear_octree(LinearOctree *__restrict octree)
{
    free(octree->particle_morton_indices_deepest_level);
    free(octree->sorted_indices);
    free(octree->tree_num_particles);
    free(octree->tree_num_internal_children);
    free(octree->tree_first_particle_sorted_idx);
    free(octree->tree_first_internal_children_idx);
    free(octree->tree_mass);
    free(octree->tree_center_of_mass_x);
    free(octree->tree_center_of_mass_y);
    free(octree->tree_center_of_mass_z);
}

WIN32DLL_API bool linear_octree_check_if_included(
    const int64 morton_index_i,
    const int64 morton_index_j,
    const int level
)
{
    return (morton_index_i >> (3 * (MORTON_MAX_LEVEL - level))) == (morton_index_j >> (3 * (MORTON_MAX_LEVEL - level)));
}
