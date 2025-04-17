/**
 * \file acceleration_barnes_hut.c
 * \brief Implementation of Barnes-Hut algorithm
 * 
 * \author Ching-Yin Ng
 */

#include <math.h>
#include <stdio.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "acceleration.h"
#include "linear_octree.h"

/**
 * \brief Helper function to compute acceleration
 * 
 * \param[out] a Array of acceleration vectors to be modified
 * \param[in] system Pointer to the gravitational system
 * \param[in] acceleration_param Pointer to the acceleration parameters
 * \param[in] octree Pointer to the linear octree
 */
IN_FILE void helper_compute_acceleration(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param,
    const LinearOctree *restrict octree
);

WIN32DLL_API ErrorStatus acceleration_barnes_hut(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param
)
{
    ErrorStatus error_status;

    /* Empty the input array */
    const int num_particles = system->num_particles;
    for (int i = 0; i < num_particles; i++)
    {
        a[i * 3 + 0] = 0.0;
        a[i * 3 + 1] = 0.0;
        a[i * 3 + 2] = 0.0;
    }

    /* Construct octree */
    LinearOctree octree = get_new_linear_octree();
    error_status = WRAP_TRACEBACK(construct_octree(
        &octree,
        system,
        acceleration_param,
        NULL,
        -1.0
    ));
    if (error_status.return_code != GRAV_SUCCESS)
    {
        return error_status;
    }

    /* Compute acceleration */
    helper_compute_acceleration(
        a,
        system,
        acceleration_param,
        &octree
    );

    /* Free memory */
    free_linear_octree(&octree);

    return make_success_error_status();
}

IN_FILE void helper_compute_acceleration(
    double *restrict a,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param,
    const LinearOctree *restrict octree
)
{
    typedef struct Stack
    {
        int node;
        int processed_children;
        struct Stack *parent;
    } Stack;

    /* Declare variables */
    const int num_particles = system->num_particles;
    const double *restrict x = system->x;
    const double *restrict m = system->m;
    const double G = system->G;
    const double softening_length = acceleration_param->softening_length;
    const double softening_length_squared = softening_length * softening_length;
    const double opening_angle = acceleration_param->opening_angle;
    const double opening_angle_squared = opening_angle * opening_angle;

    const double box_length = octree->box_width * 2.0;
    const int64 *restrict particle_morton_indices_deepest_level = octree->particle_morton_indices_deepest_level;
    const int *restrict sorted_indices = octree->sorted_indices;
    const int *restrict tree_num_particles = octree->tree_num_particles;
    const int *restrict tree_num_internal_children = octree->tree_num_internal_children;
    const int *restrict tree_first_particle_sorted_idx = octree->tree_first_particle_sorted_idx;
    const int *restrict tree_first_internal_children_idx = octree->tree_first_internal_children_idx;
    const double *restrict tree_mass = octree->tree_mass;
    const double *restrict tree_center_of_mass_x = octree->tree_center_of_mass_x;
    const double *restrict tree_center_of_mass_y = octree->tree_center_of_mass_y;
    const double *restrict tree_center_of_mass_z = octree->tree_center_of_mass_z;

#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < num_particles; i++)
    {
        const int idx_i = sorted_indices[i];    // For coalesced memory access
        const int64 morton_index_i = particle_morton_indices_deepest_level[idx_i];
        const double x_i[3] = {x[idx_i * 3 + 0], x[idx_i * 3 + 1], x[idx_i * 3 + 2]};
    
        Stack stack[MORTON_MAX_LEVEL + 1];
        Stack *current_stack = &(stack[0]);
        current_stack->processed_children = -1;
        current_stack->node = 0;
        current_stack->parent = NULL;
        double acceleration[3] = {0.0, 0.0, 0.0};

        int level = 1;

        /* Tree walk */
        while (true)
        {
            const int current_node = current_stack->node;
            for (int j = (current_stack->processed_children) + 1; j < tree_num_internal_children[current_node]; j++)
            {
                const int child_j = tree_first_internal_children_idx[current_node] + j;
                const int num_children_j = tree_num_internal_children[child_j];
                const int start_idx_j = tree_first_particle_sorted_idx[child_j];

                // If object i is included, then we need to traverse deeper
                const bool is_included = linear_octree_check_if_included(
                    morton_index_i,
                    particle_morton_indices_deepest_level[sorted_indices[start_idx_j]],
                    level
                );

                // Check Barnes-Hut criteria
                if (!is_included)
                {
                    const double R[3] = {
                        x_i[0] - tree_center_of_mass_x[child_j],
                        x_i[1] - tree_center_of_mass_y[child_j],
                        x_i[2] - tree_center_of_mass_z[child_j]
                    };
                    const double box_length_j = box_length / (2 << level);
                    const double norm_square = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];

                    // Check if box_length_j / norm < opening_angle
                    // Use squared values to avoid sqrt
                    if ((box_length_j * box_length_j) < opening_angle_squared * norm_square)
                    {
                        const double R_norm = sqrt(
                            norm_square + softening_length_squared
                        );

                        const double temp_value = G * tree_mass[child_j] / (R_norm * R_norm * R_norm);
                        acceleration[0] -= temp_value * R[0];
                        acceleration[1] -= temp_value * R[1];
                        acceleration[2] -= temp_value * R[2];
    
                        current_stack->processed_children = j;
                        continue;
                    }
                }

                /* Traverse deeper */

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
                        const double R[3] = {
                            x_i[0] - x[idx_j * 3 + 0],
                            x_i[1] - x[idx_j * 3 + 1],
                            x_i[2] - x[idx_j * 3 + 2]
                        };
                        const double R_norm = sqrt(
                            R[0] * R[0] + 
                            R[1] * R[1] + 
                            R[2] * R[2] +
                            softening_length_squared
                        );

                        // Calculate the acceleration
                        const double temp_value = G * m[idx_j] / (R_norm * R_norm * R_norm);
                        acceleration[0] -= temp_value * R[0];
                        acceleration[1] -= temp_value * R[1];
                        acceleration[2] -= temp_value * R[2];
                    }

                    current_stack->processed_children = j;
                    continue;
                }

                /* Internal node */
                else
                {
                    Stack *new_item = &(stack[level + 1]);
                    new_item->node = child_j;
                    new_item->processed_children = -1;
                    new_item->parent = current_stack;

                    current_stack = new_item;
                    level++;
                    break;
                }
            }

            if ((current_stack->processed_children + 1) >= tree_num_internal_children[current_stack->node])
            {
                Stack *parent_stack = current_stack->parent;
                if (!parent_stack)
                {
                    break;
                }

                current_stack = parent_stack;
                current_stack->processed_children += 1;
                level--;
            }
        }

        a[idx_i * 3 + 0] = acceleration[0];
        a[idx_i * 3 + 1] = acceleration[1];
        a[idx_i * 3 + 2] = acceleration[2];
    }
}
