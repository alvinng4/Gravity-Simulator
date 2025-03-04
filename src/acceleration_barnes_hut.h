#ifndef ACCELERATION_BARNES_HUT_H
#define ACCELERATION_BARNES_HUT_H

#include "gravity_sim.h"

#define MAX_NUM_PARTICLES_PER_LEAF 1 // Note: Potential optimization, fix this to 1 we can remove some loops and we may not need tree_num_particles array
#define MORTON_MAX_LEVEL 21 // Maximum level for 64-bit Morton index, don't change



/**
 * \brief Compute acceleration with Barnes-Hut algorithm
 * 
 * \param a Array of acceleration vectors to be modified
 * \param system Pointer to the gravitational system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \retval SUCCESS If the computation is successful
 * \retval error code if errors occurred
 */
int acceleration_barnes_hut(
    real *__restrict a,
    const System *__restrict system,
    const AccelerationParam *__restrict acceleration_param
);

#ifdef __cplusplus
    extern "C" {
#endif
    int barnes_hut_setup_octree(
        real *__restrict width,
        int *__restrict allocated_internal_nodes,
        int *__restrict actual_num_internal_nodes,
        const int objects_count,
        const real *__restrict x,
        const real *__restrict m,
        int64 **leaf_morton_indices_deepest_level,
        int **sorted_indices,
        int **tree_start_particle_sorted_idx,
        int **tree_num_particles,
        int **tree_num_internal_children,
        int **tree_idx_first_internal_child,
        real **tree_total_mass,
        real **tree_center_of_mass_x,
        real **tree_center_of_mass_y,
        real **tree_center_of_mass_z
    );
#ifdef __cplusplus
    }
#endif










#endif
