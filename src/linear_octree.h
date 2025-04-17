/**
 * \file linear_octree.h
 * \brief Linear octree for Barnes-Hut algorithm
 * 
 * \author Ching-Yin Ng
 */

#ifndef LINEAR_OCTREE_H
#define LINEAR_OCTREE_H

#include "common.h"

/* Maximum level for 64-bit Morton index, do not change */
#define MORTON_MAX_LEVEL 21


/**
 * \brief Linear octree structure 
 */
typedef struct LinearOctree
{
    /** \brief The width of the bounding box. */
    double box_width;

    /** \brief The number of internal nodes in the octree. */
    int num_internal_nodes;

    /** \brief Array of particles' Morton indices at the deepest level of the octree. */
    int64 *particle_morton_indices_deepest_level;

    /** \brief Array of sorted indices of particles. */
    int *sorted_indices;

    /** \brief Array storing the number of particles in each node. */
    int *tree_num_particles;

    /** \brief Array storing the number of internal children for each internal node. */
    int *tree_num_internal_children;

    /** \brief Array storing the sorted index of the first particle in each node. */
    int *tree_first_particle_sorted_idx;

    /** \brief Array storing the index of the first internal child for each node. */
    int *tree_first_internal_children_idx;

    /** \brief Array storing the total mass of particles contained within each node. */
    double *tree_mass;

    /** \brief Array storing the x-component of the center of mass for each node. */
    double *tree_center_of_mass_x;

    /** \brief Array storing the y-component of the center of mass for each node. */
    double *tree_center_of_mass_y;

    /** \brief Array storing the z-component of the center of mass for each node. */
    double *tree_center_of_mass_z;
} LinearOctree;

/**
 * \brief Get a new linear octree struct
 * 
 * \return LinearOctree
 */
LinearOctree get_new_linear_octree(void);

/**
 * \brief Construct the linear octree
 * 
 * \param octree Pointer to the linear octree
 * \param system Pointer to the system
 * \param acceleration_param Pointer to the acceleration parameters
 * 
 * \return ErrorStatus
 */
ErrorStatus construct_octree(
    LinearOctree *restrict octree,
    const System *restrict system,
    const AccelerationParam *restrict acceleration_param,
    const double *restrict box_center,
    const double box_width
);

/**
 * \brief Free the memory allocated for the linear octree
 * 
 * \param octree Pointer to the linear octree
 */
void free_linear_octree(LinearOctree *restrict octree);

/**
 * \brief Check if two Morton indices are included in the same octant
 * 
 * \param morton_index_i Morton index of the first object at the deepest level
 * \param morton_index_j Morton index of the second object at the deepest level
 * \param level Level of the Morton indices
 */
bool linear_octree_check_if_included(
    const int64 morton_index_i,
    const int64 morton_index_j,
    const int level
);

#endif
