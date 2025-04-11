#ifndef PYTHON_INTERFACE_H
#define PYTHON_INTERFACE_H

#include "common.h"


/**
 * \brief Free memory allocated for a int32 array
 * 
 * \param ptr Pointer to the int32 array
 */
void free_memory_int32(int32 *__restrict ptr);

/**
 * \brief Free memory allocated for a double array
 * 
 * \param ptr Pointer to the double array
 */
void free_memory_double(double *__restrict ptr);

/**
 * \brief Load a built-in system
 * 
 * \param system_name Name of the system
 * \param num_particles_ptr Pointer to the number of particles
 * \param particle_ids_ptr Pointer to the particle ID array
 * \param x_ptr Pointer to the coordinate array
 * \param v_ptr Pointer to the velocity array
 * \param m_ptr Pointer to the mass array
 * \param G_ptr Pointer to the gravitational constant
 * 
 * \retval 0 Success
 * \retval 1 Failure
 */
int32 load_built_in_system_python(
    const char *system_name,
    int *num_particles_ptr,
    int32 **particle_ids_ptr,
    double **x_ptr,
    double **v_ptr,
    double **m_ptr,
    double *G_ptr
);

/**
 * \brief Convert Keplerian elements to Cartesian coordinates
 * 
 * \param[out] x Output position x-coordinate
 * \param[out] y Output position y-coordinate
 * \param[out] z Output position z-coordinate
 * \param[out] v_x Output velocity x-coordinate
 * \param[out] v_y Output velocity y-coordinate
 * \param[out] v_z Output velocity z-coordinate
 * \param[in] semi_major_axis Semi-major axis
 * \param[in] eccentricity Eccentricity
 * \param[in] inclination Inclination
 * \param[in] argument_of_periapsis Argument of periapsis
 * \param[in] longitude_of_ascending_node Longitude of ascending node
 * \param[in] true_anomaly True anomaly
 * \param[in] total_mass Total mass of the system
 * \param[in] G Gravitational constant
 */
void keplerian_to_cartesian_python(
    double *__restrict x,
    double *__restrict y,
    double *__restrict z,
    double *__restrict v_x,
    double *__restrict v_y,
    double *__restrict v_z,
    const double semi_major_axis,
    const double eccentricity,
    const double inclination,
    const double argument_of_periapsis,
    const double longitude_of_ascending_node,
    const double true_anomaly,
    const double total_mass,
    const double G
);

/**
 * \brief Launch simulation from Python
 * 
 * \retval 0 Success
 * \retval 1 Failure
 */
int launch_simulation_python(
    int32 *__restrict num_particles,
    int32 *particle_ids,
    double *x,
    double *v,
    double *m,
    int32 **new_particle_ids_ptr,
    double **new_x_ptr,
    double **new_v_ptr,
    double **new_m_ptr,
    const double G,
    const int32 integrator,
    const double dt,
    const double tolerance,
    const double initial_dt,
    const bool whfast_remove_invalid_particles,
    const int32 acceleration_method,
    const double opening_angle,
    const double softening_length,
    const int32 max_num_particles_per_leaf,
    const int32 output_method,
    char *output_dir,
    const bool output_initial,
    const double output_interval,
    const int32 coordinate_output_dtype,
    const int32 velocity_output_dtype,
    const int32 mass_output_dtype,
    const int32 verbose,
    const bool enable_progress_bar,
    bool *is_exit_ptr,
    const double tf
);

/**
 * \brief Compute the energy from solution state
 * 
 * \param[out] energy Pointer to the energy array to be updated
 * \param[in] G Gravitational constant
 * \param[in] sol_state Pointer to the solution state array
 * \param[in] num_snapshots Number of time steps
 * \param[in] num_particles Number of particles in the system
 */
void compute_energy_python(
    double *__restrict energy,
    const double G,
    const double *__restrict sol_state,
    const int32 num_snapshots,
    const int32 num_particles
);

/**
 * \brief Compute the linear momentum from solution state
 * 
 * \param[out] linear_momentum Pointer to the linear momentum array to be updated
 * \param[in] sol_state Pointer to the solution state array
 * \param[in] num_snapshots Number of time steps
 * \param[in] num_particles Number of particles in the system
 */
void compute_linear_momentum_python(
    double *__restrict linear_momentum,
    const double *__restrict sol_state,
    const int32 num_snapshots,
    const int32 num_particles
);

/**
 * \brief Compute the angular momentum from solution state
 * 
 * \param[out] angular_momentum Pointer to the angular momentum array to be updated
 * \param[in] sol_state Pointer to the solution state array
 * \param[in] num_snapshots Number of time steps
 * \param[in] num_particles Number of particles in the system
 */
void compute_angular_momentum_python(
    double *__restrict angular_momentum,
    const double *__restrict sol_state,
    const int32 num_snapshots,
    const int32 num_particles
);

#endif
