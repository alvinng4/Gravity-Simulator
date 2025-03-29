#ifndef SYSTEM_H
#define SYSTEM_H

#include "error.h"

typedef struct System
{
    int objects_count;
    int *particle_ids;
    double *x;
    double *v;
    double *m;
    double G;
} System;

System get_new_system(void);


/**
 * \brief Get a new system with initialized memory for the given number of objects
 * 
 * \param[out] system Pointer to the system to be initialized
 * \param[in] objects_count Number of objects
 * 
 * \return ErrorStatus
 */
ErrorStatus get_initialized_system(
    System *__restrict system,
    const int objects_count
);

void free_system(System *__restrict system);

ErrorStatus initialize_built_in_system(
    System *__restrict system,
    const char *__restrict system_name
);







#endif
