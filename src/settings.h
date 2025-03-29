#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdbool.h>

#define GRAV_VERBOSITY_IGNORE_ALL 0
#define GRAV_VERBOSITY_IGNORE_INFO 1
#define GRAV_VERBOSITY_NORMAL 2

typedef struct Settings
{
    int verbose;
    bool enable_progress_bar;
    bool is_exit;
} Settings;

/**
 * \brief Get a new settings struct
 * 
 * \return Settings
 */
Settings get_new_settings(void);

#endif
