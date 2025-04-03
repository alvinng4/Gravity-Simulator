#include <stdbool.h>

#include "settings.h"

Settings get_new_settings(void)
{
    Settings settings;
    settings.verbose = GRAV_VERBOSITY_NORMAL;
    settings.enable_progress_bar = true;
    settings.remove_invalid_particles = true;
    settings.silence_remove_invalid_particles = false;
    settings.is_exit = false;
    return settings;
}
