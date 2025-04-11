#include <stdbool.h>

#include "settings.h"

Settings get_new_settings(void)
{
    static bool is_exit = false;

    Settings settings;
    settings.verbose = GRAV_VERBOSITY_NORMAL;
    settings.enable_progress_bar = true;
    settings.is_exit_ptr = &is_exit;
    return settings;
}
