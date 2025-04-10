#include <stdio.h>

#include "common.h"
#include "grav_sim.h"

#ifdef USE_OPENMP
    #include <omp.h>
#endif

#ifdef USE_HDF5
    #include <hdf5.h>
#endif

#ifdef USE_FFTW3
    #include <fftw3.h>
#endif

WIN32DLL_API void print_compilation_info(void)
{
    const char *new_line = "\n";
    const char *straight_line = "-----------------------------------------------------------------\n";

    fputs(straight_line, stdout);
    fputs(get_grav_sim_logo_string(), stdout);
    fputs(new_line, stdout);
    fputs(new_line, stdout);

    fputs("grav_sim\n", stdout);
    fputs(new_line, stdout);

    /* Print OS information */
#ifdef _WIN32
    fputs("Operating System: Windows\n", stdout);
#elif __APPLE__
    fputs("Operating System: MacOS\n", stdout);
#elif __linux__
    fputs("Operating System: Linux\n", stdout);
#else
    fputs("Operating System: Unknown\n", stdout);
#endif

    /* Print compilation information */
    fputs("Compilation Info:\n", stdout);

    /* OpenMP */
#ifdef USE_OPENMP
    fputs("  Compiled with OpenMP: true\n", stdout);
#else
    fputs("  Compiled with OpenMP: false\n", stdout);
#endif

    /* HDF5 */
#ifdef USE_HDF5
    fputs("  Compiled with HDF5: true\n", stdout);
    #ifdef H5_VERS_MAJOR
        printf("    Version: %d.%d.%d\n", H5_VERS_MAJOR, H5_VERS_MINOR, H5_VERS_RELEASE);
    #endif
#else
    fputs("  Compiled with HDF5: false\n", stdout);
#endif

    /* FFTW3 */
#ifdef USE_FFTW3
    fputs("  Compiled with FFTW3: true\n", stdout);
    printf("    Version: %s\n", fftw_version);;
#else
    fputs("  Compiled with FFTW3: false\n", stdout);
#endif

    fputs(new_line, stdout);

    /* Build and compiler info */
    printf("Build time: %s %s\n", __DATE__, __TIME__);

#ifdef _MSC_VER
    printf("Compiler: MSVC (version: %d)\n", _MSC_VER);
#elif defined(__clang__)
    printf("Compiler: Clang (version: %d)\n", __clang_major__);
#elif defined(__GNUC__)
    printf("Compiler: GCC (version: %d)\n", __GNUC__);
#endif
    fputs(straight_line, stdout);
}
