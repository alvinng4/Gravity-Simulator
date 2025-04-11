cmake_minimum_required(VERSION 3.10)

# Help finding the libraries
cmake_policy(SET CMP0074 NEW)

project(grav_sim LANGUAGES C)

# Set the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# List of source files
set(SRCS
    acceleration.c
    acceleration_barnes_hut.c
    cosmology.c
    error.c
    grav_sim.c
    integrator.c
    integrator_rk_embedded.c
    integrator_ias15.c
    integrator_whfast.c
    linear_octree.c
    math_functions.c
    output.c
    progress_bar.c
    python_interface.c
    settings.c
    system.c
    utils.c
)

set(FFTW3_SRCS
    acceleration_PM.c
    integrator_cosmology_leapfrog.c
)

# Create the library target
add_library(grav_sim SHARED ${SRCS})

# Use the current source directory as a public include for this target
target_include_directories(grav_sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/pcg)

# Debug mode
option(DEBUG "Debug mode" OFF)

if (DEBUG)
    target_compile_options(
        grav_sim PRIVATE
        -O0 -g -Wall -Wextra -Wpedantic -fsanitize=address -fsanitize=undefined
    )
else()
    target_compile_options(
        grav_sim PRIVATE
        -O3 -Wall -Wextra -Wpedantic
    )
endif()

# Optionally enable OpenMP
option(USE_OPENMP "Enable OpenMP" OFF)
if (USE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_C_FOUND)
        target_link_libraries(grav_sim PRIVATE OpenMP::OpenMP_C)
        target_compile_definitions(grav_sim PRIVATE USE_OPENMP)
        message(STATUS "Building with OpenMP support")
    else()
        message(WARNING "OpenMP requested but not found")
    endif()
endif()

# Optionally enable HDF5 support
option(USE_HDF5 "Enable HDF5" OFF)
if (USE_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS C)
    if (HDF5_FOUND)
        target_include_directories(grav_sim PRIVATE ${HDF5_INCLUDE_DIRS})
        target_compile_definitions(grav_sim PRIVATE USE_HDF5)
        target_link_libraries(grav_sim PRIVATE ${HDF5_LIBRARIES})
        message(STATUS "Building with HDF5 support")
        message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
    else()
        message(WARNING "HDF5 requested but not found")
    endif()
endif()

# Optionally enable FFTW support
option(USE_FFTW3 "Enable FFTW3" OFF)
if (USE_FFTW3)
    find_package(FFTW3 REQUIRED)
    if (FFTW3_FOUND)
        target_sources(grav_sim PRIVATE ${FFTW3_SRCS})
        target_include_directories(grav_sim PRIVATE ${FFTW3_INCLUDE_DIRS})
        target_compile_definitions(grav_sim PRIVATE USE_FFTW3)
        target_link_libraries(grav_sim PRIVATE ${FFTW3_LIBRARIES})
        message(STATUS "Building with FFTW3 support")
    else()
        message(WARNING "FFTW3 requested but not found")
    endif()
endif()

# Add the PCG library as a subdirectory
add_subdirectory(pcg)
target_link_libraries(grav_sim PRIVATE pcg)

if(UNIX AND NOT APPLE)
    target_link_libraries(grav_sim PRIVATE rt)
endif()

# Version info
if(DEFINED VERSION_INFO)
    target_compile_definitions(grav_sim PRIVATE VERSION_INFO="${VERSION_INFO}")
else()
    target_compile_definitions(grav_sim PRIVATE VERSION_INFO="unknown")
endif()

target_link_libraries(grav_sim PRIVATE m)

# Install rules
install(TARGETS grav_sim
    ARCHIVE DESTINATION lib
)
