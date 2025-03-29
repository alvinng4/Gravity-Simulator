cmake_minimum_required(VERSION 3.10)
project(gravity_sim LANGUAGES C)

# Set the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# List of source files
set(SRCS
    acceleration.c
    acceleration_barnes_hut.c
    error.c
    gravity_sim.c
    integrator.c
    integrator_rk_embedded.c
    integrator_ias15.c
    integrator_whfast.c
    linear_octree.c
    math_functions.c
    output.c
    progress_bar.c
    settings.c
    system.c
    utils.c
)

# Create the static library target
add_library(gravity_sim STATIC ${SRCS})

# Use the current source directory as a public include for this target
target_include_directories(gravity_sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Optionally toggle Debug mode (custom flag)
option(DEBUG "Debug mode" OFF)

if (DEBUG)
    target_compile_options(
        gravity_sim PRIVATE
        -O0 -g -Wall -Wextra -Wpedantic -fsanitize=address -fsanitize=undefined
    )
else()
    target_compile_options(
        gravity_sim PRIVATE
        -O3 -Wall -Wextra -Wpedantic
    )
endif()

# Optionally enable OpenMP
option(USE_OPENMP "Enable OpenMP" OFF)
if (USE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_C_FOUND)
        target_link_libraries(gravity_sim PRIVATE OpenMP::OpenMP_C)
        target_compile_definitions(gravity_sim PRIVATE USE_OPENMP)
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
        target_include_directories(gravity_sim PRIVATE ${HDF5_INCLUDE_DIRS})
        target_compile_definitions(gravity_sim PRIVATE USE_HDF5)
        target_link_libraries(gravity_sim PRIVATE ${HDF5_LIBRARIES})
        message(STATUS "Building with HDF5 support")
    else()
        message(WARNING "HDF5 requested but not found")
    endif()
endif()

# Install rules
install(TARGETS gravity_sim
    ARCHIVE DESTINATION lib
)
