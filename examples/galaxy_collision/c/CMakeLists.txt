cmake_minimum_required(VERSION 3.10)
project(galaxy_collision C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra -Wpedantic")

# Define source and library directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../src)
include_directories(${SRC_DIR})

# Debugging flags
option (DEBUG "Debug mode" OFF)
if (DEBUG)
    set(CMAKE_C_FLAGS "-O0 -g -Wall -Wextra -Wpedantic -fsanitize=address -fsanitize=undefined")
endif()

# Define the executable
add_executable(galaxy_collision galaxy_collision.c)

# HDF5 library
find_package(HDF5 REQUIRED COMPONENTS C)
target_include_directories(galaxy_collision PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(galaxy_collision PRIVATE ${HDF5_LIBRARIES})

# Link against the grav_sim library
add_subdirectory(${SRC_DIR} build)
target_link_libraries(galaxy_collision PRIVATE grav_sim)
