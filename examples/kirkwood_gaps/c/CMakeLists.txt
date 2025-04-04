cmake_minimum_required(VERSION 3.10)
project(kirkwood_gaps C)

# Set the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Define source and library directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../src)
include_directories(${SRC_DIR})

# Define the executable
add_executable(kirkwood_gaps kirkwood_gaps.c)

# Debug mode
option(DEBUG "Debug mode" OFF)

if (DEBUG)
    target_compile_options(
        kirkwood_gaps PRIVATE
        -O0 -g -Wall -Wextra -Wpedantic -fsanitize=address -fsanitize=undefined
    )
else()
    target_compile_options(
        kirkwood_gaps PRIVATE
        -O3 -Wall -Wextra -Wpedantic
    )
endif()

# Link against the gravity_sim library
add_subdirectory(${SRC_DIR} build)
target_link_libraries(kirkwood_gaps gravity_sim)
