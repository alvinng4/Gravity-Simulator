cmake_minimum_required(VERSION 3.14)

project(galaxy_collision
    VERSION 1.0.1
    LANGUAGES C
    DESCRIPTION "Galaxy collision example using grav_sim"
)

### Options ###
option(ENABLE_SANITIZERS
    "Enable address and undefined-behavior sanitizers"
    OFF
)

### grav_sim ###
# Directory containing grav_sim's top-level CMakeLists.txt.
set(GRAV_SIM_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../../.."
    CACHE PATH
    "Path to the grav_sim source directory"
)

if(NOT EXISTS "${GRAV_SIM_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
        "grav_sim was not found at:\n"
        "  ${GRAV_SIM_SOURCE_DIR}\n"
        "Set GRAV_SIM_SOURCE_DIR to the grav_sim source directory."
    )
endif()

add_subdirectory(
    "${GRAV_SIM_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/grav_sim"
)

### Executable ###
add_executable(galaxy_collision
    galaxy_collision.c
)

set_target_properties(galaxy_collision PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS OFF
)

### Dependencies ###
target_link_libraries(galaxy_collision PRIVATE
    grav_sim::grav_sim
)

### HDF5 ###
find_package(HDF5 REQUIRED COMPONENTS C)

if(TARGET HDF5::HDF5)
    target_link_libraries(galaxy_collision PRIVATE HDF5::HDF5)
elseif(TARGET hdf5::hdf5)
    target_link_libraries(galaxy_collision PRIVATE hdf5::hdf5)
else()
    target_include_directories(galaxy_collision PRIVATE
        ${HDF5_INCLUDE_DIRS}
    )

    target_link_libraries(galaxy_collision PRIVATE
        ${HDF5_LIBRARIES}
    )

    if(HDF5_DEFINITIONS)
        target_compile_options(galaxy_collision PRIVATE
            ${HDF5_DEFINITIONS}
        )
    endif()
endif()

### Warnings ###
if(MSVC)
    target_compile_options(galaxy_collision PRIVATE
        /W4
    )
elseif(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(galaxy_collision PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wpointer-arith
        -Wcast-qual
        -Wstrict-prototypes
        -Wmissing-prototypes
    )
endif()

### Sanitizers ###
if(ENABLE_SANITIZERS)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(galaxy_collision PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
        )

        target_link_options(galaxy_collision PRIVATE
            -fsanitize=address,undefined
        )
    else()
        message(WARNING
            "ENABLE_SANITIZERS is ON, but sanitizers are not configured "
            "for compiler '${CMAKE_C_COMPILER_ID}'"
        )
    endif()
endif()
