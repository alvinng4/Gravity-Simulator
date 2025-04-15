cmake_minimum_required(VERSION 3.10)

set(PCG_SOURCES
    pcg_basic.c
)

add_library(pcg STATIC ${PCG_SOURCES})

target_include_directories(pcg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
