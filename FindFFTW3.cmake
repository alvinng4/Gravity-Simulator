# - Try to find FFTW
#
# By default, it will look only for the serial libraries with single, double,
# and long double precision. Any combination of precision (SINGLE, DOUBLE,
# LONGDOUBLE) and library type (SERIAL, [THREADS|OPENMP], MPI) is possible by
# using the COMPONENTS keyword. For example,
#
# find_package(FFTW3 COMPONENTS SINGLE DOUBLE OPENMP MPI)
#
# Once done this will define
#  FFTW3_FOUND - System has FFTW3
#  FFTW3_INCLUDE_DIRS - The FFTW3 include directories
#  FFTW3_LIBRARIES - The libraries needed to use FFTW3
#  FFTW3_DEFINITIONS - Compiler switches required for using FFTW3
#  FFTW3_$KIND_$PARALLEL_FOUND- Set if FFTW3 exists in KIND precision format for PARALLEL mode.
#                             where KIND can be: SINGLE, DOUBLE, LONGDOUBLE
#                             and PARALLEL: SERIAL, OPENMP, MPI, THREADS.
#  FFTW3_$KIND_$PARALLEL_LIBRARY - The libraries needed to use.
#  FFTW3_INCLUDE_DIR_PARALLEL - The FFTW3 include directories for parallels mode.

cmake_policy(SET CMP0054 NEW)

if(FFTW3_FOUND)
  return()
endif()

if(FFTW3_INCLUDE_DIR AND FFTW3_LIBRARIES)
  set(FFTW3_FOUND TRUE)
  foreach(component ${FFTW3_FIND_COMPONENTS})
    if("${FFTW3_${component}_LIBRARY}" STREQUAL "")
        set(FFTW3_${component}_LIBRARY "${FFTW3_LIBRARIES}")
    endif()
  endforeach()
  return()
endif()

macro(find_specific_libraries KIND PARALLEL)
  list(APPEND FFTW3_FIND_COMPONENTS ${KIND}_${PARALLEL})
  if(NOT (${PARALLEL} STREQUAL "SERIAL") AND NOT ${PARALLEL}_FOUND)
    message(FATAL_ERROR "Please, find ${PARALLEL} libraries before FFTW")
  endif()

  find_library(FFTW3_${KIND}_${PARALLEL}_LIBRARY NAMES
    fftw3${SUFFIX_${KIND}}${SUFFIX_${PARALLEL}}${SUFFIX_FINAL} HINTS ${HINT_DIRS})
  if(FFTW3_${KIND}_${PARALLEL}_LIBRARY MATCHES fftw3)
    list(APPEND FFTW3_LIBRARIES ${FFTW3_${KIND}_${PARALLEL}_LIBRARY})
    set(FFTW3_${KIND}_${PARALLEL}_FOUND TRUE)

    STRING(TOLOWER "${KIND}" kind)
    STRING(TOLOWER "${PARALLEL}" parallel)
    if(FFTW3_${kind}_${parallel}_LIBRARY MATCHES "\\.a$")
      add_library(fftw3::${kind}::${parallel} STATIC IMPORTED GLOBAL)
    else()
      add_library(fftw3::${kind}::${parallel} SHARED IMPORTED GLOBAL)
    endif()

    # MPI Has a different included library than the others
    # FFTW3_INCLUDE_DIR_PARALLEL will change depending of which on is used.
    set(FFTW3_INCLUDE_DIR_PARALLEL ${FFTW3_INCLUDE_DIR} )
    if(PARALLEL STREQUAL "MPI")
      set(FFTW3_INCLUDE_DIR_PARALLEL ${FFTW3_${PARALLEL}_INCLUDE_DIR})
    endif()

    set_target_properties(fftw3::${kind}::${parallel} PROPERTIES
      IMPORTED_LOCATION "${FFTW3_${KIND}_${PARALLEL}_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR_PARALLEL}")

    # adding target properties to the different cases
    ##   MPI
    if(PARALLEL STREQUAL "MPI")
      if(MPI_C_LIBRARIES)
        set_target_properties(fftw3::${kind}::mpi PROPERTIES
          IMPORTED_LOCATION "${FFTW3_${KIND}_${PARALLEL}_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR_PARALLEL}"
          IMPORTED_LINK_INTERFACE_LIBRARIES ${MPI_C_LIBRARIES})
      endif()
    endif()
    ##   OpenMP
    if(PARALLEL STREQUAL "OPENMP")
      if(OPENMP_C_FLAGS)
        set_target_properties(fftw3::${kind}::${parallel} PROPERTIES
           IMPORTED_LOCATION "${FFTW3_${KIND}_${PARALLEL}_LIBRARY}"
           INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR_PARALLEL}"
           INTERFACE_COMPILE_OPTIONS "${OPENMP_C_FLAGS}")
        endif()
    endif()
    ##  THREADS
    if(PARALLEL STREQUAL "THREADS")
      if(CMAKE_THREAD_LIBS_INIT) # TODO: this is not running
        set_target_properties(fftw3::${kind}::${parallel} PROPERTIES
          IMPORTED_LOCATION "${FFTW3_${KIND}_${PARALLEL}_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR_PARALLEL}"
          INTERFACE_COMPILE_OPTIONS "${CMAKE_THREAD_LIBS_INIT}")
      endif()
    endif()
  endif()
endmacro()




if(NOT FFTW3_FIND_COMPONENTS)
  set(FFTW3_FIND_COMPONENTS SINGLE DOUBLE LONGDOUBLE SERIAL)
endif()

string(TOUPPER "${FFTW3_FIND_COMPONENTS}" FFTW3_FIND_COMPONENTS)

list(FIND FFTW3_FIND_COMPONENTS SINGLE LOOK_FOR_SINGLE)
list(FIND FFTW3_FIND_COMPONENTS DOUBLE LOOK_FOR_DOUBLE)
list(FIND FFTW3_FIND_COMPONENTS LONGDOUBLE LOOK_FOR_LONGDOUBLE)
list(FIND FFTW3_FIND_COMPONENTS THREADS LOOK_FOR_THREADS)
list(FIND FFTW3_FIND_COMPONENTS OPENMP LOOK_FOR_OPENMP)
list(FIND FFTW3_FIND_COMPONENTS MPI LOOK_FOR_MPI)
list(FIND FFTW3_FIND_COMPONENTS SERIAL LOOK_FOR_SERIAL)

# FIXME - This may fail in computers wihtout serial
# Default serial to obtain version number
set(LOOK_FOR_SERIAL 1)

# set serial as default if none parallel component has been set
if((LOOK_FOR_THREADS LESS 0) AND (LOOK_FOR_MPI LESS 0) AND
    (LOOK_FOR_OPENMP LESS 0))
  set(LOOK_FOR_SERIAL 1)
endif()

if(MPI_C_FOUND)
  set(MPI_FOUND ${MPI_C_FOUND})
endif()
unset(FFTW3_FIND_COMPONENTS)




if(WIN32)
  set(HINT_DIRS ${FFTW3_DIRECTORY} $ENV{FFTW3_DIRECTORY})
else()
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_FFTW QUIET fftw3)
    set(FFTW3_DEFINITIONS ${PC_FFTW3_CFLAGS_OTHER})
  endif()
  set(HINT_DIRS ${PC_FFTW3_INCLUDEDIR} ${PC_FFTW3_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIR} $ENV{FFTW3_INCLUDE_DIR} )
endif()

find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h HINTS ${HINT_DIRS})
if (LOOK_FOR_MPI)  # Probably is going to be the same as fftw3.h
  find_path(FFTW3_MPI_INCLUDE_DIR NAMES fftw3-mpi.h HINTS ${HINT_DIRS})
endif()

function(find_version OUTVAR LIBRARY SUFFIX)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/fftw${SUFFIX}/main.c
      # TODO: do we need to add include for mpi headers?
      "#include <fftw3.h>
       #include <stdio.h>
       int main(int nargs, char const *argv[]) {
           printf(\"%s\", fftw${SUFFIX}_version);
           return 0;
       }"
  )
if(NOT CMAKE_CROSSCOMPILING)
    try_run(RUN_RESULT COMPILE_RESULT
        "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/fftw${SUFFIX}/"
        "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/fftw${SUFFIX}/main.c"
        CMAKE_FLAGS
          -DLINK_LIBRARIES=${LIBRARY}
          -DINCLUDE_DIRECTORIES=${FFTW3_INCLUDE_DIR}
        RUN_OUTPUT_VARIABLE OUTPUT
        COMPILE_OUTPUT_VARIABLE COUTPUT
    )
  endif()
  if(RUN_RESULT EQUAL 0)
    string(REGEX REPLACE
        ".*([0-9]+\\.[0-9]+\\.[0-9]+).*"
        "\\1" VERSION_STRING "${OUTPUT}"
    )
    set(${OUTVAR} ${VERSION_STRING} PARENT_SCOPE)
  endif()
endfunction()

set(SUFFIX_DOUBLE "")
set(SUFFIX_SINGLE "f")
set(SUFFIX_LONGDOUBLE "l")
set(SUFFIX_SERIAL "")
set(SUFFIX_OPENMP "_omp")
set(SUFFIX_MPI "_mpi")
set(SUFFIX_THREADS "_threads")
set(SUFFIX_FINAL "")

if(WIN32)
  set(SUFFIX_FINAL "-3")
else()
  set(HINT_DIRS ${PC_FFTW3_LIBDIR} ${PC_FFTW3_LIBRARY_DIRS}
    $ENV{FFTW3_LIBRARY_DIR} ${FFTW3_LIBRARY_DIR} )
endif(WIN32)

unset(FFTW3_LIBRARIES)
set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR} ) # TODO what's for?
set(FFTW3_FLAGS_C "")
foreach(KIND SINGLE DOUBLE LONGDOUBLE)
  if(LOOK_FOR_${KIND} LESS 0)
    continue()
  endif()
  foreach(PARALLEL SERIAL MPI OPENMP THREADS)
    if(LOOK_FOR_${PARALLEL} LESS 0)
      continue()
    endif()
    find_specific_libraries(${KIND} ${PARALLEL})
  endforeach()
endforeach()

if(FFTW3_INCLUDE_DIR)
  list(GET FFTW3_FIND_COMPONENTS 0 smallerrun)
  string(REPLACE "_" ";" RUNLIST ${smallerrun})
  list(GET RUNLIST 0 KIND)
  list(GET RUNLIST 1 PARALLEL)
  unset(smallerrun)
  unset(RUNLIST)
  # suffix is quoted so it pass empty in the case of double as it's empty
  find_version(FFTW3_VERSION_STRING ${FFTW3_${KIND}_${PARALLEL}_LIBRARY}
    "${SUFFIX_${KIND}}")
endif()

# FIXME: fails if use REQUIRED.
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(FFTW3
    REQUIRED_VARS FFTW3_LIBRARIES FFTW3_INCLUDE_DIR
    VERSION_VAR FFTW3_VERSION_STRING
    HANDLE_COMPONENTS
)