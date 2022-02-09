# Find package module for CNPY library.
#
# The following variables are set by this module:
#
#   CNPY_FOUND: TRUE if CNPY is found.
#   CNPY_INCLUDE_DIRS: Include directories for CNPY.
#   CNPY_LIBRARIES: Libraries required to link CNPY.
#
# The following variables control the behavior of this module:
#
# CNPY_INCLUDE_DIR_HINTS: List of additional directories in which to
#                         search for CNPY includes.
# CNPY_LIBRARY_DIR_HINTS: List of additional directories in which to
#                         search for CNPY libraries.

set(CNPY_INCLUDE_DIR_HINTS "" CACHE PATH "CNPY include directory")
set(CNPY_LIBRARY_DIR_HINTS "" CACHE PATH "CNPY library directory")

unset(CNPY_FOUND)
unset(CNPY_ROOT)
unset(CNPY_LIBRARY)

find_path(CNPY_ROOT
    NAMES
    GL/CNPY.h
    cnpy.h
    PATHS
    ${CNPY_INCLUDE_DIR_HINTS}
    /usr/include
    /usr/local/include
    /sw/include
    /opt/include
    /opt/local/include
    /usr/local
    /usr/lib
    CNPY
    CNPY
    )
find_library(CNPY_LIBRARY
    NAMES
    cnpy
    PATHS
    ${CNPY_LIBRARY_DIR_HINTS}
    /usr/lib64
    /usr/lib
    /usr/local/lib64
    /usr/local/lib
    /sw/lib
    /opt/lib
    /opt/local/lib
    /usr/local
    )


if(CNPY_ROOT AND CNPY_LIBRARY)
    set(CNPY_FOUND TRUE)
    message(STATUS "Found CNPY")
    message(STATUS "  Includes : ${CNPY_ROOT}")
    message(STATUS "  Libraries : ${CNPY_LIBRARY}")
else()
    set(CNPY_FOUND FALSE)
    if(CNPY_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find CNPY")
    endif()
endif()
