# Find package module for TIFF library.
#
# The following variables are set by this module:
#
#   TIFF_FOUND: TRUE if TIFF is found.
#   TIFF_INCLUDE_DIRS: Include directories for TIFF.
#   TIFF_LIBRARIES: Libraries required to link TIFF.
#
# The following variables control the behavior of this module:
#
# TIFF_INCLUDE_DIR_HINTS: List of additional directories in which to
#                         search for TIFF includes.
# TIFF_LIBRARY_DIR_HINTS: List of additional directories in which to
#                         search for TIFF libraries.

set(TIFF_INCLUDE_DIR_HINTS "" CACHE PATH "TIFF include directory")
set(TIFF_LIBRARY_DIR_HINTS "" CACHE PATH "TIFF library directory")

unset(TIFF_FOUND)
unset(TIFF_INCLUDE_DIRS)
unset(TIFF_LIBRARIES)

find_path(TIFF_INCLUDE_DIRS
    NAMES
    GL/TIFF.h
    tiff.h
    PATHS
    ${TIFF_INCLUDE_DIR_HINTS}
    /usr/include
    /usr/local/include
    /sw/include
    /opt/include
    /opt/local/include
    /usr/local/Cellar)
find_library(TIFF_LIBRARIES
    NAMES
    tiff
    TIFF
    TIFF
    TIFF32
    PATHS
    ${TIFF_LIBRARY_DIR_HINTS}
    /usr/lib64
    /usr/lib
    /usr/local/lib64
    /usr/local/lib
    /sw/lib
    /opt/lib
    /opt/local/lib
    /usr/local/Cellar)


if(TIFF_INCLUDE_DIRS AND TIFF_LIBRARIES)
    set(TIFF_FOUND TRUE)
    message(STATUS "Found TIFF")
    message(STATUS "  Includes : ${TIFF_INCLUDE_DIRS}")
    message(STATUS "  Libraries : ${TIFF_LIBRARIES}")
else()
    set(TIFF_FOUND FALSE)
    if(TIFF_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find TIFF")
    endif()
endif()
