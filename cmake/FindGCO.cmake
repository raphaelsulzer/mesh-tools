# Find package module for GCO library.
#
# The following variables are set by this module:
#
#   GCO_FOUND: TRUE if GCO is found.
#   GCO_INCLUDE_DIRS: Include directories for GCO.
#   GCO_LIBRARIES: Libraries required to link GCO.
#
# The following variables control the behavior of this module:
#
# GCO_INCLUDE_DIR_HINTS: List of additional directories in which to
#                         search for GCO includes.
# GCO_LIBRARY_DIR_HINTS: List of additional directories in which to
#                         search for GCO libraries.

set(GCO_INCLUDE_DIR_HINTS "" CACHE PATH "GCO include directory")
set(GCO_LIBRARY_DIR_HINTS "" CACHE PATH "GCO library directory")

unset(GCO_FOUND)
unset(GCO_ROOT)
unset(GCO_LIBRARY)

find_path(GCO_ROOT
    NAMES
    GL/GCO.h
    gco.h
    GCoptimization.h
    PATHS
    ${GCO_INCLUDE_DIR_HINTS}
    /usr/include
    /usr/local/include
    /sw/include
    /opt/include
    /opt/local/include
    /usr/local
    /usr/lib
    /usr/lib/gco-v3.0-master/build
    /usr/local/include/gco-v3.0-master/build
    /usr/local/include/gco-v3.0/build
    /usr/local/include/gco-v3.0
    /usr/local/include/gco-v3.0-master
    /usr/local/gco-v3.0-master/build
    /usr/lib/gco-v3.0-master
    /usr/local/gco-v3.0-master
    /usr/local/gco-v3.0
    )
find_library(GCO_LIBRARY
    NAMES
    gco
    GCO
    GCO32
    PATHS
    ${GCO_LIBRARY_DIR_HINTS}
    /usr/lib64
    /usr/lib
    /usr/local/lib64
    /usr/local/lib
    /sw/lib
    /opt/lib
    /opt/local/lib
    /usr/local
    /usr/lib/gco-v3.0-master/build
    /usr/local/gco-v3.0-master/build
    /usr/local/include/gco-v3.0-master/build
    /usr/local/include/gco-v3.0-master
    /usr/lib/gco-v3.0-master
    /usr/local/gco-v3.0-master
    /usr/lib/gco-v3.0/build
    /usr/local/gco-v3.0/build
    /usr/local/include/gco-v3.0/build
    /usr/local/include/gco-v3.0
    /usr/lib/gco-v3.0
    /usr/local/gco-v3.0
    )


if(GCO_ROOT AND GCO_LIBRARY)
    set(GCO_FOUND TRUE)
    message(STATUS "Found GCO")
    message(STATUS "  Includes : ${GCO_ROOT}")
    message(STATUS "  Libraries : ${GCO_LIBRARY}")
else()
    set(GCO_FOUND FALSE)
    if(GCO_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find GCO")
    endif()
endif()
