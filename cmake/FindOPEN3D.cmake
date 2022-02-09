# Find package module for OPEN3D library.
#
# The following variables are set by this module:
#
#   OPEN3D_FOUND: TRUE if OPEN3D is found.
#   OPEN3D_INCLUDE_DIRS: Include directories for OPEN3D.
#   OPEN3D_LIBRARIES: Libraries required to link OPEN3D.
#
# The following variables control the behavior of this module:
#
# OPEN3D_INCLUDE_DIR_HINTS: List of additional directories in which to
#                         search for OPEN3D includes.
# OPEN3D_LIBRARY_DIR_HINTS: List of additional directories in which to
#                         search for OPEN3D libraries.

set(OPEN3D_INCLUDE_DIR_HINTS "" CACHE PATH "OPEN3D include directory")
set(OPEN3D_LIBRARY_DIR_HINTS "" CACHE PATH "OPEN3D library directory")

unset(OPEN3D_FOUND)

find_path(OPEN3D_INCLUDE_DIRS
    NAMES
    GL/OPEN3D.h
    OPEN3D.h
    OPEN3Dptimization.h
    PATHS
    ${OPEN3D_INCLUDE_DIR_HINTS}
    /usr/include
    /usr/local/include
    /sw/include
    /opt/include
    /opt/local/include/Open3D
    /usr/local/include/Open3D
    /usr/local
    /usr/lib)
find_library(OPEN3D_LIBRARIES
    NAMES
    OPEN3D
    OPEN3D
    OPEN3D32
    PATHS
    ${OPEN3D_LIBRARY_DIR_HINTS}
    /usr/lib64
    /usr/lib
    /usr/local/lib64
    /usr/local/lib
    /sw/lib
    /opt/lib
    /opt/local/lib
    /usr/local
    /usr/lib/OPEN3D-v3.0-master/build
    /usr/local/OPEN3D-v3.0-master/build
    /usr/lib/OPEN3D-v3.0-master
    /usr/local/OPEN3D-v3.0-master)

if(OPEN3D_INCLUDE_DIRS AND OPEN3D_LIBRARIES)
    set(OPEN3D_FOUND TRUE)
endif()
