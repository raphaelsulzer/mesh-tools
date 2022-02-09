# Try to find the MPFR library
# See http://www.mpfr.org/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(MPFR 2.3.0)
# to require version 2.3.0 to newer of MPFR.
#
# Once done this will define
#
#  MPFR_FOUND - system has MPFR lib with correct version
#  MPFR_INCLUDES - the MPFR include directory
#  MPFR_LIBRARIES - the MPFR library
#  MPFR_VERSION - MPFR version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2010 Jitse Niesen, <jitse@maths.leeds.ac.uk>
# Copyright (c) 2015 Jack Poulson, <jack.poulson@gmail.com>
# Redistribution and use is allowed according to the terms of the BSD license.

find_path(MPFR_INCLUDES NAMES mpfr.h PATHS $ENV{GMPDIR} $ENV{MPFRDIR}
  ${INCLUDE_INSTALL_DIR})

# Set MPFR_FIND_VERSION to 1.0.0 if no minimum version is specified
if(NOT MPFR_FIND_VERSION)
  if(NOT MPFR_FIND_VERSION_MAJOR)
    set(MPFR_FIND_VERSION_MAJOR 1)
  endif()
  if(NOT MPFR_FIND_VERSION_MINOR)
    set(MPFR_FIND_VERSION_MINOR 0)
  endif()
  if(NOT MPFR_FIND_VERSION_PATCH)
    set(MPFR_FIND_VERSION_PATCH 0)
  endif()
  set(MPFR_FIND_VERSION
    "${MPFR_FIND_VERSION_MAJOR}.${MPFR_FIND_VERSION_MINOR}.${MPFR_FIND_VERSION_PATCH}")
endif()

if(MPFR_INCLUDES)
  # Query MPFR_VERSION
  file(READ "${MPFR_INCLUDES}/mpfr.h" _mpfr_version_header)

  string(REGEX MATCH "define[ \t]+MPFR_VERSION_MAJOR[ \t]+([0-9]+)"
    _mpfr_major_version_match "${_mpfr_version_header}")
  set(MPFR_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+MPFR_VERSION_MINOR[ \t]+([0-9]+)"
    _mpfr_minor_version_match "${_mpfr_version_header}")
  set(MPFR_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+MPFR_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
    _mpfr_patchlevel_version_match "${_mpfr_version_header}")
  set(MPFR_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

  set(MPFR_VERSION
    ${MPFR_MAJOR_VERSION}.${MPFR_MINOR_VERSION}.${MPFR_PATCHLEVEL_VERSION})

  # Check whether found version exceeds minimum required
  if(${MPFR_VERSION} VERSION_LESS ${MPFR_FIND_VERSION})
    set(MPFR_VERSION_OK FALSE)
    message(STATUS "MPFR version ${MPFR_VERSION} found in ${MPFR_INCLUDES}, "
                   "but at least version ${MPFR_FIND_VERSION} is required")
  else()
    set(MPFR_VERSION_OK TRUE)
  endif()
endif()

find_library(MPFR_LIBRARIES mpfr
  PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG
                                  MPFR_INCLUDES MPFR_LIBRARIES MPFR_VERSION_OK)
mark_as_advanced(MPFR_INCLUDES MPFR_LIBRARIES)

if(MPFR_FOUND)
   if(NOT MPFR_FIND_QUIETLY)
      MESSAGE(STATUS "Found MPFR: ${MPFR_LIBRARIES}")
   endif()
elseif(MPFR_FOUND)
   if(MPFR_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find MPFR")
   endif()
endif()


## Find mpfr with API version ?.?
##
## Usage:
##   find_package(mpfr [REQUIRED] [QUIET])
##
## Sets the following variables:
##   - MPFR_FOUND          .. true if library is found
##   - MPFR_LIBRARY      .. full path to library
##   - MPFR_INCLUDE_DIR    .. full path to include directory
##
## Honors the following optional variables:
##   - MPFR_INCLUDE_LOC    .. include directory path, to be searched before defaults
##   - MPFR_LIBRARY_LOC    .. the library's directory path, to be searched before defaults
##   - MPFR_STATIC_LIBRARY .. if true, find the static library version
##
## Copyright 2015 Joachim Coenen, Forschungszentrum Jülich.
## Redistribution permitted.

## find the mpfr include directory
#find_path(MPFR_INCLUDE_DIR mpfr.h
#    PATH_SUFFIXES include mpfr/include mpfr
#    PATHS
#    ${MPFR_INCLUDE_LOC}
#    "C:/Program Files/mpfr/"
#    ~/Library/Frameworks/
#    /Library/Frameworks/
#    /usr/local/
#    /usr/
#    /sw/ # Fink
#    /opt/local/ # DarwinPorts
#    /opt/csw/ # Blastwave
#    /opt/
#    )

#set(CMAKE_REQUIRED_INCLUDES ${MPFR_INCLUDE_DIR})
#set(CMAKE_REQUIRED_QUIET False)

## attempt to find static library first if this is set
#if(MPFR_STATIC_LIBRARY)
#    set(MPFR_STATIC mpfr.a)
#    #set(MPFR_STATIC mpfr.lib)
#endif()

## find the mpfr library
#find_library(MPFR_LIBRARY_RELEASE
#    NAMES
#    #${MPFR_STATIC}
#    mpfr
#    PATH_SUFFIXES lib64 lib
#    PATHS
#    ${MPFR_LIBRARY_LOC}
#    "C:/Program Files/mpfr/Release"
#    ~/Library/Frameworks
#    /Library/Frameworks
#    /usr/local
#    /usr
#    /sw
#    /opt/local
#    /opt/csw
#    /opt
#    )
#if (MPFR_USE_DEBUG_LIBRARY)
#    # find the mpfr library
#    find_library(MPFR_LIBRARY_DEBUG
#        NAMES
#        #${MPFR_STATIC}
#        mpfr_d
#        PATH_SUFFIXES lib64 lib
#        PATHS
#        ${MPFR_LIBRARY_LOC}
#        "C:/Program Files/mpfr/Debug"
#        ~/Library/Frameworks
#        /Library/Frameworks
#        /usr/local
#        /usr
#        /sw
#        /opt/local
#        /opt/csw
#        /opt
#        )
#else()
#    set(MPFR_LIBRARY_DEBUG ${MPFR_LIBRARY_RELEASE})
#endif()

##message(STATUS "Found FindMPFR (release): ${MPFR_LIBRARY_RELEASE}")
##message(STATUS "Found FindMPFR (release): ${MPFR_LIBRARY_RELEASE}")
#set(MPFR_LIBRARY ${MPFR_LIBRARY_RELEASE})
#set(MPFR_LIBRARIES ${MPFR_LIBRARY_RELEASE})
#message(STATUS "Found MPFR: \n\t-include: ${MPFR_INCLUDE_DIR}\n\t-lib: ${MPFR_LIBRARY}")
#mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY_RELEASE)
#mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY_DEBUG)
