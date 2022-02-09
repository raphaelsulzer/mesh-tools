set(GMP_PREFIX "" CACHE PATH "path ")


find_path(GMP_INCLUDE_DIR gmp.h gmpxx.h
    PATHS ${GMP_PREFIX}/include /usr/include /usr/local/include )

find_library(GMP_LIB NAMES gmp libgmp
    PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)

find_library(GMPXX_LIB NAMES gmpxx libgmpxx
    PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)



if(GMP_INCLUDE_DIR AND GMP_LIB AND GMPXX_LIB)
    #get_filename_component(GMP_LIBRARY_DIR ${GMP_LIBRARY} PATH)
    set(GMP_FOUND TRUE)
    set(GMP_LIBRARY ${GMP_LIB} ${GMPXX_LIB})
endif()

if(GMP_FOUND)
   if(NOT GMP_FIND_QUIETLY)
      MESSAGE(STATUS "Found GMP: ${GMP_LIBRARY}")
   endif()
elseif(GMP_FOUND)
   if(GMP_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find GMP")
   endif()
endif()


## from here: https://github.com/seahorn/crab/blob/master/cmake/FindGMP.cmake

# GMP_USE_STATIC_LIBS   - Set to ON to force the use of static libraries.

## To switch between static and dynamic without resetting the cache
# if ("${GMP_USE_STATIC_LIBS}" STREQUAL "${GMP_USE_STATIC_LIBS_LAST}")
#   set(GMP_USE_STATIC_LIBS_CHANGED OFF)
# else ()
#   set(GMP_USE_STATIC_LIBS_CHANGED ON)
# endif()
# set(GMP_USE_STATIC_LIBS_LAST "${GMP_USE_STATIC_LIBS}")
# if (GMP_USE_STATIC_LIBS_CHANGED)
#   set(GMP_LIB "GMP_LIB-NOTFOUND")
# endif()

## Support preference of static libs by adjusting CMAKE_FIND_LIBRARY_SUFFIXES
#if(GMP_USE_STATIC_LIBS )
#  # save CMAKE_FIND_LIBRARY_SUFFIXES
#  set(_GMP_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
#  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
#endif()

#set(GMP_SEARCH_PATH "" CACHE PATH "Search path for gmp.")
#find_path(GMP_INCLUDE_DIR NAMES gmp.h PATHS ${GMP_SEARCH_PATH}/include)
#find_library(GMP_LIB NAMES gmp gmpxx PATHS ${GMP_SEARCH_PATH}/lib)

#mark_as_advanced(GMP_SEARCH_PATH GMP_INCLUDE_DIR GMP_LIB)

#include (FindPackageHandleStandardArgs)
#find_package_handle_standard_args(GMP
#  REQUIRED_VARS GMP_INCLUDE_DIR GMP_LIB)

#if(GMP_USE_STATIC_LIBS )
#  # restore CMAKE_FIND_LIBRARY_SUFFIXES
#  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_GMP_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
#endif()


#if(GMP_FOUND)
#  set(GMP_LIBRARY ${GMP_LIB} ${GMPXX_LIB})
#  set(GMP_LIBRARIES ${GMP_LIBRARY})
#   if(NOT GMP_FIND_QUIETLY)
#      MESSAGE(STATUS "Found GMP: ${GMP_LIBRARY}")
#   endif()
#elseif(GMP_FOUND)
#   if(GMP_FIND_REQUIRED)
#      message(FATAL_ERROR "Could not find GMP")
#   endif()
#endif()
