# Findtclap.cmake
#
# A simple Find module for `find_package(tclap)`, for the TCLAP project:
#
#     http://tclap.sourceforge.net/
#
# This Find module respects ``${tclap_ROOT}`` and ``$ENV{tclap_ROOT}`` variables.  See
# CMake Policy CMP0074 for more information.
#
# The following variables are set:
#
# ``tclap_FOUND``
#     ``TRUE`` if package was found, ``FALSE`` if not found.
#
# ``tclap_INCLUDE_DIR``
#     The include directory if needed.  Only set if ``tclap_FOUND``.
#
# If ``tclap_FOUND``, then the imported interface library ``tclap::tclap`` is created
# with the appropriate include directory.  Example usage:
#
# .. code-block:: cmake
#
#     find_package(tclap REQUIRED)
#     target_link_libraries(mylib PUBLIC tclap::tclap)
#     # ... or ...
#     target_include_directories(mylib PUBLIC ${tclap_INCLUDE_DIR})
#
# .. warning::
#
#     At least at this time, TCLAP only includes version information in the installed
#     ``pkg-config`` file ``tclap.pc``.  This means that version information cannot
#     be extracted when ``pkg-config`` is not available.  Since TCLAP has a very stable
#     interface, in the event that version information cannot be extracted then
#
#     1. ``find_package(tclap X.Y.Z)`` will warn that no version information was found.
#     2. ``find_package(tclap X.Y.Z EXACT)`` will fail, since no version information was
#        recovered.
#
#     Unless you truly need to force exact version information, it is encouraged to
#     *avoid* using ``EXACT`` in ``find_package`` calls for ``tclap``.
if (POLICY CMP0074)
  cmake_policy(PUSH)
  cmake_policy(SET CMP0074 NEW)
endif()

# tclap uses Autotools and by default will install `tclap.pc`.  Try and use pkg-config
# to find the package details (potentially overidden by `tclap_ROOT` below).
find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(pc_tclap QUIET tclap)
  if (pc_tclap_VERSION)
    set(tclap_VERSION ${pc_tclap_VERSION})
  endif()
else()
  set(tclap_VERSION "tclap_VERSION-NOTFOUND")
endif()

# NOTE: order matters, search the `tclap_ROOT` variables first.  That is, `tclap_ROOT`
# takes precedence over any potential results from `pkg-config`.
find_path(
  tclap_INCLUDE_DIR /usr/local/tclap-1.2.2/include/tclap/CmdLine.h
  HINTS ${tclap_ROOT} $ENV{tclap_ROOT} ${pc_tclap_INCLUDE_DIRS}
)

# If tclap/CmdLine.h is not found, `if (tclap_INCLUDE_DIR)` evaluates to FALSE.
if (tclap_INCLUDE_DIR)
  set(tclap_FOUND TRUE)
  add_library(tclap::tclap IMPORTED INTERFACE)
  target_include_directories(tclap::tclap INTERFACE ${tclap_INCLUDE_DIR})
else()
  set(tclap_FOUND FALSE)
  if (tclap_FIND_REQUIRED)
    message(FATAL_ERROR
      "tclap could not be found!  Try setting the environment variable `tclap_ROOT` "
      "such that the file `$tclap_ROOT/tclap/CmdLine.h` can be found."
    )
  endif()
endif()

# Strategy: the tclap API is very stable / does not change often.  Since it is a header
# only library, being able to wield `tclap_ROOT` is very convenient (e.g., Windows users
# who do not want to setup `pkg-config`).  However, tclap only includes version
# information in the `tclap.pc` file.  As such, if tclap_VERSION is not able to be found
#
# - ``find_package(tclap X.Y.Z)`` => warn.
# - ``find_package(tclap X.Y.Z EXACT)`` => fail.
#
# In practice, users should only be doing ``find_package(tclap [REQUIRED] [QUIET])``,
# since enforcing a version constraint is likely not necessary.
if (tclap_FIND_VERSION)
  if (NOT tclap_VERSION)
    if (tclap_FIND_VERSION_EXACT)
      message(FATAL_ERROR
        "EXACT tclap version of ${tclap_FIND_VERSION} requested, but `tclap_VERSION` "
        "was unable to be extracted."
      )
    else()
      message(WARNING
        "tclap_VERSION was unable to be extracted, so it is not guaranteed to match "
        "specified tclap version of ${tclap_FIND_VERSION}"
      )
    endif()
  endif()
endif()

if (POLICY CMP0074)
  cmake_policy(POP)
endif()