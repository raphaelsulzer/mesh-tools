cmake_minimum_required(VERSION 3.10)
project(mesh-tools)

################################################################################
# Include CMake dependencies
################################################################################

# CMAKE_MODULE_PATH needs to be set for cmake to find the other CMakeLists.txt in ./cmake/
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
message(STATUS "cmake location: ${CMAKE_INSTALL_PREFIX}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# print which compiler is used
message(STATUS "used compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "used compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

################################################################################
# Build Options
################################################################################
option(OPENMP_ENABLED "Whether to enable OpenMP support" ON)
option(O3D_ENABLED "Whether to enable OPEN3D" OFF)
option(PCL_ENABLED "Whether to enable PCL" OFF)
option(OMVS_ENABLED "Whether to enable OpenMVS, also enables OpenCV" OFF)
option(COLMAP_ENABLED "Wether to enable colmap support" OFF)
option(RECONBENCH_ENABLED "Wether to enable reconbench support" OFF)
option(STATIC "Whether to compile for standalone use" OFF)

################################################################################
# Programs to build
################################################################################
option(LABATUT "Whether to build labatut exe" ON)
option(OCC2MESH "Whether to build occ2mesh exe" OFF)
option(SCAN "Whether to build scan exe" ON)
option(FEAT "Whether to build feat exe" ON)
option(OMVS2NPZ "Whether to build omvs2npz exe" OFF)


if(STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

if(CMAKE_BUILD_TYPE)
    message(STATUS "Build type specified as ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "Build type not specified, using Debug")
    set(CMAKE_BUILD_TYPE Debug)
    set(IS_DEBUG ON)
endif()

### surpress some warnings when make
### taken from the cmake file here: http://www.open3d.org/docs/release/tutorial/C++/cplusplus_interface.html
if(UNIX)
    add_definitions(-DUNIX)
    add_compile_options(-Wno-deprecated-declarations)
    add_compile_options(-Wno-unused-result)
    add_compile_options(-Wno-all -Wno-attributes -Wno-enum-compare -Wno-delete-incomplete)
    add_definitions(-O3)
endif(UNIX)
### to surpress ROOT variable is set warning, taken from: https://gitlab.kitware.com/vtk/vtk/commit/e6f4ab044acba1363c16eda552ecb691557cdfd5
foreach(policy
    CMP0074 # CMake 3.12
    )
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

################################################################################
# Find packages
################################################################################

if(COLMAP_ENABLED)
    find_package(Glew REQUIRED)
endif()
if(STATIC)
    set(GMP_USE_STATIC_LIBS ON)
    set(MPFR_STATIC_LIBRARY TRUE)
    set(OpenCV_STATIC ON)
    set(OpenCV_SHARED OFF)
endif()
#find_package(GMP REQUIRED) # is automatically linked by cgal
#find_package(MPFR REQUIRED)
find_package(GCO REQUIRED)
find_package(xtl REQUIRED) # for xtensor
find_package(ZLIB REQUIRED) # for xtensor
find_package(xtensor REQUIRED)
find_package(xtensor-io REQUIRED)


####################################
############### EIGEN ##############
####################################
FIND_PACKAGE(Eigen REQUIRED)
#if(EIGEN_FOUND)
#        INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
#        ADD_DEFINITIONS(${EIGEN_DEFINITIONS} -D_USE_EIGEN)
#        SET(_USE_EIGEN TRUE)
#endif()
if(IS_GNU)
    # Hide incorrect warnings for uninitialized Eigen variables under GCC.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-maybe-uninitialized")
    #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized")
    list(APPEND CMAKE_CXX_FLAGS "-Wno-maybe-uninitialized")
endif()
if(IS_DEBUG)
    add_definitions("-DEIGEN_INITIALIZE_MATRICES_BY_NAN")
endif()

####################################
################ GLOG ##############
####################################
#find_package(Glog REQUIRED)
if(IS_MSVC)
    # Some fixes for the Glog library.
    add_definitions("-DGLOG_NO_ABBREVIATED_SEVERITIES")
    add_definitions("-DGL_GLEXT_PROTOTYPES")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
    list(APPEND CMAKE_CXX_FLAGS "/EHsc")
#    # Enable object level parallel builds in Visual Studio.
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

####################################
############ RECONBENCH ############
####################################
if(RECONBENCH_ENABLED)
    # what I am actually doing here is linking the two libraries OpenMesh and ANN that are needed for the reconbench code
    # instead of linking the versions provided by reconbench I should install them on the machine and link to that
    # but probably then I should also do the same for the reconbench code itself.
    # question is just if the recent versions of these libs are compatible with the reconbench code
    add_definitions(-DRECONBENCH)

    # LAPACK (global), OpenMesh and ANN is needed for reconbench
    find_package(LAPACK)
    list(APPEND RECONBENCH_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})
    list(APPEND RECONBENCH_INCLUDE_DIRS /home/rsulzer/cpp/reconbench-CMake/OpenMesh)
    list(APPEND RECONBENCH_EXTERNAL_LIBRARIES /home/rsulzer/cpp/reconbench-CMake/OpenMesh/build/libOpenMesh.a)
    list(APPEND RECONBENCH_INCLUDE_DIRS /home/rsulzer/cpp/reconbench-CMake/ann/include/ANN)
    list(APPEND RECONBENCH_EXTERNAL_LIBRARIES /home/rsulzer/cpp/reconbench-CMake/ann/build/libann.a)

    # TODO: make the includes for this as a releative path, like this:
    #    set(RECONBENCH_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../reconbench-CMake/)

    # this just collects all the cpp files in the modeling directory
    file(GLOB_RECURSE RECONBENCH_CPP_FILES external/modeling/*.cpp)
endif()

####################################
################ PCL ###############
####################################
if(PCL_ENABLED)
    add_definitions(-DPCL)
    find_package(PCL 1.7 REQUIRED)
    if(PCL_FOUND)
        list(APPEND PCL_INCLUDE_DIRS ${PCL_INCLUDE_DIRS})
        list(APPEND PCL_EXTERNAL_LIBRARIES ${PCL_LIBRARIES})
        link_directories(${PCL_LIBRARY_DIRS})
#        list(APPEND PCL_CPP_FILES
#            src/IO/ethIO.cpp
#            src/IO/tinyxml2.cpp
#            src/IO/meshlab_project.cc)
        add_definitions(${PCL_DEFINITIONS})
    endif()
endif()

####################################
############## OPENMVS #############
####################################
if(OMVS_ENABLED)
    add_definitions(-DOpenMVS)
    #find_package(OpenMVS REQUIRED)
    # it is important that cmake finds this file:
    # ${CMAKE_INSTALL_PREFIX}/lib/CMake/OpenMVS/OpenMVSConfig.cmake, which will only be there if you do make install
    # after compiling OpenMVS. Afterwards the ${OpenMVS_INCLUDE_DIRS}" should point to the installed headers in
    # /usr/local/lib/CMake/OpenMVS/../../../include/OpenMVS and not to your github OpenMVS folder, because then the paths are not
    # correct. Also I found that you have to build OpenMVS with shared libraries activated.
    # UPDATE: building OpenMVS with -DSHARED_LIBS_ON=TRUE does not work anymore. Problem is OpenCV, could try to compile OpenCV from source
    # or try with OpenMVS 2.0
    find_package(OpenMVS HINTS ${CMAKE_INSTALL_PREFIX}/lib/CMake REQUIRED)
    if(OpenMVS_FOUND)
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${cxx_default}")
        # warning: OpenMVS does not work with CXX_FLAG -D_GLIBCXX_USE_CXX11_ABI=0
        # which e.g. the default setting for Open3D compilation includes
        list(APPEND CMAKE_CXX_FLAGS ${cxx_default})
        ADD_DEFINITIONS(${OpenMVS_DEFINITIONS})
        message(STATUS "Found OpenMVS ${OpenMVS_VERSION}")
        message(STATUS "  Includes : ${OpenMVS_INCLUDE_DIRS}")
        message(STATUS "  Libraries : ${OpenMVS_LIBRARIES}")
        list(APPEND OMVS_INCLUDE_DIRS ${OpenMVS_INCLUDE_DIRS})
        list(APPEND OMVS_EXTERNAL_LIBRARIES ${OpenMVS_LIBRARIES} dl MVS)
        list(APPEND OMVS_CPP_FILES src/IO/mvIO.cpp)
    else()
        message(SEND_ERROR "OpenMVS not found")
    endif()
endif()


####################################
############## COLMAP ##############
####################################
set(COLMAP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../colmap/build)
if(COLMAP_ENABLED)
    set(COLMAP_DIR /usr/local/share/colmap)
    find_package(COLMAP REQUIRED)
    add_definitions(-DCOLMAP)
    if(COLMAP_FOUND)
        message(STATUS "Found COLMAP")
        message(STATUS "Be aware that colmap works best with eigen version 3.3.3!!!")
        list(APPEND COLMAP_INCLUDE_DIRS ${COLMAP_INCLUDE_DIRS})
        list(APPEND COLMAP_EXTERNAL_LIBRARIES ${COLMAP_LIBRARIES})
        list(APPEND CPP_FILES src/IO/colmapIO.cpp)
        link_directories(${COLMAP_LINK_DIRS})
    else()
        message(SEND_ERROR "COLMAP enabled but not found")
    endif()
endif()

####################################
############## OPEN3D ##############
####################################
if(O3D_ENABLED)
    message(STATUS "open3d is enabled")
    add_definitions(-DOpen3D)
    find_package(Open3D HINTS ${CMAKE_INSTALL_PREFIX}/lib/CMake)
    if(Open3D_FOUND)


        # IMPORTANT: OPEN3D NEEDS TO BE BUILD WITH THIS CMAKE COMMAND
        # cmake -DBUILD_EIGEN3=ON -DBUILD_GLEW=ON -DBUILD_GLFW=ON -DBUILD_JSONCPP=OFF -DBUILD_PNG=OFF -DGLIBCXX_USE_CXX11_ABI=ON -DBUILD_UNIT_TESTS=ON ..
        # on Ubuntu 16.04 also this needs to be added BUILD_FILAMENT_FROM_SOURCE=ON
        # the important thing is probably -DGLIBCXX_USE_CXX11_ABI=ON
        # see here: https://github.com/intel-isl/Open3D/issues/2286


#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Open3D_C_FLAGS}")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Open3D_CXX_FLAGS}")
#        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Open3D_EXE_LINKER_FLAGS}")

#        set(GLIBCXX_USE_CXX11_ABI ON)

        list(APPEND Open3D_LIBRARIES dl)
        message(STATUS "Found Open3D ${Open3D_VERSION}")
        message(STATUS "  Includes : ${Open3D_INCLUDE_DIRS}")
        message(STATUS "  Libraries : ${Open3D_LIBRARIES}")
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Open3D_C_FLAGS}")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Open3D_CXX_FLAGS}")


#        list(APPEND CMAKE_CXX_FLAGS ${Open3D_CXX_FLAGS})
#        list(APPEND CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0 -fPIC")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC")
#        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${Open3D_EXE_LINKER_FLAGS}")
#        list(APPEND O3D_INCLUDE_DIRS ${Open3D_INCLUDE_DIRS})
#        list(APPEND O3D_EXTERNAL_LIBRARIES ${Open3D_LIBRARIES})
        link_directories(${Open3D_LIBRARY_DIRS})
    else()
        message(SEND_ERROR "Open3D not found")
    endif()
endif()

####################################
############## OPENMP ##############
####################################
if(OPENMP_ENABLED)
    # see here how to link it since cmake 3.9: https://cliutils.gitlab.io/modern-cmake/chapters/packages/OpenMP.html
    find_package(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        message(STATUS "Enabling OpenMP support")
##        target_link_libraries(sure PUBLIC OpenMP::OpenMP_CXX)
#        add_definitions("-DOPENMP_ENABLED")
##        list(APPEND CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
###        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
##        list(APPEND CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
##        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    else()
        message(STATUS "OpenMP enabled but not found!")
    endif()
else()
    if(OMVS_ENABLED)
        message(STATUS "CAREFUL! OPENMP IS PROBABLY NEEDED FOR OPEN_MVS!!")
    endif()
    message(STATUS "Disabling OpenMP support")
endif()

####################################
############### CNPY ###############
####################################
#find_package(CNPY REQUIRED)
## not really necessary, was just for maybe replacing xtensor to load NPY files.

####################################
############### CGAL ###############
####################################

set(CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG 1)
set(CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT 1)


if(STATIC)
    set(CGAL_Boost_USE_STATIC_LIBS   ON)
endif()
# to use local CGAL one has to set Xpackage_DIR variable to the path where the XpackageConfig.cmake file is located
set(CGAL_DIR /home/adminlocal/PhD/cpp/CGAL/lib/cmake/CGAL)
find_package(CGAL 5.1.1 REQUIRED CONFIG)
if(CGAL_FOUND)
    set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
    message(STATUS "Found CGAL ${CGAL_VERSION}")
    message(STATUS "  Includes : ${CGAL_INCLUDE_DIRS}")
    message(STATUS "  Libraries : ${CGAL_LIBRARY}")
else()
    message(STATUS "Disabling CGAL support")
endif()
#set(CGAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../CGAL/)
####################################
############## BOOST ###############
####################################
if(STATIC)
    set(Boost_USE_STATIC_LIBS   ON)
endif()
# adding boost explicitly after CGAL because CGAL also calls find_package(boost)
# but without the required components, and thus sets Boost_Libraries to empty.
FIND_PACKAGE(Boost REQUIRED COMPONENTS iostreams program_options system serialization filesystem)
#if(Boost_FOUND)
#        ADD_DEFINITIONS(${Boost_DEFINITIONS} -D_USE_BOOST)
#        SET(_USE_BOOST TRUE)
#endif()

####################################
############## OPENCV ##############
####################################
if(OMVS_ENABLED)
    FIND_PACKAGE(OpenCV REQUIRED)
    if(OpenCV_FOUND)
    #            INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
            ADD_DEFINITIONS(${OpenCV_DEFINITIONS})
            SET(_USE_OPENCV TRUE)
            MESSAGE(STATUS "OpenCV ${OpenCV_VERSION} found (include: ${OpenCV_INCLUDE_DIRS})")
            list(APPEND OCV_EXTERNAL_LIBRARIES ${OpenCV_LIBS})
            message("ocv libs ${OpenCV_LIBS}")
    else()
            MESSAGE("-- Can't find OpenCV. Please specify OpenCV directory using OpenCV_DIR variable")
    endif()
endif()



################################################################################
############################### BUILD PROJECT ##################################
################################################################################

include_directories(
    src
    external
    ${EIGEN_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${GLOG_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
#    ${GLEW_INCLUDE_DIRS}
    ${GMP_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
    ${GCO_ROOT}
    ${CGAL_INCLUDE_DIRS}
    ${xtensor_INCLUDE_DIRS}
    ${xtensor-io_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    ${OMVS_INCLUDE_DIRS}
    ${RECONBENCH_INCLUDE_DIRS}
#    ${O3D_INCLUDE_DIRS}
    ${Open3D_INCLUDE_DIRS}
    ${PCL_INCLUDE_DIRS}
)

# link directories should not be necessary, see:
# https://cmake.org/cmake/help/v3.0/command/link_directories.html
#if(COLMAP_ENABLED)
#    link_directories(${COLMAP_LINK_DIRS})
#endif()
#if(O3D_ENABLED)
#    link_directories(${Open3D_LIBRARY_DIRS})
#endif()

message("CMAKE_CXX_FLAGS: " "${CMAKE_CXX_FLAGS}")
message("CMAKE_CXX_LINKER_FLAGS: " "${CMAKE_LINKER_FLAGS}")

add_library(rPLY
    external/rPLY/mesh.cpp
    external/rPLY/ply.cpp
    external/rPLY/rply.c)

add_library(util
    src/util/geometricOperations.cpp
    src/util/helper.cpp
    )
target_link_libraries(util
    ${Boost_LIBRARIES}
    CGAL::CGAL
    )

add_library(fileIO
    src/IO/fileIO.cpp
#    src/IO/ttIO.cpp
    src/IO/simplificationIO.cpp
    ${OpenMVS_INCLUDE_DIRS}
    ${OMVS_CPP_FILES})

target_link_libraries(fileIO
#    ${CNPY_LIBRARY}
#    ${CMAKE_DL_LIBS}
    ${OCV_EXTERNAL_LIBRARIES}
    ${OMVS_EXTERNAL_LIBRARIES}
    ${Open3D_LIBRARIES}
    )

add_library(learning
    src/learning/learning.cpp
#    src/learning/learningIO.cpp
    src/learning/learningIO_bin.cpp
    src/learning/learningMath.cpp
    src/learning/learningRayTracing.cpp
#    src/learning/learningRayTracingGroundTruth.cpp
    )
target_link_libraries(learning
    util
    ZLIB::ZLIB)

add_library(processing
    src/processing/evaluation.cpp
    src/processing/rayTracingFacet.cpp
    src/processing/meshProcessing.cpp
    src/processing/graphCut.cpp
    src/processing/pointSetProcessing.cpp
    src/processing/normalAndSensorProcessing.cpp
    src/processing/field.cpp
    )
target_link_libraries(processing
    util)

if(LABATUT)
    add_executable(labatut
        src/exe/labatut.cpp
        src/exe/inputParser.cpp
        )
    target_link_libraries(labatut
                        fileIO
                        rPLY
                        util
                        processing
                        ${OCV_EXTERNAL_LIBRARIES}
                        ${OMVS_EXTERNAL_LIBRARIES} # this has to be written before OpenMP for some reason
                        ${RECONBENCH_EXTERNAL_LIBRARIES}
                        ${PCL_EXTERNAL_LIBRARIES}
                        ${Boost_LIBRARIES}
                        ${CMAKE_DL_LIBS}
                        ${GLOG_LIBRARIES}
                        ${CERES_LIBRARIES}
                        ${GCO_LIBRARY}
                        ${GMP_LIBRARY}
                        ${MPFR_LIBRARY}
                        ${Open3D_LIBRARIES}
                        CGAL::CGAL
                        OpenMP::OpenMP_CXX
                        ZLIB::ZLIB
                        xtensor
                        xtensor-io
                        )
endif()

if(OCC2MESH)
    add_executable(occ2mesh
        src/exe/occ2mesh.cpp
        src/exe/inputParser.cpp
        )
    target_link_libraries(occ2mesh
                        fileIO
                        rPLY
                        learning
                        util
                        processing
                        ${OCV_EXTERNAL_LIBRARIES}
                        ${OMVS_EXTERNAL_LIBRARIES} # this has to be written before OpenMP for some reason
                        ${RECONBENCH_EXTERNAL_LIBRARIES}
                        ${PCL_EXTERNAL_LIBRARIES}
                        ${Boost_LIBRARIES}
                        ${CMAKE_DL_LIBS}
                        ${GLOG_LIBRARIES}
                        ${CERES_LIBRARIES}
                        ${GCO_LIBRARY}
                        ${GMP_LIBRARY}
                        ${MPFR_LIBRARY}
                        ${Open3D_LIBRARIES}
                        CGAL::CGAL
                        OpenMP::OpenMP_CXX
                        ZLIB::ZLIB
                        xtensor
                        xtensor-io
                        )
endif()

if(FEAT)
    add_executable(feat
        src/exe/feat.cpp
        src/exe/inputParser.cpp
        ${OMVS_CPP_FILES}
        ${OMVS_INCLUDE_DIRS}
        )
    target_link_libraries(feat
                        fileIO
                        learning
                        rPLY
                        util
                        processing
                        ${OCV_EXTERNAL_LIBRARIES}
                        ${OMVS_EXTERNAL_LIBRARIES} # this has to be written before OpenMP for some reason
                        ${RECONBENCH_EXTERNAL_LIBRARIES}
                        ${PCL_EXTERNAL_LIBRARIES}
                        ${Boost_LIBRARIES}
                        ${GLOG_LIBRARIES}
                        ${CERES_LIBRARIES}
                        ${GCO_LIBRARY}
                        ${GMP_LIBRARY}
                        ${MPFR_LIBRARY}
                        ${Open3D_LIBRARIES}
                        CGAL::CGAL
                        OpenMP::OpenMP_CXX
                        ZLIB::ZLIB
                        xtensor
                        xtensor-io
                        )
endif()

if(OMVS2NPZ)
    add_executable(omvs2npz
        src/exe/omvs2npz.cpp
        src/exe/inputParser.cpp
        ${OMVS_CPP_FILES}
        ${OMVS_INCLUDE_DIRS}
        )
    target_link_libraries(omvs2npz
                        fileIO
                        learning
                        rPLY
                        util
                        processing
                        ${OCV_EXTERNAL_LIBRARIES}
                        ${OMVS_EXTERNAL_LIBRARIES} # this has to be written before OpenMP for some reason
                        ${RECONBENCH_EXTERNAL_LIBRARIES}
                        ${PCL_EXTERNAL_LIBRARIES}
                        ${Boost_LIBRARIES}
                        ${GLOG_LIBRARIES}
                        ${CERES_LIBRARIES}
                        ${GCO_LIBRARY}
                        ${GMP_LIBRARY}
                        ${MPFR_LIBRARY}
                        ${Open3D_LIBRARIES}
                        CGAL::CGAL
                        OpenMP::OpenMP_CXX
                        xtensor
                        xtensor-io
    #                    OpenMVS::MVS
                        )
endif()

if(SCAN)
    add_executable(scan
        src/exe/scan.cpp
        src/exe/inputParser.cpp
        )
    target_link_libraries(scan
                        fileIO
                        rPLY
                        util
                        processing
                        learning
                        ${Boost_LIBRARIES}
                        ${GMP_LIBRARY}
                        ${MPFR_LIBRARY}
                        CGAL::CGAL
                        ${OMVS_EXTERNAL_LIBS}
                        OpenMP::OpenMP_CXX
                        ${CMAKE_DL_LIBS}
                        )
endif()



#add_executable(vsa
#    src/exe/vsa.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(vsa
#                    fileIO
#                    rPLY
#                    util
#                    processing
#                    ${Boost_LIBRARIES}
#                    ${GMP_LIBRARY}
#                    ${MPFR_LIBRARY}
#                    CGAL::CGAL
#                    ${OMVS_EXTERNAL_LIBS}
#                    OpenMP::OpenMP_CXX
#                    )



#add_executable(collapse
#    src/exe/collapse.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(collapse
#                    fileIO
#                    rPLY
#                    util
#                    processing
#                    ${Boost_LIBRARIES}
#                    ${GMP_LIBRARY}
#                    ${MPFR_LIBRARY}
#                    CGAL::CGAL
#                    ${OMVS_EXTERNAL_LIBS}
#                    OpenMP::OpenMP_CXX
#                    )



#add_executable(eval
#    src/exe/eval.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(eval
#                    fileIO
#                    rPLY
#                    util
#                    processing
#                    ${Boost_LIBRARIES}
#                    ${GMP_LIBRARY}
#                    ${MPFR_LIBRARY}
#                    CGAL::CGAL
#                    ${OMVS_EXTERNAL_LIBS}
#                    OpenMP::OpenMP_CXX
#                    )



#add_executable(iso
#    src/exe/iso.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(iso
#                    fileIO
#                    rPLY
#                    util
#                    processing
#                    learning
#                    ${Boost_LIBRARIES}
#                    ${GMP_LIBRARY}
#                    ${MPFR_LIBRARY}
#                    CGAL::CGAL
#                    ${OMVS_EXTERNAL_LIBS}
#                    OpenMP::OpenMP_CXX
#                    )

#add_executable(sample
#    src/exe/sample.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(sample
#                    fileIO
#                    rPLY
#                    util
#                    processing
#                    ${Boost_LIBRARIES}
#                    ${GMP_LIBRARY}
#                    ${MPFR_LIBRARY}
#                    CGAL::CGAL
#                    OpenMP::OpenMP_CXX
#                    )


#add_executable(boundingBox
#    src/exe/boundingBox.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(boundingBox
#                    fileIO
#                    rPLY
#                    util
#                    ${Boost_LIBRARIES}
#                    ${GMP_LIBRARY}
#                    ${MPFR_LIBRARY}
#                    CGAL::CGAL
#                    OpenMP::OpenMP_CXX
#                    )

#add_executable(simplification
#    src/exe/simplification.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(simplification
#                    fileIO
#                    rPLY
#                    util
#                    processing
#                    ${GCO_LIBRARY}
#                    ${Boost_LIBRARIES}
#                    ${GMP_LIBRARY}
#                    ${MPFR_LIBRARY}
#                    ${Open3D_LIBRARIES}
#                    CGAL::CGAL
#                    ${OMVS_EXTERNAL_LIBS}
#                    OpenMP::OpenMP_CXX
#                    )


#add_executable(normal
#    src/exe/normalEstimation.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(normal
#                        fileIO
#                        rPLY
#                        processing
#                        util
#                        ${PCL_EXTERNAL_LIBRARIES}
#                        ${Boost_LIBRARIES}
#                        ${GMP_LIBRARY}
#                        ${MPFR_LIBRARY}
#                        CGAL::CGAL
#                        OpenMP::OpenMP_CXX
#                    )


#add_executable(eth3d
#    src/exe/eth3d.cpp
#    src/exe/inputParser.cpp
#    src/IO/ethIO.cpp
#    src/IO/tinyxml2.cpp
#    src/IO/meshlab_project.cc
#    )
#target_link_libraries(eth3d
#                        fileIO
#                        rPLY
#                        util
#                        processing
#                        ${OCV_EXTERNAL_LIBRARIES}
#                        ${OMVS_EXTERNAL_LIBRARIES} # this has to be written before OpenMP for some reason
#                        ${PCL_EXTERNAL_LIBRARIES}
#                        ${Boost_LIBRARIES}
#                        ${GMP_LIBRARY}
#                        ${MPFR_LIBRARY}
#                        CGAL::CGAL
#                        OpenMP::OpenMP_CXX
#                    )


#add_executable(stats
#    src/exe/stats.cpp
#    src/exe/inputParser.cpp
#    )
#target_link_libraries(stats
#                        fileIO
#                        rPLY
#                        processing
#                        util
#                        ${Boost_LIBRARIES}
#                        ${GMP_LIBRARY}
#                        ${MPFR_LIBRARY}
#                        CGAL::CGAL
#                        OpenMP::OpenMP_CXX
#                    )

#add_executable(poisson
#    src/poisson/poisson.cpp
#    )
#target_link_libraries(poisson
#                    fileIO
#                    rPLY
#                    util
#                    processing
#                    ${Boost_LIBRARIES}
#                    ${GMP_LIBRARY}
#                    ${MPFR_LIBRARY}
#                    ${CGAL_LIBRARY}
#                    )




