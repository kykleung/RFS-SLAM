
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#
# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
#
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

PROJECT(RFSSLAM)
SET(RFSSLAM_MAJOR_VERSION 1)
SET(RFSSLAM_MINOR_VERSION 1)
SET(RFSSLAM_PATCH_VERSION 0)
SET(RFSSLAM_VERSION ${RFSSLAM_MAJOR_VERSION}.${RFSSLAM_MINOR_VERSION}.${RFSSLAM_PATCH_VERSION})
message(STATUS "RFSSLAM Version ${RFSSLAM_VERSION}")

# Compiler settings
# Note: to change compiler, use for example...
# cmake -D CMAKE_C_COMPILER=gcc CMAKE_CXX_COMPILER=c++ .
SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native ")
#MESSAGE("CMAKE_CXX_COMPILTER_ID = ${CMAKE_CXX_COMPILER_ID}")
IF(APPLE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-qualifiers ") # For Eigen const warnings
  IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    MESSAGE("NOTE: Libraries being linked to RFSSLAM must also be compiled using GCC")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-enum-compare ") # For Eigen enum compare warnings
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q ") #required for march=native errors. This forces clang assembler to be used, and will produce some argument unused warnings. 
  ENDIF()
ENDIF()

# OpenMP
IF(APPLE AND ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"))
  SET(USE_OPENMP OFF) # OpenMP does not yet work with Clang, but hopefully soon
ELSE()
  SET(USE_OPENMP ON CACHE BOOL "Use OpenMP to compile multithreaded versions of SLAM algorithms.")
ENDIF()
IF(USE_OPENMP)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
ENDIF()

# Boost
FIND_PACKAGE(Boost COMPONENTS timer chrono system filesystem graph program_options REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME ON)

# For cmake custom library searches
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Eigen
FIND_PACKAGE(Eigen REQUIRED)

# Google perftools
FIND_PACKAGE(Perftools)
IF(Perftools_FOUND)
  IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(USE_CPU_PROFILER OFF CACHE BOOL "Use Google Perftools for CPU profiling.")
    SET(USE_HEAP_PROFILER OFF CACHE BOOL "Use Google Perftools for heap profiling.")
  ELSE()
    SET(USE_CPU_PROFILER OFF) # Does not play nice with Clang
    SET(USE_HEAP_PROFILER OFF)
  ENDIF()
  IF(USE_CPU_PROFILER)
    ADD_DEFINITIONS(-D_PERFTOOLS_CPU)
    IF(APPLE)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-no_pie ") # Required due to ASLR
    ENDIF()
  ENDIF()
  IF(USE_HEAP_PROFILER)
    ADD_DEFINITIONS(-D_PERFTOOLS_HEAP)
    IF(APPLE)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-no_pie ") # Required due to ASLR
    ENDIF()
  ENDIF()
ENDIF()

# GTest
# If using Ubuntu apt-get to install, go to /usr/src/gtest,
# run cmake . and make  . After compiling , 
# move libgtest_main.a and libgtest.a to /usr/lib
enable_testing()
find_package(GTest REQUIRED)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)
#MESSAGE("CMAKE_RUNTIME_OUTPUT_DIRECTORY = " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ./lib)
#MESSAGE("CMAKE_LIBRARY_OUTPUT_DIRECTORY = " ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./lib)
#MESSAGE("CMAKE_ARCHIVE_OUTPUT_DIRECTORY = " ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

FILE(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples")
FILE(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/bin/test")

#ADD_SUBDIRECTORY()
INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(test)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
IF(GTEST_FOUND)
  INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
ENDIF()
#LINK_DIRECTORIES()

ADD_LIBRARY(rfsslam 
  src/TimeStamp.cpp 
  src/Timer.cpp
  src/Frame.cpp
  src/ProcessModel_Odometry1D.cpp 
  src/ProcessModel_Odometry2D.cpp
  src/ProcessModel_Ackerman2D.cpp
  src/MeasurementModel_Rng1D.cpp
  src/MeasurementModel_RngBrg.cpp
  src/MeasurementModel_XY.cpp
  src/MeasurementModel_VictoriaPark.cpp
  src/KalmanFilter_RngBrg.cpp
  src/HungarianMethod.cpp
  src/MurtyAlgorithm.cpp
  src/BruteForceAssignment.cpp
  src/CostMatrix.cpp
  src/PermutationLexicographic.cpp
  src/MatrixPermanent.cpp
  src/misc/memProfile.cpp
  )
TARGET_LINK_LIBRARIES(rfsslam ${Boost_LIBRARIES})

ADD_EXECUTABLE(rbphdslam2dSim src/rbphdslam2dSim.cpp)
TARGET_LINK_LIBRARIES(rbphdslam2dSim ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(fastslam2dSim src/fastslam2dSim.cpp)
TARGET_LINK_LIBRARIES(fastslam2dSim ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY}) 

ADD_EXECUTABLE(rbphdslam_VictoriaPark src/rbphdslam_VictoriaPark.cpp)
TARGET_LINK_LIBRARIES(rbphdslam_VictoriaPark ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(fastslam_VictoriaPark src/fastslam_VictoriaPark.cpp)
TARGET_LINK_LIBRARIES(fastslam_VictoriaPark ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(analysis2dSim src/analysis2dSim.cpp)
TARGET_LINK_LIBRARIES(analysis2dSim ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(spatialIndexTree src/examples/spatialIndexTree.cpp)
TARGET_LINK_LIBRARIES(spatialIndexTree ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(convertLogFiles src/convertLogFiles.cpp)
TARGET_LINK_LIBRARIES(convertLogFiles ${Boost_LIBRARIES})

IF(GTEST_FOUND)
  
  ADD_EXECUTABLE(spatialIndexTreeTest test/SpatialIndexTreeTest.cpp)
  TARGET_LINK_LIBRARIES(spatialIndexTreeTest ${GTEST_LIBRARIES} ${Boost_LIBRARIES} rfsslam pthread)
  ADD_TEST(unitTest spatialIndexTreeTest)

  ADD_EXECUTABLE(matrixPermanentTest test/MatrixPermanentTest.cpp)
  TARGET_LINK_LIBRARIES(matrixPermanentTest ${GTEST_LIBRARIES} rfsslam pthread)
  ADD_TEST(unitTest matrixPermanentTest)

ENDIF()

# Examples

ADD_EXECUTABLE(linearAssignment_CostMatrixPartitioning src/examples/linearAssignment_CostMatrixPartitioning.cpp)
TARGET_LINK_LIBRARIES(linearAssignment_CostMatrixPartitioning
${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(linearAssignment_LexicographicOrdering src/examples/linearAssignment_LexicographicOrdering.cpp)
TARGET_LINK_LIBRARIES(linearAssignment_LexicographicOrdering ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(linearAssignment_MurtyAlgorithm src/examples/linearAssignment_MurtyAlgorithm.cpp)
TARGET_LINK_LIBRARIES(linearAssignment_MurtyAlgorithm ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(ospaError src/examples/ospaError.cpp)
TARGET_LINK_LIBRARIES(ospaError ${Boost_LIBRARIES} rfsslam)

set_target_properties(linearAssignment_CostMatrixPartitioning
	linearAssignment_LexicographicOrdering
	linearAssignment_MurtyAlgorithm 
	ospaError
	spatialIndexTree
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples"
)


# Installation

INSTALL(TARGETS rfsslam
  DESTINATION /usr/local/lib/rfsslam
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
  EXPORT rfsslam-targets)
INSTALL(EXPORT rfsslam-targets DESTINATION /usr/local/lib/rfsslam)
CONFIGURE_FILE(rfsslam-config.cmake.in "rfsslam-config.cmake" @ONLY)
CONFIGURE_FILE(rfsslam-config-version.cmake.in "rfsslam-config-version.cmake" @ONLY)
INSTALL(FILES ${CMAKE_BINARY_DIR}/rfsslam-config.cmake ${CMAKE_BINARY_DIR}/rfsslam-config-version.cmake
  DESTINATION /usr/local/lib/rfsslam)

FILE(GLOB headerFiles "include/*.hpp")
INSTALL(FILES ${headerFiles}
  DESTINATION /usr/local/include/rfsslam
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
