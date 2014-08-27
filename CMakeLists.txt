# Keith Leung 2013

# Note: to change compiler, use...
# cmake -D CMAKE_C_COMPILER=gcc CMAKE_CXX_COMPILER=c++ .

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE Release)
PROJECT(RFSSLAM)
SET(RFSSLAM_MAJOR_VERSION 1)
SET(RFSSLAM_MINOR_VERSION 1)
SET(RFSSLAM_PATCH_VERSION 0)
SET(RFSSLAM_VERSION ${RFSSLAM_MAJOR_VERSION}.${RFSSLAM_MINOR_VERSION}.${RFSSLAM_PATCH_VERSION})
message(STATUS "RFSSLAM Version ${RFSSLAM_VERSION}")

# Explicitly set for colorgcc. By default, they point to gcc and c++ with absolute path
SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "c++")

# Disable of edit the following line according to your processor's support for SSE
# For X86-64 architecture, SSE2 is generally enabled by default
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-march=native")

FILE(MAKE_DIRECTORY "bin")
FILE(MAKE_DIRECTORY "bin/examples")
FILE(MAKE_DIRECTORY "bin/test")
FILE(MAKE_DIRECTORY "lib")
FILE(MAKE_DIRECTORY "data")

FIND_PACKAGE(Boost COMPONENTS math_c99 timer system thread filesystem graph REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
FIND_PACKAGE(Eigen 3.0.0 REQUIRED)

# Change the following include directories depending on where gtest is installed
# If using Ubuntu apt-get to install, go to /usr/src/gtest,
# run cmake . and make  . After compiling , 
# move libgtest_main.a and libgtest.a to /usr/lib
enable_testing()
find_package(GTest REQUIRED)


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)
MESSAGE("CMAKE_RUNTIME_OUTPUT_DIRECTORY = " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./lib)
MESSAGE("CMAKE_ARCHIVE_OUTPUT_DIRECTORY = " ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

#ADD_SUBDIRECTORY()
INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(test)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
#LINK_DIRECTORIES()

ADD_LIBRARY(rfsslam 
  src/TimeStamp.cpp 
  src/Timer.cpp
  src/Pose.cpp 
  src/Measurement_Odometry2D.cpp 
  src/ProcessModel_Odometry1D.cpp 
  src/ProcessModel_Odometry2D.cpp
  src/MeasurementModel_Rng1D.cpp
  src/MeasurementModel_RngBrg.cpp
  src/KalmanFilter_RngBrg.cpp
  src/HungarianMethod.cpp
  src/MurtyAlgorithm.cpp
  src/BruteForceAssignment.cpp
  src/CostMatrix.cpp
  src/PermutationLexicographic.cpp
  )

ADD_EXECUTABLE(rbphdslam2dSim src/rbphdslam2dSim.cpp)
TARGET_LINK_LIBRARIES(rbphdslam2dSim config++ ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(fastslam2dSim src/fastslam2dSim.cpp)
TARGET_LINK_LIBRARIES(fastslam2dSim config++ ${Boost_LIBRARIES} rfsslam)

#ADD_EXECUTABLE(unitTest test/test.cpp src/Measurement.cpp src/Pose.cpp)
#TARGET_LINK_LIBRARIES(unitTest ${GTEST_LIBRARIES} ${Boost_LIBRARIES} pthread TimeStamp MeasurementModel ProcessModel LinearAssignment)
#ADD_TEST(unitTest bin/unitTest)

ADD_EXECUTABLE(analysis2dSim src/analysis2dSim.cpp src/Pose.cpp)
TARGET_LINK_LIBRARIES(analysis2dSim ${Boost_LIBRARIES} rfsslam)

#ADD_EXECUTABLE(hmTest src/test.cpp)
#TARGET_LINK_LIBRARIES(hmTest LinearAssignment)

ADD_EXECUTABLE(convertLogFiles src/convertLogFiles.cpp)
TARGET_LINK_LIBRARIES(convertLogFiles ${Boost_LIBRARIES})

# Examples

ADD_EXECUTABLE(examples/linearAssignment_CostMatrixPartitioning src/examples/linearAssignment_CostMatrixPartitioning.cpp)
TARGET_LINK_LIBRARIES(examples/linearAssignment_CostMatrixPartitioning ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(examples/linearAssignment_LexicographicOrdering src/examples/linearAssignment_LexicographicOrdering.cpp)
TARGET_LINK_LIBRARIES(examples/linearAssignment_LexicographicOrdering ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(examples/linearAssignment_MurtyAlgorithm src/examples/linearAssignment_MurtyAlgorithm.cpp)
TARGET_LINK_LIBRARIES(examples/linearAssignment_MurtyAlgorithm ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(examples/ospaError src/examples/ospaError.cpp)
TARGET_LINK_LIBRARIES(examples/ospaError ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(test/cpu_timer_test src/tests/cpu_timer_info.cpp)
TARGET_LINK_LIBRARIES(test/cpu_timer_test ${Boost_LIBRARIES})

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
