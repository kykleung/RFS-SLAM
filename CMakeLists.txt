# Keith Leung 2013

# Note: to change compiler, use...
# cmake -D CMAKE_C_COMPILER=gcc CMAKE_CXX_COMPILER=c++ .

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE Release)
PROJECT(RFS_SLAM)

# Explicitly set for colorgcc. By default, they point to gcc and c++ with absolute path
SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "c++")

# Disable of edit the following line according to your processor's support for SSE
# For X86-64 architecture, SSE2 is generally enabled by default
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-msse4.2")

FILE(MAKE_DIRECTORY "bin")
FILE(MAKE_DIRECTORY "lib")
FILE(MAKE_DIRECTORY "data")

FIND_PACKAGE(Boost COMPONENTS math_c99 timer system thread REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME ON)

FIND_PACKAGE(Eigen REQUIRED)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)
MESSAGE("CMAKE_RUNTIME_OUTPUT_DIRECTORY = " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./lib)
MESSAGE("CMAKE_ARCHIVE_OUTPUT_DIRECTORY = " ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

#ADD_SUBDIRECTORY()
INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(test)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Eigen_INCLUDE_DIRS})
#LINK_DIRECTORIES()

ADD_LIBRARY(ProcessModel src/ProcessModel.cpp)
ADD_LIBRARY(MeasurementModel src/MeasurementModel.cpp)

ADD_EXECUTABLE(sim2d src/simulation2d.cpp src/Measurement.cpp src/Pose.cpp)
TARGET_LINK_LIBRARIES(sim2d config++ ${Boost_LIBRARIES} ProcessModel MeasurementModel)

