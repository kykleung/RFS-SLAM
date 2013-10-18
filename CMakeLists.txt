# Keith Leung 2013

# Note: to change compiler, use...
# cmake -D CMAKE_C_COMPILER=gcc CMAKE_CXX_COMPILER=c++ .

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE Release)
PROJECT(RB_PHD_FILTER)

# Explicitly set for colorgcc. By default, they point to gcc and c++ with absolute path
SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "c++")

# Disable of edit the following line according to your processor's support for SSE
# For X86-64 architecture, SSE2 is generally enabled by default
#SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-msse4.2 -O0")
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-msse4.2")

FILE(MAKE_DIRECTORY "bin")
FILE(MAKE_DIRECTORY "lib")

FIND_PACKAGE(Boost COMPONENTS math_c99 timer system thread REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME ON)

FIND_PACKAGE(Eigen REQUIRED)

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
INCLUDE_DIRECTORIES(${Eigen_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
#LINK_DIRECTORIES()

ADD_EXECUTABLE(laTest src/test.cpp)

ADD_LIBRARY(ProcessModel src/ProcessModel.cpp)
ADD_LIBRARY(MeasurementModel src/MeasurementModel.cpp)

ADD_EXECUTABLE(unitTest test/test.cpp src/Measurement.cpp src/Pose.cpp)
TARGET_LINK_LIBRARIES(unitTest ${GTEST_LIBRARIES} ${Boost_LIBRARIES} pthread MeasurementModel ProcessModel)
ADD_TEST(unitTest bin/unitTest)

ADD_EXECUTABLE(sim2d src/simulation2d.cpp src/Measurement.cpp src/Pose.cpp)
TARGET_LINK_LIBRARIES(sim2d config++ ${Boost_LIBRARIES} ProcessModel MeasurementModel)

ADD_EXECUTABLE(sim1d src/simulation1d.cpp src/Measurement.cpp src/Pose.cpp)
TARGET_LINK_LIBRARIES(sim1d config++ ${Boost_LIBRARIES} ProcessModel MeasurementModel)

