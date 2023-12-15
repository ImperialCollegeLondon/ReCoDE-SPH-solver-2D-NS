cmake_minimum_required(VERSION 3.26)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(SPH_SOLVER)

# Set binary path
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(SOURCES
    src/class.cpp
    src/class.h)

find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    add_executable(SPH_SOLVER src/SPH-main.cpp ${SOURCES})
    target_link_libraries(SPH_SOLVER PRIVATE Boost::program_options)
endif()
