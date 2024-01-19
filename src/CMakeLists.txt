cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD 20)

project(SPH_SOLVER)

# Set binary path
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../exec/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(SOURCES
    sph.cpp
    sph.h
    initial_conditions.cpp
    initial_conditions.h
    main_prog_funcs.h)

find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    add_executable(SPH_SOLVER SPH-main.cpp ${SOURCES})
    target_link_libraries(SPH_SOLVER PRIVATE Boost::program_options)
endif()
