cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)

project(SPH_SOLVER)

# Set binary path
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../exec/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# If the user doesn't specify, will build release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# General Flags
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -Og")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

# Setting profiling linking
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -pg)
    set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} -pg)
endif()

set(SOURCES
    particles.cpp
    fluid.cpp
    sph_solver.cpp
    initial_conditions.cpp)

set(HEADERS
    particles.h
    fluid.h
    sph_solver.h
    initial_conditions.h
    main_prog_funcs.h)

find_package(Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    add_executable(SPH_SOLVER SPH-main.cpp ${SOURCES} ${HEADERS})
    target_link_libraries(SPH_SOLVER PRIVATE Boost::program_options)
endif()
