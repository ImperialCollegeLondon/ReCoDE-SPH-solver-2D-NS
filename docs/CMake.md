# Building with CMake

CMake is a really useful tool that helps in building large projects written in C++. It takes care of all the linking and dependency fetching for the user. It also provides a more user-friendly scripting language format in order to basically "code" the building steps. More information can be found in the official page: https://cmake.org/cmake/help/latest/guide/tutorial/index.html

## CMakeLists file

Everything related to CMake should be inside the `CMakeLists.txt` file. Without it, `cmake` cannot start the procedure of building the project. There are many features in modern CMake and massive documentation about it but we're only going to mention some of them here because the current project is not long enough to need additional functionality.

## CMake Flags

In the beginning of the `CMakeLists.txt` we can add the flags. These flags could be related to CMake (e.g `cmake_minimum_required(VERSION 3.26)`) or they can be related to C++ (e.g `set(CMAKE_CXX_STANDARD 20)`). For example, the latter flag will set the C++ Standard for the CMake build. 

## Project

You start the build (or `project` as called inside `CMakeLists.txt` by naming the final outcome of the build using `project(<name>)`). Because the building procedure will create an output, it is good practise to specify where this binary will reside
```
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
```
At this point we should also mention that CMake has some internal variables that can be accessed using `${}` as it can be seen above with `${CMAKE_SOURCE_DIR}`. There are many variables like this, more information can be found here: https://cmake.org/cmake/help/latest/manual/cmake-variables.7.html

## Sources

The `SOURCES` determine all the source files that will be used in the building of the project. In the case that the project contained many source and many header files, it would be a good practise to break those two appart into something like this:
```
set(SOURCES
    <all_of_sources>)
set(INCLUDES
    <all_of_headers>)

add_executable(<project_name> <main_cpp> ${SOURCES} ${INCLUDES})
```

For our purpose that is unesseccary since we only have two files.

## External Dependencies (Boost)

CMake is great in fetching and linking any external dependency that the project might contain. It's especially helpful with Boost since it can use `find_package` and fetch the existing installation. When finding packages, it is best practise to include only the headers that the user would require (in our case `program_options`). In the case that Boost package is found, we can add the output executable (`SPH_SOLVER`) and link it to the `SOURCE` files we have declared earlier (CMake will go and fetch them automatically). For Boost to work though we have another step, we need to link Boost's libraries to the target (`SPH_SOLVER`). This can be done with the `target_link_libraries` command.

## Usage

The most used way to build with CMake and then compile a project is the following.
1. Create a `build` directory inside the project root directory and then `cd` in it:
```
mkdir build; cd build
```
2. Build with cmake the target on the root directory (utilising `CMakeLists.txt`):
```
cmake ../src
```
3. Compile inside `build` directory using CMake
```
cmake --build .
```

That's it! The project binary can then be accessed (as long as the building and compilation finished smoothly) and executed inside `build` directory.