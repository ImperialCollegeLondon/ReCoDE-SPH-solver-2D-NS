# Building with CMake

CMake is a tool that facilitates the build of large projects written in C++. It takes care of all the linking and dependency fetching for the user. Moreover, it provides a user-friendly scripting language format, enabling the user to specify the building steps using code. More information about how to use CMake can be found in the [official CMake page](https://cmake.org/cmake/help/latest/guide/tutorial/index.html).

## CMakeLists file

Everything related to CMake should be inside a file called `CMakeLists.txt`. Without it, `cmake` cannot start the procedure of building the project. There is a plethora of features in modern CMake and extensive documentation about it. For the purpose of this project, the CMake functionality that is required and used is limited, thus only some of these features are mentioned and described.

## CMake Flags

In the beginning of the `CMakeLists.txt`, it's possible to incorporate various flags. These flags can be related to CMake (e.g `cmake_minimum_required(VERSION 3.26)`) or to C++ (e.g `set(CMAKE_CXX_STANDARD 20)`). For example, the latter flag sets the C++ Standard for the CMake build.

## Project

You start the build (or project as called inside CMakeLists.txt by naming the final outcome of the build using project(<name>))

When the CMake building process completes, a binary output is produced. It is good practice to specify where this binary will reside when the build is done. This is done with the following lines:

```
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
```

Here, we store the generated binary within a folder of the CMake source directory, called _build_.

CMake has some internal variables that can be accessed using `${}`, as seen in the code above with `${CMAKE_SOURCE_DIR}`. More information about these variables can be found [here](https://cmake.org/cmake/help/latest/manual/cmake-variables.7.html).

## Sources

The `SOURCES` attribute specifies all the source files that are used in the build of the project. In the case that the C++ project contains many source and many header files, it is considered a good practice to specify them separately, as follows:

```
set(SOURCES
    <all_of_sources>)
set(INCLUDES
    <all_of_headers>)

add_executable(<project_name> <main_cpp> ${SOURCES} ${INCLUDES})
```

For this project, this is not necessary, as there are only five files.

## External Dependencies (Boost)

CMake is great in fetching and linking any external dependency that a C++ project might require. It's especially helpful with _[Boost](https://www.boost.org/)_, which is utilised by this project, as it can use `find_package` and fetch the existing installation. When finding packages, it is best practice to include only the headers that the user requires (in our case `program_options`).

In the case the _Boost_ package is found, we can add the output executable (`SPH_SOLVER`) and link it to the `SOURCE` files we have declared earlier (CMake will go and fetch them automatically).

For Boost to work, one more step is required. We need to link Boost's libraries to the target (`SPH_SOLVER`). This is done using the `target_link_libraries` command.
