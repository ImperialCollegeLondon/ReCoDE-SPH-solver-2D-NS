<!-- Your Project title, make it sound catchy! -->
<div align="center">

# RECODE-SPH-SOLVER-2D-NS

[![C++](https://img.shields.io/badge/C++-00599C?logo=c%2B%2B&logoColor=white)](https://cplusplus.com/)
[![MPI](https://img.shields.io/badge/MPI-003366?logo=mpi&logoColor=white)](https://www.mpi-forum.org/)
[![Boost](https://img.shields.io/badge/Boost-00599C?logo=boost&logoColor=white)](https://www.boost.org/)
[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)](https://www.python.org/)


<!-- Provide a short description to your project -->
</div>

## Description

In this project we present a numerical code in C++ which solves the two-dimensional Navier-Stokes equations using the smoothed-particle hydrodynamics (SPH) approach. The focus lies in the implementation (and documentation) of good C++ practices and the development of skills related to efficient, robust, extensible and readable scientific code. The learning process regarding this project can be twofold:

1) The student can study the material provided in the `main` branch of the present repository. It is independent of all the other branches and can be used as a standalone educational resource. In this the implemented SPH methodology is explained as well as the structure of the source code and the post-processing scripts.

2) The student can start by studying progressively the branches `v0` - `v5` in order to experience the process which was followed in order to improve and optimize the herein code. Several comments have been added for each individual version in the corresponding branch to highlight the improvements which were implemented compared to its ancestors.

<!-- What should the students going through your exemplar learn -->

## Learning Outcomes

- Advanced I/O (input/output)
- OOP (object-oriented programming)
- C++ containers
- Performance and memory optimization tools and skills
- Parallel programming using MPI

<!-- How long should they spend reading and practising using your Code.
Provide your best estimate -->

| Task       | Time    |
| ---------- | ------- |
| Reading    | 10 hours |
| Practicing | 4 hours |

## Requirements





<!--
If your exemplar requires students to have a background knowledge of something
especially this is the place to mention that.

List any resources you would recommend to get the students started.

If there is an existing exemplar in the ReCoDE repositories link to that.
-->

### Academic

- Experience with basic programming concepts (for loops, functions, reading and writing files etc.).

- Some experience with C++ (familiarity with pointers, C++ classes and the use of external libraries).

- Basic understanding of numerical analysis concepts (time marching, temporal integration etc.)




<!-- List the system requirements and how to obtain them, that can be as simple
as adding a hyperlink to as detailed as writing step-by-step instructions.
How detailed the instructions should be will vary on a case-by-case basis.

Here are some examples:

- 50 GB of disk space to hold Dataset X
- Anaconda
- Python 3.11 or newer
- Access to the HPC
- PETSc v3.16
- gfortran compiler
- Paraview
-->

### System

<!-- Instructions on how the student should start going through the exemplar.

Structure this section as you see fit but try to be clear, concise and accurate
when writing your instructions.

For example:
Start by watching the introduction video,
then study Jupyter notebooks 1-3 in the `intro` folder
and attempt to complete exercise 1a and 1b.

Once done, start going through through the PDF in the `main` folder.
By the end of it you should be able to solve exercises 2 to 4.

A final exercise can be found in the `final` folder.

Solutions to the above can be found in `solutions`.
-->

## Getting Started

<!-- An overview of the files and folder in the exemplar.
Not all files and directories need to be listed, just the important
sections of your project, like the learning material, the code, the tests, etc.

A good starting point is using the command `tree` in a terminal(Unix),
copying its output and then removing the unimportant parts.

You can use ellipsis (...) to suggest that there are more files or folders
in a tree node.

-->

## Project Structure

```log
.
├── docs
│   ├── index.md
│   ├── 1.SPH.md
│   ├── 2.Code_overview.md
│   ├── 3.Code_execution.md
│   ├── 4.IO_overview.md
│   ├── 5.OOP_concepts.md
│   ├── A1.ClangFormat.md
│   └── A2.CMake.md
├── exec
│   ├── build
│   └── input
│       ├── case.txt
│       ├── constants.txt
│       ├── domain.txt
│       ├── ic-block-drop.txt
│       ├── ic-droplet.txt
│       ├── ic-one-particles.txt
│       ├── ic-two-particles.txt
│       ├── ic-three-particles.txt
│       └── ic-four-particles.txt
├── notebooks
├── post
│   ├── plot_energies.ipynb
│   ├── plot_energies.py
│   ├── visualise_particles.ipynb
│   └── visualise_particles.py
└── src
    ├── CMakeLists.txt
    ├── fluid.cpp
    ├── fluid.h
    ├── initial_conditions.cpp
    ├── initial_conditions.h
    ├── main_prog_funcs.h
    ├── particles.cpp
    ├── particles.h
    ├── SPH-main.cpp
    ├── sph_solver.cpp
    └── sph_solver.h
```

<!-- Change this to your License. Make sure you have added the file on GitHub -->

## License

This project is licensed under the [BSD-3-Clause license](LICENSE.md)
