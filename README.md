<!-- Your Project title, make it sound catchy! -->

# RECODE-SPH-SOLVER-2D-NS

[![C++](https://img.shields.io/badge/C++-00599C?logo=c%2B%2B&logoColor=white)](https://cplusplus.com/)
[![Boost](https://img.shields.io/badge/Boost-00599C?logo=boost&logoColor=white)](https://www.boost.org/)
[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)](https://www.python.org/)

<!-- Provide a short description to your project -->

## Description

In this project we present a numerical code in C++ which solves the two-dimensional Navier-Stokes equations using the smoothed-particle hydrodynamics (SPH) approach. The focus lies in the implementation (and documentation) of good C++ practices and the development of skills related to efficient, robust, extensible and readable scientific code. The learning process regarding this project can be twofold:

1) The student can study the material provided in the `main` branch of the present repository. It is independent of all the other branches and can be used as a standalone educational resource. In this the implemented SPH methodology is explained as well as the structure of the source code and the post-processing scripts.

2) The student can start by studying progressively the branches `v0` - `v5` in order to experience the process which was followed in order to improve and optimize the herein code. Several comments have been added for each individual version in the corresponding branch to highlight the improvements which were implemented compared to its ancestors.

<!-- What should the students going through your exemplar learn -->

## Learning Outcomes

- I/O (input/output)
- OOP (object-oriented programming)
- C++ containers
- Performance and memory optimization tools and skills

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
│   ├── 2.Code_Overview.md
│   ├── 3.IO_Overview.md
│   ├── 4.OOP_Concepts.md
│   ├── 5.STL.md
│   ├── 6.Efficient_Programming.md
│   ├── 7.Neighbour_Search.md
│   ├── 8.Adaptive_Timestep.md
│   ├── 9.Code_Execution.md
│   ├── 10.GH_Codespaces.md
│   ├── 11.Profiling.md
│   ├── 12.Post_Analysis.md
│   ├── A1.ClangFormat.md
│   ├── A2.CMake.md
│   └── B1.Exercises.md
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
│   ├── simulation_animation.ipynb
│   ├── simulation_animation.py
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
<!-- readme: collaborators -start -->
<table>
<tr>
    <td align="center">
        <a href="https://github.com/cmcooling">
            <img src="https://avatars.githubusercontent.com/u/12444767?v=4" width="100;" alt="cmcooling"/>
            <br />
            <sub><b>Chris Cooling</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/ChrisPetalotis">
            <img src="https://avatars.githubusercontent.com/u/17175412?v=4" width="100;" alt="ChrisPetalotis"/>
            <br />
            <sub><b>Christos Petalotis</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/kmichali">
            <img src="https://avatars.githubusercontent.com/u/25386647?v=4" width="100;" alt="kmichali"/>
            <br />
            <sub><b>Katerina Michalickova</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/dc2917">
            <img src="https://avatars.githubusercontent.com/u/45606273?v=4" width="100;" alt="dc2917"/>
            <br />
            <sub><b>Dan Cummins</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/ByronAV">
            <img src="https://avatars.githubusercontent.com/u/54893888?v=4" width="100;" alt="ByronAV"/>
            <br />
            <sub><b>Vyron Avramidis</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Vasilis421">
            <img src="https://avatars.githubusercontent.com/u/62700983?v=4" width="100;" alt="Vasilis421"/>
            <br />
            <sub><b>Vasileios Christou</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/geoefsta">
            <img src="https://avatars.githubusercontent.com/u/62747725?v=4" width="100;" alt="geoefsta"/>
            <br />
            <sub><b>George Efstathiou</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/ge420">
            <img src="https://avatars.githubusercontent.com/u/116084464?v=4" width="100;" alt="ge420"/>
            <br />
            <sub><b>George Efstathiou</b></sub>
        </a>
    </td></tr>
</table>
<!-- readme: collaborators -end -->

<!-- readme: contributors -start -->
<table>
<tr>
    <td align="center">
        <a href="https://github.com/ge420">
            <img src="https://avatars.githubusercontent.com/u/116084464?v=4" width="100;" alt="ge420"/>
            <br />
            <sub><b>Null</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Vasilis421">
            <img src="https://avatars.githubusercontent.com/u/62700983?v=4" width="100;" alt="Vasilis421"/>
            <br />
            <sub><b>Vasileios Christou</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/ByronAV">
            <img src="https://avatars.githubusercontent.com/u/54893888?v=4" width="100;" alt="ByronAV"/>
            <br />
            <sub><b>Vyron Avramidis</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/ChrisPetalotis">
            <img src="https://avatars.githubusercontent.com/u/17175412?v=4" width="100;" alt="ChrisPetalotis"/>
            <br />
            <sub><b>Christos Petalotis</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/geoefsta">
            <img src="https://avatars.githubusercontent.com/u/62747725?v=4" width="100;" alt="geoefsta"/>
            <br />
            <sub><b>Null</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/cmcooling">
            <img src="https://avatars.githubusercontent.com/u/12444767?v=4" width="100;" alt="cmcooling"/>
            <br />
            <sub><b>Chris Cooling</b></sub>
        </a>
    </td></tr>
</table>
<!-- readme: contributors -end -->

<!-- Change this to your License. Make sure you have added the file on GitHub -->

## License

This project is licensed under the [BSD-3-Clause license](LICENSE.md)
