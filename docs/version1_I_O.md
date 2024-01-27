# Description of the code

## Overview

The material in `v1` builds upon the foundation of `v0` by introducing a refined approach to handle input and output (I/O) procedures within the code. We remove hardcoded values and instead incorporate data retrieved from external input files. Furthermore, the code responsible for managing output values and files has been restructured for improved organisation and maintainability. These enhancements enable greater flexibility in parameter manipulation, paving the way for a more generalisable code design, suitable for various domain boundaries and fluid parameters.

### Inputs

---

To facilitate this enhanced I/O functionality, the `initialise()` function has been extended to incorporate the process of reading additional parameters from newly introduced input files. The `<boost/program_options.hpp>` library is employed once again, enabling seamless mapping of read parameters to their corresponding values and assigning them to variables with identical names. Additionally, error handling is integrated into the input file reading process to guarantee that the provided values conform to the constraints imposed by the underlying mathematical models.

A set of new input files are introduced in `v1`, located within the `/exec/input/` directory, that empower users to tailor the execution of the program to their specific requirements. By accessing and modifying these files, users can manipulate various parameters, such as domain boundaries, fluid properties, and simulation settings. User-provided values are subjected to validation against the constraints imposed by the underlying mathematical models. If any discrepancies are detected, the implemented error handling mechanism gracefully handles the situation, preventing the program from proceeding with invalid configurations.

The newly introduced input files are the following:

- `constants.txt`: values that don't change throughout the program execution and define the behaviour of the simulation. These include:
  - radius of influence (`h`)
  - gas constant (`gas_constant`)
  - resting density (`density_resting`)
  - viscosity (`viscosity`)
  - acceleration due to gravity (`acceleration_gravity`)
  - coefficient of restituion (`coeff_restitution`)
- `domain.txt`: defines the dimensions of the domain utilised in the simulation
- `ic-one-particle.txt`: sets the initial positions when the selected initial condition is "ic-one-particle"
- `ic-{two, three, four}-particles.txt`: sets the initial positions for the corresponding cases
- `ic-block-drop.txt`: sets initial conditions for SPH, including the number of particles, the length and width of the block, and the initial axes positions for the center of the block
- `ic-droplet.txt`: sets initial conditions for SPH, including the number of particles, the size of the radius of the droplet, and the initial axes positions for the center of the droplet

### Outputs

---

`v1` introduces significant improvements to the management of program outputs, enhancing both data organisation and post-processing efficiency.

One key change is the adoption of CSV format for output files, replacing the previous .txt format. This switch not only optimises data storage, but also improves readability and facilitates data manipulation. Additionally, output files are now stored in a centralised location, specifically within the `/exec/output/` directory. This centralisation simplifies data organisation and retrieval, making it easier for users to access and analyse output data.

Another notable improvement is the restructuring of the code related to output handling. The code has been divided into more granular functions, promoting code modularity and maintainability. This refactoring also streamlines the flow of logic in the output handling process, making it more organised and easier to understand.

With the transition to .csv output format, the scripts for post-processing and visualisation were refactored to incorporate less, more efficient code, reducing code redundancy and improving the efficiency of the code, consequently making post-processing operations more streamlined and resource-efficient.

Furthermore, a new output file is, now, created during program execution, specifically to store the initial positions of the particles as defined in the input files. The name of this file is `initial-position.txt`. This addition provides a valuable reference point for analysing the simulation's evolution, allowing users to accurately visualise the initial state of the system.
