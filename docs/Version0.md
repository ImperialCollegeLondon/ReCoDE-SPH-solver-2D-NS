# Description of the code

## Overview

The code in branch `v0` contains a serial C++ implementation of the algorithm described in `SPH.md`. The variables associated with the particles' positions, velocities and forces, are stored as data members of an object called `SPH`. The member functions of the `SPH` class manifest the functions that embody the steps of the algorithm. While the current version of the code produces accurate results, subsequent chapters will enhance its capabilities and improve its structure. The functionality and rationale behind each improvement will be analysed in detail.

## Compiling and executing the code

The list of requirements for the v0 code is:

- A `C++20` version
- The `Boost` library

To compile the code, the user has to change the working directory to `exec/build` and then run the following in the terminal:

```
cmake ../../src
```

and then:

```
cmake --build .
```

This produces an executable file, called `SPH-SOLVER`, in the `exec/build` folder. To execute the code, the user needs to run the following:

```
./SPH-SOLVER
```

To clean the `build` directory, the user can use the following command:

```
cmake --build . --target clean
```

This will effectively delete the binary from the `build` directory.

## Files

This version of the code displays a serial implementation of the SPH algorithm in C++. It comprises three `*.cpp` files and their corresponding header (`*.h`) files. The code is accompanied by two input text files, one for the input variables of the executed case and the input parameters and one for the domain boundaries.

- `src/SPH-main.cpp`
- `src/sph.{h, cpp}`
- `src/initial_conditions.{h, cpp}`
- `src/main_prog_func.h`
- `src/CMakeLists.txt`
- `exec/inputs/case.txt`

## Main program

The main program of the code (i.e.`SPH-main.cpp`) is used to read the input files, export the output files, initialise the `SPH` class and perform the time integration by calling `SPH`'s member functions.

## Structure of the class

The `SPH` class is initialised by using the number of particles (`N`) which is required to determine the size of the arrays in the constructor. The operator `()` has been overloaded to place the particles in their initial conditions and to set their initial velocities to the corresponding arrays. The class has three main functions for the temporal integration:

- `SPH::calculateParticleDistance()`: Calculates the distances between the particles.

- `SPH::calculateDensity()`: Updates the density.

- `SPH::particleIterations()`: Calculates all the forces acting on the particles and updates their positions based on the leapfrog scheme.

## Input Parameters

### Initial conditions

As stated earlier, the `SPH` class is initialised in the `SPH-main.cpp` file, where one of the following initial conditions can be specified by the user from the `exec/inputs/case.txt` file:

- A single particle (`ic-one-particle`) at : (0.5, 0.5) to test the correctness of time integration and gravity forcing, as well as the bottom boundary condition.

- Two particles (`ic-two-particles`) at : (0.5, 0.5) and (0.5, h) to assess the pressure force and viscous terms.

- Three particles (`ic-three-particles`) at : (0.5, 0.5), (0.495, h) and (0.505, h) to assess left and right boundary conditions.

- Four particles (`ic-four-particles`) at : (0.505, 0.5), (0.515, 0.5), (0.51, 0.45) and (0.5, 0.45) to assess multiple particle interaction.

- A Dam break (`ic-dam-break`): a grid of particles occupying the region $[0,0.2]^2$.

- A Block drop (`ic-block-drop`): a grid of particles occupying the region $[0.1,0.3]\times[0.3,0.6]$.

- A Droplet (`ic-droplet`): particles occupying a circle of radius 0.1, centered at the point $[0.5,0.7]$

### Simulation Duration

The real simulated time can be set in units of seconds with the input key `T`. For instance, a simulation of 10 seconds can be set by typing `T = 10` in the `exec/inputs/case.txt` file.

### Time-step

The time-step `dt` can be set in units of seconds with the input key `dt`. For instance, a time-step of 0.0001 seconds can be set by typing `dt = 1e-4` in the `exec/inputs/case.txt` file.

### Radius of Influence

The last input parameter that the user can specify in the `exec/inputs/case.txt` file is the radius of influence `h` in units of metres. This parameter dictates the maximum distance in which an entity (particle or domain boundary) has an influence in the behaviour of another entity.

### Reading Inputs

The program expects the aforementioned parameters. When the function `initialize()` is invoked by the main program to read the `exec/inputs/case.txt` file, the `<boost/program_options.hpp>` library is utilised. This library facilitates the mapping of these parameters to their corresponding values, which are then stored in their respective variables. This approach enhances the flexibility and robustness of the input reading process. Users can specify input parameters in the `exec/inputs/case.txt` file in any order, provided they are presented as `key = value` pairs.

```cpp
// Process to obtain the directions provided by the user
po::options_description desc("Allowed options");
desc.add_options()("init_condition", po::value<std::string>(),
                    "take an initial condition")("T", po::value<double>(),
                                                "take integration time")(
    "dt", po::value<double>(), "take time-step")("h", po::value<double>(),
                                                "take radius of influence");

po::variables_map vm;
std::ifstream inputFile;
inputFile.open("../exec/inputs/case.txt");

if (inputFile.is_open()) {
po::store(po::parse_config_file(inputFile, desc), vm);
inputFile.close();
} else {
std::cerr << "Error opening file: inputs.txt" << std::endl;
}

po::notify(vm);

```

## Class initialisation

Once the input values are stored, the provided initial condition is used to determine the number of particles. It is also utilized to declare containers within the `SPH` object, responsible for storing information related to particle properties, and to allocate memory. This process takes place in the constructor of the class. There, the containers are declared as `new` raw pointers, dynamically allocating memory proportional on the number of particles.

```cpp
// User defined constructor
SPH::SPH(const unsigned n_new) : nb_particles(n_new) {

  position_x = new double[nb_particles];
  position_y = new double[nb_particles];
  velocity_x = new double[nb_particles];
  velocity_y = new double[nb_particles];

  distance = new double[nb_particles * nb_particles];
  distance_q = new double[nb_particles * nb_particles];

  particle_density = new double[nb_particles];

  particle_pressure = new double[nb_particles];

  particle_speed_sq = new double[nb_particles];
}

```

To avoid the use of multiple `if` statements, two `std::map` objects are used to map the different conditions to their corresponding number of particles and their corresponding initialisation function.

```cpp
// Create map to associate initial condition names with number of particles
std::map<std::string, int> initConditionToParticlesMap = {
    {"ic-one-particle", 1},      {"ic-two-particles", 2},
    {"ic-three-particles", 3},   {"ic-four-particles", 4},
    {"ic-dam-break", n3},        {"ic-block-drop", n1 * n2},
    {"ic-droplet", dropletn(n3)}};

nb_particles =
    initConditionToParticlesMap[vm["init_condition"].as<std::string>()];

// Define the solver object (called sph)
// In its definition, the number of particles is required
SPH sph(nb_particles);

// Create map to associate function names with function pointers
std::map<std::string, std::function<void(int, SPH &)>> functionMap = {
    {"ic-one-particle", ic_one_particle},
    {"ic-two-particles", ic_two_particles},
    {"ic-three-particles", ic_three_particles},
    {"ic-four-particles", ic_four_particles},
    {"ic-dam-break", ic_dam_break},
    {"ic-droplet", ic_droplet}};

// Get the function pointer from the map
auto initFunc = functionMap.find(vm["init_condition"].as<std::string>());
if (initFunc != functionMap.end()) {
    int n_particles = nb_particles;

    // The ic-droplet case requires a different n argument.
    if (vm["init_condition"].as<std::string>() == "ic-droplet") {
        n_particles = n3;
    }
    // Retrieves and runs the provided function object
    initFunc->second(n_particles, sph);

} else if (vm["init_condition"].as<std::string>() == "ic-block-drop") {
    /**The ic-block-drop case is not in the map because it has two
    * additional parameters, so it requires a different case.
    **/
    ic_block_drop(nb_particles, n1, n2, sph);
} else {
    std::cerr << "Error: Function not found!" << std::endl;
}

```

## Time integration

Following the initialisation of the class and the output files, the function `time_integration()` is invoked. Within this function, the aforementioned `SPH::` functions are executed at each timestep.

```cpp
for (int t = 0; t < total_iter; t++) {

    sph.calculateParticleDistance();
    sph.calculateDensity();
    sph.particleIterations();
}
```

The total number of timesteps is determined by the total integration time and the prescribed timestep as follows:

```
total_iter =
      ceil(total_time / dt); // Transform time in seconds to iterations
```

## Outputs

Upon successful execution, the program generates two files:

- _Energies File_: This file, containing Total, Kinetic, and Potential energies, is updated at each timestep. The results can be visualized by using the script `post/plot_energies.ipynb`.

- _Particle Positions File_: This file captures the final positions of the particles, and it can be visualized using the script `post/visualize_particles.ipynb`.

```cpp
// Write energies on the Energy-File
vOut2 << t * dt << "  " << sph.getKineticEnergy() << "  "
        << sph.getPotentialEnergy() << "  "
        << sph.getPotentialEnergy() + sph.getKineticEnergy()
        << "\n";

// Get the positions after integration is completed
if (t == total_iter - 1) {

    for (int k = 0; k < nb_particles; k++) {

    vOut << sph.updatePosition_x(k) << " " << sph.updatePosition_y(k)
            << "\n";
    }
}
```
