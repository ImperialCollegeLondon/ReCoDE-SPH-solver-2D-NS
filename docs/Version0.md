# Description of the code

## Overview

The code in `v0` contains a serial C++ implementation of the algorithm described in `SPH.md`. The variables which are associated with the particles' positions, velocities and forces, are stored as members of an object called `SPH`. Furthermore, the functions which manifest the steps of the aforementioned algorithm are declared as the methods of this `SPH` class. It must be noted that although the produced results are correct, this version will be improved upon in chapters which follow. At each step, the improvements will be explained in terms of how they work and the motivation for their implementation.

## Compiling and executing the code

The list of requirements for the v0 code is:

- A `C++17` version 
- The `boost_program_options` library

To compile the user has to simply type the following commands in the terminal:

- `make clean`
- `make`

This will produce an executable called `SPH-SOLVER` in the src folder and the user needs to type:

- `./SPH-SOLVER`

## Files
This version of the code displays a serial implementation of the SPH algorithm in C++. It comprises three `*.cpp` files and their corresponding header (`*.h`) files. The code is accompanied by two input text files, one for the input variables of the executed case and the input parameters and one for the domain boundaries. 

- `src/main-SPH.cpp`
- `src/sph.cpp`
- `src/initial_conditions.cpp`
- `src/main_prog_func.h`
- `src/sph.h`
- `src/initial_conditions.h`
- `src/Makefile`
- `exec/inputs/case.txt`

## Main program

The main program of the code (i.e.`main-SPH.cpp`) is used to read the input files, export the output files, initialise the SPH class and perform the time integration by calling the SPH object's methods.

## Structure of the class

The SPH class is initialised by using the number of particles (`N`) which is required to determine the size of the arrays in the constructor. The operator `()` has been overloaded to place the particles in their initial conditions and to set their initial velocities to the corresponding arrays. The class has three main functions for the temporal integration:

- `SPH::calc_particle_distance()`: Calculates the distances between the particles.

- `SPH::calc_density()`: Updates the density.

- `SPH::particle_iterations()`: Calculates all the forces acting on the particles and updates their positions based on the leapfrog scheme. 

## Input Parameters

### Initial conditions
As stated earlier, the SPH class is initialised in the SPH-main.cpp file where one of the following initial conditions can be specified by the user from the `exec/inputs/case.txt` file:

- A single particle (`ic-one-particle`) at : (0.5, 0.5) to test the correctness of time integration and gravity forcing, as well as the bottom boundary condition.

- Two particles (`ic-two-particles`) at : (0.5, 0.5) and (0.5, h) to assess the pressure force and viscous terms. 

- Three particles (`ic-three-particles`) at : (0.5, 0.5), (0.495, h) and (0.505, h) to assess left and right boundary conditions. 

-  Four particles (`ic-four-particles`) at : (0.505, 0.5), (0.515, 0.5), (0.51, 0.45) and (0.5, 0.45) to assess multiple particle interaction.

-  A Dam break (`ic-dam-break`): a grid of particles occupying the region $[0,0.2]^2$.

-  A Block drop (`ic-block-drop`): a grid of particles occupying the region $[0.1,0.3]\times[0.3,0.6]$.

-  A Droplet (`ic-droplet`): particles occupying a circle of radius 0.1, centered at the point $[0.5,0.7]$


### Simulation Duration
The real simulated time can be set in units of seconds with the input key `T`. For instance, a simulation of 10 seconds can be set by typing `T = 10` in the `exec/inputs/case.txt` file.

### Time-step

The time-step `dt` can be set in units of seconds with the input key `dt`. For instance, a time-step of 0.0001 seconds can be set by typing `dt = 1e-4` in the `exec/inputs/case.txt` file.

### Radius of Influence
The last input parameter that the user can specify in the `exec/inputs/case.txt` file is the radius of influence `h` in units of m. This parameter dictates the maximum distance in which an entity (particle or domain boundary) has an influence in the behavior of another entity.

### Reading Inputs
The aforementioned parameters are expected by the program, and therefore, while reading the `exec/inputs/case.txt` file in the function `initialise()` which is called by the main program, the `<boost/program_options.hpp>` library is used to map those parameters to their values, which are finally stored in their corresponding variables. This practice constitutes in making the input reading process more flexible and error-proof. The user can specify the input parameters in the `exec/inputs/case.txt` file in any order, as long as they are given as `key = value` pairs.

``` 
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
After storing the input values, the initial condition is used to determine the number of particles, as well as to declare the containers which store the information related to the particles' properties in the `sph` object and allocate memory. This is done in the constructor of the class where the containers are declared as `new` raw pointers and occupy memory that depends on the number of particles. 

```
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

```
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
initFunc->second(n_particles, sph);

} else {
/**The ic-block-drop case is not in the map because it has two
    * additional parameters, so it requires a different case.
    **/
if (vm["init_condition"].as<std::string>() == "ic-block-drop") {
    ic_block_drop(nb_particles, n1, n2, sph);

} else {
    std::cerr << "Error: Function not found!" << std::endl;
}
}

```

## Time integration

After initializing the class and the output files, the function `time_integration()` is called where the aforementioned `SPH::` functions are being executed during every timestep.

```
for (int t = 0; t < total_iter; t++) {

    sph.calc_particle_distance();
    sph.calc_density();
    sph.particle_iterations();
}
```
The total number of timesteps is determined by the total integration time and the prescribed timestep as follows:

```
total_iter =
      ceil(total_time / dt); // Transform time in seconds to iterations
```
## Outputs

Uppon succesful execution the program will result in two files:

- One for the Energies (Total, Kinetic, Potential) which is written during every time step and the results be plotted by using the script `post/plot_energies.ipynb`.

- One for the final positions of the particles which can be plotted by using the script `post/visualise_particles.ipynb`.

```
// Write energies on the Energy-File
vOut2 << t * dt << "  " << sph.return_kinetic_energy() << "  "
        << sph.return_potential_energy() << "  "
        << sph.return_potential_energy() + sph.return_kinetic_energy()
        << "\n";

// Get the positions after integration is completed
if (t == total_iter - 1) {

    for (int l = 0; l < nb_particles; l++) {

    vOut << sph.return_position_x(l) << " " << sph.return_position_y(l)
            << "\n";
    }
}
```


