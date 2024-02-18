# Code overview

The code herein contains a serial C++ implementation of the SPH methodology described in `SPH.md`. The variables associated with the particles' positions, velocities and forces, are stored as members of an object called `fluid`, while the methods of another class called `SphSolver` manifest the steps of the algorithm.

## Files

The `src` directory comprises five `*.cpp` files and their corresponding header (`*.h`) files, as well as the files to build the code.

- `src/SPH-main.cpp`
- `src/sph.{h, cpp}`
- `src/initial_conditions.{h, cpp}`
- `src/main_prog_func.h`
- `src/particles.{h, cpp}`
- `src/fluid.{h, cpp}`
- `src/sph_solver.{h, cpp}`
- `src/CMakeLists.txt`

## Reading Inputs

The program expects the parameters which are specified in the `exec/input` directory. When the function `initialize()` is invoked by the main program to read the `.txt` files, the `<boost/program_options.hpp>` library is utilised. This library facilitates the mapping of these parameters to their corresponding values, which are then stored in their respective variable maps in the `retrieveInputsFromFile()` function. This approach enhances the flexibility and robustness of the input reading process. Users can specify input parameters in the `.txt` files in any order, provided they are presented as `key = value` pairs.

```cpp
/* **************************** SPH_main.cpp **************************** */

void initialise(fluid** fluidPtr, SphSolver& sphSolver) {

  // Process to obtain the inputs provided by the user
  po::options_description desc("Allowed options");
  desc.add_options()("init_condition", po::value<std::string>(),
                      "take an initial condition")("T", po::value<double>(),
                                                  "take integration time")(
      "dt", po::value<double>(), "take time-step")("h", po::value<double>(),
                                                  "take radius of influence")(
      "gas_constant", po::value<double>(), "take gas constant")(
      "density_resting", po::value<double>(), "take resting density")(
      "viscosity", po::value<double>(), "take viscosity")(
      "acceleration_gravity", po::value<double>(), "take acc due to gravity")(
      "coeff_restitution", po::value<double>(), "take coeff of restitution")(
      "left_wall", po::value<double>(), "take left wall position")(
      "right_wall", po::value<double>(), "take right wall position")(
      "bottom_wall", po::value<double>(), "take bottom wall position")(
      "top_wall", po::value<double>(), "take top wall position")(
      "length", po::value<double>(), "take length of the block")(
      "width", po::value<double>(), "take width of the block")(
      "radius", po::value<double>(), "take radius of the droplet")(
      "n", po::value<int>(), "take number of particles")(
      "center_x", po::value<double>(), "take center of the particle mass in x")(
      "center_y", po::value<double>(), "take center of the particle mass in y")(
      "init_x_1", po::value<double>(), "take x_1")(
      "init_y_1", po::value<double>(), "take y_2")(
      "init_x_2", po::value<double>(), "take x_2")(
      "init_y_2", po::value<double>(), "take y_2")(
      "init_x_3", po::value<double>(), "take x_3")(
      "init_y_3", po::value<double>(), "take y_3")(
      "init_x_4", po::value<double>(), "take x_4")(
      "init_y_4", po::value<double>(), "take y_4")(
      "outputFrequency", po::value<int>(),
      "take frequency that output will be written to file");
  ...

  icCase = caseVm["init_condition"].as<std::string>();
  std::string fileName = icCase + ".txt";
  retrieveInputsFromFile(fileName, icCase, desc, icVm);

  handleInputErrors(caseVm, domainVm, constantsVm, icVm);
  ...
}


void retrieveInputsFromFile(std::string fileName, std::string icCase,
                            po::options_description desc,
                            po::variables_map& vm) {
  std::ifstream caseFile;
  std::string errorMessage = "Error opening file: " + fileName;
  if (fileName == icCase + ".txt") {
    errorMessage +=
        " Make sure that the value of the init_condition in the case.txt "
        "file is one of the following: ic-one-particle, ic-two-particles, "
        "ic-three-particles, ic-four-particles, ic-droplet, ic-block-drop.";
  }
  // Try to open the file
  try {
    caseFile.open("../input/" + fileName);
    // Throw an exception if the file cannot be opened
    if (!caseFile.is_open()) {
      throw std::runtime_error(errorMessage);
    }
    po::store(po::parse_config_file(caseFile, desc), vm);
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  po::notify(vm);
}
```

Additionally, using the `handleInputErrors()` function, error handling is integrated into the input file reading process to guarantee that the provided values conform to the constraints imposed by the underlying mathematical models and the physical meaning of each variable. For example, if the user attempts to set a negative value for the timestep, or a value that is greater than the integration time, the program will throw an error and instruct the user to choose a more suitable value.

```cpp
/* **************************** SPH_main.cpp **************************** */

void handleInputErrors(po::variables_map caseVm, po::variables_map domainVm,
                       po::variables_map constantsVm, po::variables_map icVm) {
  try {
    // Error handling for the total integration time
    if (caseVm["T"].as<double>() <= 0) {
      throw std::runtime_error(
          "Error: Total integration time must be positive!");
      // Error handling for the time step
    } else if (caseVm["dt"].as<double>() <= 0 or
               caseVm["dt"].as<double>() > caseVm["T"].as<double>()) {
      throw std::runtime_error(
          "Error: Time step must be positive and lower than the total "
          "integration time!");
      // Error handling for the output frequency
    } else if (caseVm["output_frequency"].as<int>() <= 0 or
               caseVm["output_frequency"].as<int>() >
                   ceil(caseVm["T"].as<double>() / caseVm["dt"].as<double>())) {
      throw std::runtime_error(
          "Error: Output frequency must be positive and lower than the total "
          "number of iterations!");
      // Error handling for the domain boundaries input
    } else if (domainVm["left_wall"].as<double>() >=
                   domainVm["right_wall"].as<double>() ||
               domainVm["bottom_wall"].as<double>() >=
                   domainVm["top_wall"].as<double>()) {
      throw std::runtime_error(
          "Error: Please adjust your domain boundaries so that left_wall < "
          "right wall and bottom_wall < top_wall.");
      // Error handling for the number of particles
    } else if (icVm["n"].as<int>() <= 0) {
      throw std::runtime_error("Error: Number of particles must be positive!");
    }
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }
}
```

## Class initialisation

The code makes use of three different classes which are purposed to represent the fluid and the SPH algorithm deployed in this project. More details regarding the classes and the design choices can be found in the `docs/OOP_concepts.md` and the reader is advised to study it before proceeding with this chapter.

Firstly, one SphSolver object and one fluid pointer to an object are being declared in the main program. The pointer declaration is used for the `fluid`, because to initialise the object properly the number of particles is required in the user defined constructor and this information is not yet available since the input files have not been read. These objects are passed as a reference to the `initialise()` function. 

Once the input values are read and stored, the provided IC is used to determine the number of particles. This means that although the user has already provided a number of particles, this is just an indication, since the IC (droplet and block drop) require specific formation and the particles to be distributed uniformly. These two conditions cannot be satisfied simultaneously by any number of particles and therefore several adjustments need to be made. The functions `closestIntegerSqrt()` and `rectangleN()` from `initial_conditions.h` are functions suitable for this purpose. 

The IC functions are being called within the `setInitialConditions()` function and a reference to the pointer of the fluid object is passed as an argument<b>*</b>, as well as the updated number of particles. Inside these functions the user defined constructor of the `fluid` is being called and the memory allocation process for the object's containers is invoked. In this, the containers are declared as `new` raw pointers to arrays, dynamically allocating memory proportional to the number of particles. The function used to initialise the `fluid` class for the simple cases of 1,2,3 and 4 particles is demonstrated below.

<b>\*</b> The rationale behind passing the pointer to the fluid object as a reference is identical to the reason an object is passed as a reference to a function. In these functions we are allocating new memory that the pointer should point to. If the pointer was passed by value (`Fluid *fluidPtr`) then a **copy** of the pointer would be pointing to the new memory and our original fluid pointer will still be a `nullptr`. Of course, since the allocation is happening inside the function, the caller is responsible to `delete` the object manually.

```cpp
/* **************************** initial_conditions.cpp **************************** */

void icBasic(std::unique_ptr<Fluid> &fluidPtr, int nbParticles,
             std::vector<double> &positionX, std::vector<double> &positionY) {
  fluidPtr = std::make_unique<Fluid>(nbParticles);

  Fluid &fluid = *fluidPtr;  // Use a reference to the object

  for (int i = 0; i < nbParticles; i++) {
    fluid.setPositionX(i, positionX[i]);
    fluid.setPositionY(i, positionY[i]);
    fluid.setVelocityX(i, 0.0);
    fluid.setVelocityY(i, 0.0);
  }
}
```

Inside the `setInitialConditions()` function, an `std::map` object is used to map the initial condition case with a fixed number of particles to their corresponding number of particles, to avoid the use of multiple `if` statements and make the code cleaner. Then, the initial condition case is set based on the user's input, including additional error handling for the ic-case specific input parameters. The workflow for the droplet case is demonstrated below.

```cpp
/* **************************** SPH_main.cpp **************************** */

void setInitialConditions(std::string icCase, std::unique_ptr<Fluid>& fluidPtr,
                          po::variables_map icVm, po::variables_map domainVm) {
  // Fixed nbParticles ic cases map
  std::map<std::string, int> initConditionToParticlesMap = {
      {"ic-one-particle", 1},
      {"ic-two-particles", 2},
      {"ic-three-particles", 3},
      {"ic-four-particles", 4}};

  // Get the number of particles based on the ic case (for the more complex ic)
  int nbParticles;
  if (icCase == "ic-droplet" || icCase == "ic-block-drop") {
    nbParticles = icVm["n"].as<int>();
  } else {
    nbParticles = initConditionToParticlesMap[icCase];
  }

  // Block of code ...

  } else if (icCase == "ic-droplet") {
    // Get the droplet radius and center coordinates from the ic file
    double radius = icVm["radius"].as<double>();
    // Error handling for the droplet radius
    try {
      if (radius <= 0) {
        throw std::runtime_error("Error: Radius must be positive!");
      }
    } catch (std::runtime_error& e) {
      // Handle the exception by printing the error message and exiting the
      // program
      std::cerr << e.what() << std::endl;
      exit(1);
    }
    double centerX = icVm["center_x"].as<double>();
    double centerY = icVm["center_y"].as<double>();
    // Error handling for the droplet initial position (centerX, centerY)
    try {
      if (centerX - radius < domainVm["left_wall"].as<double>() ||
          centerX + radius > domainVm["right_wall"].as<double>() ||
          centerY - radius < domainVm["bottom_wall"].as<double>() ||
          centerY + radius > domainVm["top_wall"].as<double>()) {
        throw std::runtime_error(
            "Error: The droplet must be within the domain boundaries! Please "
            "adjust the center coordinates.");
      }
    } catch (std::runtime_error& e) {
      // Handle the exception by printing the error message and exiting the
      // program
      std::cerr << e.what() << std::endl;
      exit(1);
    }
    icDroplet(fluidPtr, nbParticles, radius, centerX, centerY);
  } else {
    std::cerr << "Error: Initial condition function not found! Make sure "
              << "that the value of the init_condition in the case.txt file is "
              << "one of the following: ic-one-particle, ic-two-particles, "
              << "ic-three-particles, ic-four-particles, ic-droplet, "
              << "ic-block-drop." << std::endl;
    exit(1);
  }
}

```

Finally, after the object initialisation, the rest of the parameters which are required by the `SphSolver` and the `Fluid` objects are being set with the use of setter functions, the Neighbour Searching grid is created, the initialised particles are placed in their corresponding grid cells, and their neighbours are identified.

```cpp
  /* **************************** SPH_main.cpp **************************** */

  // void initialise(fluid** fluidPtr, SphSolver& sphSolver) { ...

  // Set the parameters of the solver for the specific simulation
  sphSolver.setTimestep(caseVm["dt"].as<double>());
  sphSolver.setTotalIterations(
      ceil(caseVm["T"].as<double>() / caseVm["dt"].as<double>()));
  sphSolver.setOutputFrequency(caseVm["output_frequency"].as<int>());
  sphSolver.setCoeffRestitution(constantsVm["coeff_restitution"].as<double>());
  sphSolver.setLeftWall(domainVm["left_wall"].as<double>());
  sphSolver.setRightWall(domainVm["right_wall"].as<double>());
  sphSolver.setTopWall(domainVm["top_wall"].as<double>());
  sphSolver.setBottomWall(domainVm["bottom_wall"].as<double>());
  sphSolver.setPrecalculatedValues(constantsVm["h"].as<double>());

  // Define the fluid based on the inputs
  fluidPtr->setRadInfl(constantsVm["h"].as<double>());
  fluidPtr->setGasConstant(constantsVm["gas_constant"].as<double>());
  fluidPtr->setDensityResting(constantsVm["density_resting"].as<double>());
  fluidPtr->setViscosity(constantsVm["viscosity"].as<double>());
  fluidPtr->setAccelerationGravity(
      constantsVm["acceleration_gravity"].as<double>());

  // Calculate the mass of the particles
  sphSolver.createGrid(*fluidPtr);
  sphSolver.neighbourParticlesSearch(*fluidPtr);
  std::vector<std::vector<std::pair<int, double>>> neighboursPerParticle =
      sphSolver.getNeighbourParticles();
  fluidPtr->calculateMass(neighboursPerParticle);
```

## Output files

The output files are being initialised with the use of the `initOutputFiles()`. The outputs are exported in `.csv` format which displays good readability and facilitates data manipulation compared to `.txt` files. They are stored in a centralised location, specifically within the `/exec/output/` directory. This centralisation simplifies data organisation and retrieval, making it easier for users to access and analyse output data.

Upon successful execution, the program generates two types of files:

- _Energies File_: This file, containing Total, Kinetic, and Potential energies, is updated at each timestep. The results can be visualized by using the script `post/plot_energies.py`.

- _Particle Positions File_: This file captures the positions of the particles at a given timestep, and it can be visualized using the script `post/visualize_particles.py`.

```cpp
/* **************************** SPH_main.cpp **************************** */

void storeToFile(fluid& fluid, int nbParticles, std::string type,
                 std::ofstream& targetFile, double dt, int currentIteration) {
  if (type == "energy") {
    // Write energies on the Energy-File
    targetFile << currentIteration * dt << "," << fluid.getKineticEnergy()
               << "," << fluid.getPotentialEnergy() << ","
               << fluid.getPotentialEnergy() + fluid.getKineticEnergy()
               << "\n";
  } else if (type == "position") {
    for (int k = 0; k < nbParticles; k++) {
      targetFile << fluid.getPositionX(k) << "," << fluid.getPositionY(k)
                 << "\n";
    }
  }
}
```

<div style="text-align: center;">
    <img src="images/energies_README.png" alt="Alt Text 1" style="display: inline-block; width: 400px;">
    <img src="images/positions_README.png" alt="Alt Text 2" style="display: inline-block; width: 400px;">
</div>

 ***Energy plots (left) and initial position (right) for a droplet of 60 particles.***

## Time integration

Following the initialisation of the class and the output files, the function `SphSolver::timeIntegration()` is invoked. Within this function, the steps of the SPH algorithm are executed, and the outputs are exported.

```cpp
  /* ***************************** SPH-main.cpp ****************************** */

  // Time integration loop
  SphSolver.timeIntegration(*sphFluid, finalPositionsFile, energiesFile);

  /* **************************** sph_solver.cpp **************************** */

  void SphSolver::timeIntegration(fluid &data,
                                  std::ofstream &finalPositionsFile,
                                  std::ofstream &energiesFile) {
  std ::cout << "Time integration started -- OK"
             << "\n";

  numberOfParticles = data.getNumberOfParticles();

  for (int time = 0; time < totalIterations; time++) {
    t = time;
    // In each iteration the distances between the particles are recalculated,
    // as well as their density and pressure
    data.calculateParticleDistance();
    data.calculateDensity();
    data.calculatePressure();
    particleIterations(data);

    if (time % outputFrequency == 0) {
      storeToFile(data, "energy", energiesFile, dt, t);
    }
  }
  // Store particles' positions after integration is completed
  storeToFile(data, "position", finalPositionsFile, dt, totalIterations);

  std ::cout << "Time integration finished -- OK"
             << "\n";
}

```
