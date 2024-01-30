# Compiling and executing the SPH-SOLVER

The list of requirements for the source code is:

- A `C++20` version
- The `Boost` library

To compile the code, the user has to change the working directory to `exec/build` and then run the following instructions in the terminal:

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

# Setting up a case

To set up a case the user has to set the parameters of the problem by using the `.txt` files which can be found in the `exec/input` directory. In `case.txt` the user can specify the type of initial condition, the total simulated time (in seconds), the timestep, and the desired output frequency. The initial condition (IC) can be one of the following

- A single particle (`ic-one-particle`) : to test the correctness of time integration and gravity forcing, as well as the bottom boundary condition.

- Two particles (`ic-two-particles`) : to assess the pressure force and viscous terms.

- Three particles (`ic-three-particles`) : to assess left and right boundary conditions.

- Four particles (`ic-four-particles`) : to assess multiple particle interaction.

- A Block drop (`ic-block-drop`): a grid of particles occupying a rectangular region.

- A Droplet (`ic-droplet`): particles occupying a circular region.

After selecting the desired IC, the user has to specify its parameters (number and position of the particles) in the homonymous to the IC `.txt` file. The domain is rectangular and two-dimensional with corners which have coordinates that can be specified in `domain.txt`. Finally, the constant parameters of the problem which characterize the fluid and the solver set-up can be specified in `constants.txt`.

- `case.txt`
  - initial condition (`init_condition`)
  - simulation time (`T`)
  - time-step (`dt`)
  - -output frequency (`output_frequency`)
- `constants.txt`:
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