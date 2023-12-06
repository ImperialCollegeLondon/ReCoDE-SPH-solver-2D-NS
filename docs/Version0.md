# Smooth particle hydrodynamics (SPH)

# The algorithm

# Description of the code

This version of the code displays a serial implementation of the SPH algorithm in C++. It comprises two .cpp files and their corresponding header (.h) files. 

•main-SPH.cpp
•class.cpp
•class.h
•Makefile

The algorithm is implemented in the class called SPH which also stores information (position, velocity and neighbouring particles) about the particles themselves. The SPH class is initialised in the SPH-main.cpp file where one of the following initial conditions can be specified by the user from the "exec/inputs/case.txt" file:

• A single particle at (0.5, 0.5) to test the correctness of time integration and gravity forcing, as well as the bottom boundary condition.

• Two particles at (0.5, 0.5) and (0.5, h) to assess the pressure force and viscous terms. 

• Three particles at (0.5, 0.5), (0.495, h) and (0.505, h) to assess left and right boundary conditions. 

• Four particles at: (0.505, 0.5), (0.515, 0.5), (0.51, 0.45) and (0.5, 0.45) to assess multiple particle interaction.

• Dam break: a grid of particles occupying the region $[0,0.2]^2$.

• Block drop: a grid of particles occupying the region $[0.1,0.3]\times[0.3,0.6]$.

• Droplet: particles occupying a circle of radius 0.1, centred at the point $[0.5,0.7]$

In the same file the user can specify the time of integration (T) as well as the timestep (dt).

The time integration loop is being executed in the main program where the SPH functions are being called and the output files are being written. Uppon succesful execution the program will result in two files:

• One for the Energies (Total, Kinetic, Potential) at each time step.

• One for the final positions of the particles.

# Structure of the class

The SPH class is initalised by using the number of particles (N) which is required to determine the size of the arrays in the constructor. Several operators have been overloaded so that the corresponding variables can be set inside the class. Specifically the operator "()" has been overloaded to place the particles in their initial conditions and to set their initial velocities, which are stored in a $4\times N$ matrix. The operator ">" to set the time of integration, the operator ">>" to set the time-step and "<" to set the radius of influence. The class has three main functions for the temporal integration:

• rVec: Calculates the distances between the particels .

• den: Updates the density.

• spatial: Calculates all the forces between the particles and updates their positions based on the leapfrog scheme.
