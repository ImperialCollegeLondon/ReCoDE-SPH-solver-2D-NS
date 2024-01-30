# Object Oriented Programming

Object Oriented Programming (OOP) is a style of programming adopted by many programming languages (e.g. C++,Java, Python) and allows the developers to define a user defined data type known as a class. A user may create instances of these classes, known as objects, and use them to store and manipulate data. The contents of the class will define the data that is stored in instances of the class in member variables, and functions known as methods which can access, manipulate, and perform calculations on the data.

The complexity of a class can be very simple, such as a variable type (i.e. integer, double, array etc.) or very complicated, such as a solver for complex non-linear problems. Every class comprises several components which are the members, the methods, the constructors and the destructor.

## Inheritance

Depending on the application, a developer can create families of objects which usually have conceptual dependency and coherence among them. A class-family consists of a base class and its derived classes, which may further have their own subclasses. The base class serves as a representation of a general concept, such as a "human" with attributes like name, age, and height. The derived classes, in turn, embody more specific instances of that concept. For example, a child class might represent a "student" with additional attributes like grades and the number of extracurricular activities. These supplementary details, which may be irrelevant when characterizing another type of human, such as a retired professional, are encapsulated as member variables within the student class. The process of deriving a new child class based on a base class is called inheritance.

## Object orientation in the SPH solver

In this project we chose object orientation which is widely supported by C++ in order to facilitate the implementation of the SPH algorithm and exploit all the benefits that accompany the use of OOP. The basic design idea is built around three classes.

## Class-particles
The first class is used to represent the building blocks of the SPH approach which are the particles. A particle however is not only relevant in the context of SPH, but it can also be used in other applications such as in the representation of multiphase flows where the particles can be droplets whose motion is two-way coupled with a carrier gas phase. Therefore, the "particles" class was decided to be kept as simple (and generic) as possible and to encapsulate only the information which can be applicable to every application that incorporates the use of cluster of particles. Therefore, the members of the "particles" class are only related to the particles' positions and velocities.

It is important to emphasize that an instance of the "particles" class contains data pertaining to all the particles involved in the simulation (therefore represents the entirety of the cluster), rather than information exclusive to a single particle.

```cpp
class particles {
 protected:
  unsigned int nb_particles;  // size of the Matrix

  // Positions
  double *position_x;
  double *position_y;

  // Velocities
  double *velocity_x;
  double *velocity_y;

  double *particle_speed_sq;

  // Distances
  double *distance;    // Array to store the distances between the particles
  double *distance_q;  // Array to store the values of the normalised distance q

  ...


};
```

The particles class is initialised in the `user defined constructor` by using the number of particles (`nb_particles`) which is required to determine the size of the arrays. 

```cpp
// User defined constructor
particles::particles(const unsigned n_new) : nb_particles(n_new) {
  position_x = new double[nb_particles];
  position_y = new double[nb_particles];
  velocity_x = new double[nb_particles];
  velocity_y = new double[nb_particles];

  distance = new double[nb_particles * nb_particles];
  distance_q = new double[nb_particles * nb_particles];

  particle_speed_sq = new double[nb_particles];
}
```

The operator `()` has been overloaded to place the particles in their initial conditions and to set their initial velocities to the corresponding arrays.

```cpp
// Overloading of ()
double &particles::operator()(unsigned row, unsigned col) {
  switch (row) {
    case 0:
      return this->position_x[col];
      break;
    case 1:
      return this->position_y[col];
      break;
    case 2:
      return this->velocity_x[col];
      break;
    case 3:
      return this->velocity_y[col];
      break;
    default:
      std::cerr << "ERROR: Out of bounds on row selection" << std::endl;
      abort();
  }
}
```

## Class-sph-fluid
A cluster of particles, depending on its behavior can represent a variety of concepts. One of these concepts is what we will refer to herein as an SPH-fluid, which apart from being discretized in particles also has other attributes such as density, mass and pressure. Therefore, the class sph-fluid is a child class of the particles class which is extended in order to encapsulate more member variables and methods for their estimation in order to fully describe the simulated fluid. In the main program the base class will never be invoked explicitly, but only implicitly through the instantiation of a sph-fluid object. However, from the developer's perspective, by using this approach the code becomes more modular and allows for the derivation of multiple children from the base class if needed.

```cpp
class fluid : public particles {
 private:
  // Constants of the problem
  double gas_constant;
  double density_resting;
  double viscosity;
  double acceleration_gravity;

  double h;  // Radius of influence

  // Mass
  double mass_assumed = 1.0;

  // Density
  double *density;

  // Pressure
  double *pressure;

 ...
 
 }
```


## Class-sph_solver
In the `sph_solver` class, are implemented the steps of the algorithm described in `SPH.md`. The main function which is called by the main program is the `sph_solver::time_integration()`