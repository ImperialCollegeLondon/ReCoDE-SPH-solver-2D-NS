/**In this header file, the SPH class is going to be defined
 * with the appropriate constructors and destructors
 * and the member functions. The SPH class is going to admit a matrix
 * as an input, as well as the number of particles.
 * Details on the functions, the arrays and the overloadings
 * can be found in the file : class.cpp
 **/
#ifndef  PARTICLES_H
#define PARTICLES_H

//#include "sph_calc.h" // TODO(Vyron): Cyclic dependencies here
class particles {

protected:
  unsigned int nb_particles; //number of particles

  // Constants of the problem
  const double gas_constant = 2000.0;
  const double density_resting = 1000.0;
  const double viscosity = 1.0;
  const double acceleration_gravity = 9.81;

  double h; // Radius of influence

  // Positions
  double *position_x;
  double *position_y;

  // Velocities
  double *velocity_x;
  double *velocity_y;
  double *particle_speed_sq;

  // Distances
  double *distance;   // Array to store the distances between the particles
  double *distance_q; // Array to store the values of the normalized distance q

  // Mass
  double mass_assumed = 1.0;

  // Density
  double *particle_density;

  // Pressure
  double *particle_pressure;


  public: 
  /******** CONSTRUCTORS/DESTRUCTOR********/

  particles() = default; // Default constructor

  ~particles(); // Destructor

  particles(const unsigned n_new); // User defined constructor for allocating the
                             // dimensions of the Matrix

  /**********OVERLOADINGS**********/

  double &operator()(unsigned row, unsigned col);

  /**********MEMBER-FUNCTIONS*********/

    // TODO(Vyron): Many of these below need to be inlined and if it's possible constexpr
  // Setter Functions.

  // Function to return the position x
  double get_position_x(int l) const;

  // Function to return the position y
  double get_position_y(int l) const;

  // Function to calculate the kinetic energy
  double get_kinetic_energy();

  // Function to calculate the potential energy
  double get_potential_energy();
};

#endif