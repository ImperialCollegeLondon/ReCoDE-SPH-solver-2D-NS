/**In this header file, the SPH class is going to be defined
 * with the appropriate constructors and destructors
 * and the member functions. The SPH class is going to admit a matrix
 * as an input, as well as the number of particles.
 * Details on the functions, the arrays and the overloadings
 * can be found in the file : class.cpp
 **/
#ifndef SPH_H
#define SPH_H
class SPH {
 private:
  unsigned int nb_particles;  // size of the Matrix

  int t;  // time

  double dt;  // timestep

  double h;  // Radius of influence

  // Constants of the problem
  double gas_constant;
  double density_resting;
  double viscosity;
  double coeff_restitution;
  double acceleration_gravity;

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

  // Mass
  double mass_assumed = 1.0;

  // Density
  double *particle_density;

  // Pressure
  double *particle_pressure;

  // Forces
  double force_pressure, force_viscous, force_gravity;
  double force_pressure_x, force_pressure_y;
  double force_viscous_x, force_viscous_y;
  double force_gravity_y;
  double force_gravity_x = 0.0;

 public:
  /******** CONSTRUCTORS/DESTRUCTOR********/

  SPH() = default;  // Default constructor

  ~SPH();  // Destructor

  SPH(const unsigned n_new);  // User defined constructor for allocating the
                              // dimensions of the Matrix

  /**********OVERLOADINGS**********/

  double &operator()(unsigned row, unsigned col);

  /**********MEMBER-FUNCTIONS*********/

  // Setter Functions.

  // Assign value to dt
  void set_timestep(double dt);

  // Assign value to h
  void set_rad_infl(double h);

  // Assign value to gas_constant
  void set_gas_constant(double gas_constant);

  // Assign value to density_resting
  void set_density_resting(double density_resting);

  // Assign value to viscosity
  void set_viscosity(double viscosity);

  // Assign value to coeff_restitution
  void set_coeff_restitution(double coeff_restitution);

  // Assign value to acceleration_gravity
  void set_acceleration_gravity(double acceleration_gravity);

  // Assign value to gas_constant
  void set_gas_constant(double gas_constant);

  // Assign value to density_resting
  void set_density_resting(double density_resting);

  // Assign value to viscosity
  void set_viscosity(double viscosity);

  // Assign value to coeff_restitution
  void set_coeff_restitution(double coeff_restitution);

  // Assign value to acceleration_gravity
  void set_acceleration_gravity(double acceleration_gravity);

  // Function to calculate the mass of the particles before the simulation starts
  void calc_mass();

  // Function to calculate the matrix with rij
  void calc_particle_distance();

  // Function to calculate the density
  void calc_density();

  // Function to calculate the pressure
  void calc_pressure();

  // Function to perform the particle iterations
  void particle_iterations();

  // Function to calculate the pressure force
  double calc_pressure_force(int particle_index, double *position);

  // Function to calculate the viscous force
  double calc_viscous_force(int particle_index, double *velocity);

  // Function to calculate the gravity force
  double calc_gravity_force(int particle_index);

  // Function to update the positions of the particles
  void update_position(int particle_index);

  // Function to initialise the time integration scheme - velocity
  double scheme_init(int particle_index, double *velocity,
                     double &force_pressure, double &force_viscous,
                     double &force_gravity);

  // Function for time integration - velocity
  double velocity_integration(int particle_index, double *velocity,
                              double &force_pressure, double &force_viscous,
                              double &force_gravity);

  // Function to treat the boundaries
  void boundaries(int particle_index);

  // Function to return the position x
  double return_position_x(int l);

  // Function to return the position y
  double return_position_y(int l);

  // Function to calculate the kinetic energy
  double return_kinetic_energy();

  // Function to calculate the potential energy
  double return_potential_energy();
};

#endif