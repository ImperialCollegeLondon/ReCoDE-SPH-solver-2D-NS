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
  unsigned int nb_particles; // size of the Matrix

  int t; // time

  double dt; // timestep

  double h; // Radius of influence

  // Constants of the problem
  const double gas_constant = 2000.0;
  const double density_resting = 1000.0;
  const double viscosity = 1.0;
  const double coeff_restitution = 0.5;
  const double acceleration_gravity = 9.81;

  // Positions
  double *position_x;
  double *position_y;

  // Velocities
  double *velocity_x;
  double *velocity_y;
  double *particle_speed;

  // Distances
  double *distance;   // Array to store the distances between the particles
  double *distance_q; // Array to store the values of the normalised distance q

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

  SPH() = default; // Default constructor

  ~SPH(); // Destructor

  SPH(const unsigned n_new); // User defined constructor for allocating the
                             // dimensions of the Matrix

  /**********OVERLOADINGS**********/

  double &operator()(unsigned row, unsigned col);

  /**********MEMBER-FUNCTIONS*********/

  // Setter Functions.
  void set_time(double t);

  void set_timestep(double dt);

  void set_rad_infl(double h);

  void calc_particle_distance();

  void calc_density();

  void calc_pressure();

  double calc_pressure_force(int index_i, double *x_y);

  double calc_viscous_force(int index_i, double *velocity);

  double calc_gravity_force(int index_i);

  double scheme_init(int index_i, double *velocity, double &force_pressure,
                     double &force_viscous, double &force_gravity);

  double velocity_integration(int index_i, double *velocity,
                              double &force_pressure, double &force_viscous,
                              double &force_gravity);

  void spatial();

  double return_position_x(int l);

  double return_position_y(int l);

  double return_kinetic_energy();

  double return_potential_energy();

  void calc_mass();
};

#endif