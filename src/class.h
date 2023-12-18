/**In this header file, the SPH class is going to be defined
 * with the appropriate constructors and destructors
 * and the member functions. The SPH class is going to admit a matrix
 * as an input, as well as the number of particles.
 * Details on the functions, the arrays and the overloadings
 * can be found in the file : class.cpp
 **/
#ifndef CLASS_H
#define CLASS_H
class SPH {

private:
  unsigned int nb_particles; // size of the Matrix

  int t; // time

  double dt; // timestep

  double h; // Radius of influence

  double **vMatrix = new double *[4]; // AlloCating space memory for the Matrix

  // Constants of the problem
  double gas_constant = 2000.0;
  double density_resting = 1000.0;
  double mass_assumed = 1.0;
  double viscosity = 1.0;
  double coeff_restitution = 0.5;
  const double acceleration_gravity = 9.81;

public:
  double *position_x;
  double *position_y;
  double *velocity_x;
  double *velocity_y;
  double *distance;
  double *particle_density;
  double *particle_pressure;
  double *particle_speed;
  double *q;
  double force_pressure, force_viscous, force_gravity;
  double force_pressure_x, force_pressure_y;
  double force_viscous_x, force_viscous_y;
  double force_gravity_y;
  double force_gravity_x = 0.0;
  int i, j;

  /******** CONSTRUCTORS/DESTRUCTOR********/

  SPH() = default; // Default constructor

  ~SPH(); // Destructor

  SPH(const unsigned n_new); // User defined constructor for allocating the
                             // dimensions of the Matrix

  /**********OVERLOADINGS**********/

  double &operator()(unsigned row, unsigned col);

  int &operator>(unsigned t);

  double &operator>>(double dt);

  double &operator<(double h_new);

  /**********MEMBER-FUNCTIONS*********/

  void position_x_init();

  void position_y_init();

  void velocity_x_init();

  void velocity_y_init();

  void calc_particle_distance();

  void calc_density();

  void calc_pressure();

  double calc_pressure_force(double *x_y);

  double calc_viscous_force(double *velocity);

  double calc_gravity_force();

  double scheme_init(double *velocity, double &force_pressure,
                     double &force_viscous, double &force_gravity);

  double velocity_integration(double *velocity, double &force_pressure,
                              double &force_viscous, double &force_gravity);

  void spatial();

  double return_position_x(int l);

  double return_position_y(int l);

  double return_kinetic_energy();

  double return_potential_energy();

  void calc_mass();
};

#endif