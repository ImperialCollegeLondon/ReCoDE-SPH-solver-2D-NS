/**In this header file, the SPH class is going to be defined
 * with the appropriate constructors and destructors
 * and the member functions. The SPH class is going to admit a matrix
 * as an input, as well as the number of particles.
 * Details on the functions, the arrays and the overloadings
 * can be found in the file : class.cpp
 **/
#ifndef PARTICLES_H
#define PARTICLES_H
class particles {
 protected:
  unsigned int nb_particles;  // number of particles and characteristic size of the class arrays

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

 public:
  /******** CONSTRUCTORS/DESTRUCTOR********/

  particles() =
      delete;  // Constructor without number of particles shouldn't exist

  ~particles();  // Destructor

  particles(const unsigned n_new);  // User defined constructor for allocating
                                    // the dimensions of the Matrix

  /**********OVERLOADINGS**********/

  double &operator()(unsigned row, unsigned col);

  particles &operator=(const particles &particles);

  /**********MEMBER-FUNCTIONS*********/

  // Function to get the number of particles
  int get_number_of_particles();

  // Function to calculate the matrix with rij
  void calc_particle_distance();

  // Function to get the position x
  double get_position_x(int k);

  // Function to get the position y
  double get_position_y(int k);

  // Function to get the normalised distance between the particles
  double get_distance_q(int k);
};

#endif