#include "class.h"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>

// Defining the user defined constructor
SPH::SPH(const unsigned n_new) {

  nb_particles = n_new;

  for (int i = 0; i < 4; i++) {

    vMatrix[i] =
        new double[nb_particles]; // Array that will be used to
                                  // pass and store the values of the initial
                                  // positions and velocities inside the class
  }

  position_x = new double[nb_particles]; // Array to store the positions in x of
                                         // the particles
  position_y = new double[nb_particles]; // Array to store the positions in y of
                                         // the particles
  velocity_x = new double[nb_particles]; // Array to store the values of the
                                         // velocity in x of the particles
  velocity_y = new double[nb_particles]; // Array to store the values of the
                                         // velocity in y of the particles
  distance = new double[nb_particles *
                        nb_particles]; // Array to store the values of the
                                       // distances between the particles
  q = new double[nb_particles * nb_particles]; // Array to store the values of q
  particle_density = new double[nb_particles]; // Array to store the densitites
                                               // of the particles
  particle_pressure =
      new double[nb_particles]; // Array to store the pressure on each particle
  particle_speed =
      new double[nb_particles]; // Array to store the norm of the velocity, used
                                // in the calculation of the kinetic energy
}

// Destructor
SPH::~SPH() {
  delete[] vMatrix;
  delete[] particle_density;
  delete[] particle_pressure;
  delete[] distance;
  delete[] particle_speed;
  delete[] position_x;
  delete[] position_y;
  delete[] velocity_x;
  delete[] velocity_y;
}

// Overloading of ()
double &SPH::operator()(unsigned row, unsigned col) {

  return vMatrix[row][col];
}

// Overloading of > to pass the number of the current time
// iteratiion to the class
int &SPH::operator>(unsigned t_new) {

  t = t_new;
  return t;
}

// Overloading of >> to pass the time-step inside tha class
double &SPH::operator>>(double dt_new) {

  dt = dt_new;
  return dt;
}

// Overloading of < to pass the radius of indfluence inside the
// class
double &SPH::operator<(double h_new) {

  h = h_new;
  return h;
}

/**Assig values to x: The values of the positions on the
 * x-axis are stored in the first row of the vMatrix, and they are now assigned
 * to a single vector **/
void SPH::position_x_init() {

  for (int i = 0; i < nb_particles; i++) {

    position_x[i] = vMatrix[0][i];
  }
}

/**Assign values to y: The values of the positions on the
 * y-axis are stored int he second row of the vMatrix, and they are now assigned
 * to a single vector **/
void SPH::position_y_init() {

  for (int i = 0; i < nb_particles; i++) {

    position_y[i] = vMatrix[1][i];
  }
}

/**Assign values to vx: The values of the initial velocities in the x
 * direction are stored int the third row of the vMatrix, and they are now
 * assigned to a single vector **/
void SPH::velocity_x_init() {

  for (int i = 0; i < nb_particles; i++) {

    velocity_x[i] = vMatrix[2][i];
  }
}

/**Assign values to vy: The values of the initial velocities in the x
 * direction are stored int the fourth row of the vMatrix, and they are now
 * assigned to a single vector **/
void SPH::velocity_y_init() {

  for (int i = 0; i < nb_particles; i++) {

    velocity_y[i] = vMatrix[3][i];
  }
}

// Function to calculate the matrix with rij
void SPH::calc_particle_distance() {

  double dx;
  double dy;

  for (int i = 0; i < nb_particles; i++) {

    for (int j = 0; j < nb_particles; j++) {

      dx = position_x[i] - position_x[j];
      dy = position_y[i] - position_y[j];

      distance[i * nb_particles + j] = sqrt(dx * dx + dy * dy);
    }
  }
}

// Function to calculate the density
void SPH::calc_density() {

  double phi;
  double pre = (4.0 / (M_PI * h * h)); // Precalculated value
  double hinv = 1.0 / h;               // This is to avoid many divisions

  // find φ
  for (int i = 0; i < nb_particles; i++) {

    particle_density[i] = 0;

    for (int j = 0; j < nb_particles; j++) {

      q[i * nb_particles + j] = std::abs(distance[i * nb_particles + j] * hinv);

      if (q[i * nb_particles + j] < 1) {

        phi = pre * (1.0 - q[i * nb_particles + j] * q[i * nb_particles + j]) *
              (1.0 - q[i * nb_particles + j] * q[i * nb_particles + j]) *
              (1.0 - q[i * nb_particles + j] * q[i * nb_particles + j]);

      }

      else {

        phi = 0.0;
      }

      particle_density[i] += mass_assumed * phi;
    }
  }
}

// Function for the pressure
void SPH::calc_pressure() {

  for (int i = 0; i < nb_particles; i++) {

    particle_pressure[i] =
        gas_constant * (particle_density[i] - density_resting);
  }
}

double SPH::calc_pressure_force(double *x_y) {

  double sum = 0.0;                          // Initializing the sumation
  double pre = (-30.0 / (M_PI * h * h * h)); // precalculated value

  for (int j = 0; j < nb_particles; j++) {

    if (i == j) {
    } else {

      if (q[i * nb_particles + j] < 1) {

        sum += (mass_assumed / particle_density[j]) *
               ((particle_pressure[i] + particle_pressure[j]) / 2.0) *
               (pre * (x_y[i] - x_y[j])) *
               (((1.0 - q[i * nb_particles + j]) *
                 (1.0 - q[i * nb_particles + j])) /
                q[i * nb_particles + j]);
      }

      else {
      }
    }
  }

  return -sum;
}

double SPH::calc_viscous_force(double *v) {

  double phisq;

  double sum = 0.0;                             // Initializing the sumation
  double pre = (40.0 / (M_PI * h * h * h * h)); // precalculated value

  for (int j = 0; j < nb_particles; j++) {

    if (i == j) {
    }

    else {

      if (q[i * nb_particles + j] < 1) {

        sum += (mass_assumed / particle_density[j]) * (v[i] - v[j]) *
               (pre * (1.0 - q[i * nb_particles + j]));
      }
    }
  }

  return -viscosity * sum;
}

// Function for the gravity force
double SPH::calc_gravity_force() {
  return -particle_density[i] * acceleration_gravity;
}

// Function for initialisation of the time integration scheme - velocity
double SPH::scheme_init(double *velocity, double &force_pressure,
                        double &force_viscous, double &force_gravity) {

  double acceleration;

  acceleration =
      (force_pressure + force_viscous + force_gravity) / particle_density[i];

  return velocity[i] + acceleration * dt * 0.5;
}

// Function for time integration - velocity
double SPH::velocity_integration(double *velocity, double &force_pressure,
                                 double &force_viscous, double &force_gravity) {

  double acceleration;
  acceleration =
      (force_pressure + force_viscous + force_gravity) / particle_density[i];

  return velocity[i] + acceleration * dt;
}

// Function to find the mass of the particles before the simulation starts
void SPH::calc_mass() {

  calc_particle_distance();
  calc_density();
  double sumden = 0.0;
  for (int i = 0; i < nb_particles; i++) {

    sumden += particle_density[i];
  }

  mass_assumed = nb_particles * density_resting / sumden;
}

// Function to perform the spatial iterations
void SPH::spatial() {

  for (i = 0; i < nb_particles; i++) {

    calc_pressure();

    // Gathering the forces calculated by the processors
    force_pressure_x = calc_pressure_force(position_x);

    force_viscous_x = calc_viscous_force(velocity_x);

    force_pressure_y = calc_pressure_force(position_y);

    force_viscous_y = calc_viscous_force(velocity_y);

    force_gravity_y = calc_gravity_force();

    // First step to initialise the scheme
    if (t == 0) {

      velocity_x[i] = scheme_init(velocity_x, force_pressure_x, force_viscous_x,
                                  force_gravity_x);
      position_x[i] = position_x[i] + velocity_x[i] * dt; // inlined
      velocity_y[i] = scheme_init(velocity_y, force_pressure_y, force_viscous_y,
                                  force_gravity_y);
      position_y[i] = position_y[i] + velocity_y[i] * dt; // inlined

    }

    // Leap frog scheme
    else {

      velocity_x[i] = velocity_integration(velocity_x, force_pressure_x,
                                           force_viscous_x, force_gravity_x);
      position_x[i] = position_x[i] + velocity_x[i] * dt; // inlined
      velocity_y[i] = velocity_integration(velocity_y, force_pressure_y,
                                           force_viscous_y, force_gravity_y);
      position_y[i] = position_y[i] + velocity_y[i] * dt; // inlined
    }

    // Boundary Conditions
    if (position_x[i] < h) {

      position_x[i] = h;
      velocity_x[i] = -coeff_restitution * velocity_x[i];
    }

    if (position_x[i] > 1.0 - h) {

      position_x[i] = 1.0 - h;
      velocity_x[i] = -coeff_restitution * velocity_x[i];
    }

    if (position_y[i] < h) {

      position_y[i] = h;
      velocity_y[i] = -coeff_restitution * velocity_y[i];
    }

    if (position_y[i] > 1.0 - h) {

      position_y[i] = 1.0 - h;
      velocity_y[i] = -coeff_restitution * velocity_y[i];
    }
  }
}

// Function to return the position x
double SPH::return_position_x(int l) { return position_x[l]; }

// Function to return the position y
double SPH::return_position_y(int l) { return position_y[l]; }

// Function to calculate the kinetic energy
double SPH::return_kinetic_energy() {

  double sum = 0;
  for (int i = 0; i < nb_particles; i++) {

    particle_speed[i] =
        velocity_x[i] * velocity_x[i] + velocity_y[i] * velocity_y[i];

    sum += particle_speed[i];
  }

  return 0.5 * mass_assumed * sum;
}

// Function to calculate the potential energy
double SPH::return_potential_energy() {

  double sum = 0;
  for (int i = 0; i < nb_particles; i++) {

    sum += position_y[i];
  }

  return mass_assumed * acceleration_gravity * sum;
}
