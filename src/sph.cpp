#include "sph.h"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// User defined constructor
SPH::SPH(const unsigned n_new) : nb_particles(n_new) {

  position_x = new double[nb_particles];
  position_y = new double[nb_particles];
  velocity_x = new double[nb_particles];
  velocity_y = new double[nb_particles];

  distance = new double[nb_particles * nb_particles];
  distance_q = new double[nb_particles * nb_particles];

  particle_density = new double[nb_particles];

  particle_pressure = new double[nb_particles];

  particle_speed = new double[nb_particles];
}

// Destructor
SPH::~SPH() {
  delete[] particle_density;
  delete[] particle_pressure;
  delete[] distance;
  delete[] distance_q;
  delete[] particle_speed;
  delete[] position_x;
  delete[] position_y;
  delete[] velocity_x;
  delete[] velocity_y;
}

// Overloading of ()
double &SPH::operator()(unsigned row, unsigned col) {
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

// Assign value to t
void SPH::set_time(double t) { this->t = t; }

// Assign value to dt
void SPH::set_timestep(double dt) { this->dt = dt; }

// Assign value to h
void SPH::set_rad_infl(double h) { this->h = h; }

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

      distance_q[i * nb_particles + j] =
          std::abs(distance[i * nb_particles + j] * hinv);

      if (distance_q[i * nb_particles + j] < 1) {

        phi = pre *
              (1.0 - distance_q[i * nb_particles + j] *
                         distance_q[i * nb_particles + j]) *
              (1.0 - distance_q[i * nb_particles + j] *
                         distance_q[i * nb_particles + j]) *
              (1.0 - distance_q[i * nb_particles + j] *
                         distance_q[i * nb_particles + j]);

      }

      else {

        phi = 0.0;
      }

      particle_density[i] += mass_assumed * phi;
    }
  }
}

// Function to calculate the pressure
void SPH::calc_pressure() {

  for (int i = 0; i < nb_particles; i++) {

    particle_pressure[i] =
        gas_constant * (particle_density[i] - density_resting);
  }
}

// Function to calculate the pressure force
double SPH::calc_pressure_force(int index_i, double *x_y) {

  double sum = 0.0;                          // Initializing the sumation
  double pre = (-30.0 / (M_PI * h * h * h)); // precalculated value

  for (int j = 0; j < nb_particles; j++) {

    if (index_i == j) {
    } else {

      if (distance_q[index_i * nb_particles + j] < 1) {

        sum += (mass_assumed / particle_density[j]) *
               ((particle_pressure[index_i] + particle_pressure[j]) / 2.0) *
               (pre * (x_y[index_i] - x_y[j])) *
               (((1.0 - distance_q[index_i * nb_particles + j]) *
                 (1.0 - distance_q[index_i * nb_particles + j])) /
                distance_q[index_i * nb_particles + j]);
      }

      else {
      }
    }
  }

  return -sum;
}

// Function to calculate the viscous force
double SPH::calc_viscous_force(int index_i, double *v) {

  double phisq;

  double sum = 0.0;                             // Initializing the sumation
  double pre = (40.0 / (M_PI * h * h * h * h)); // precalculated value

  for (int j = 0; j < nb_particles; j++) {

    if (index_i == j) {
    }

    else {

      if (distance_q[index_i * nb_particles + j] < 1) {

        sum += (mass_assumed / particle_density[j]) * (v[index_i] - v[j]) *
               (pre * (1.0 - distance_q[index_i * nb_particles + j]));
      }
    }
  }

  return -viscosity * sum;
}

// Function to calculate the gravity force
double SPH::calc_gravity_force(int index_i) {
  return -particle_density[index_i] * acceleration_gravity;
}

// Function to initialise the time integration scheme - velocity
double SPH::scheme_init(int index_i, double *velocity, double &force_pressure,
                        double &force_viscous, double &force_gravity) {

  double acceleration;

  acceleration = (force_pressure + force_viscous + force_gravity) /
                 particle_density[index_i];

  return velocity[index_i] + acceleration * dt * 0.5;
}

// Function for time integration - velocity
double SPH::velocity_integration(int index_i, double *velocity,
                                 double &force_pressure, double &force_viscous,
                                 double &force_gravity) {

  double acceleration;
  acceleration = (force_pressure + force_viscous + force_gravity) /
                 particle_density[index_i];

  return velocity[index_i] + acceleration * dt;
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

  int i;

  for (i = 0; i < nb_particles; i++) {

    calc_pressure();

    // Gathering the forces calculated by the processors
    force_pressure_x = calc_pressure_force(i, position_x);

    force_viscous_x = calc_viscous_force(i, velocity_x);

    force_pressure_y = calc_pressure_force(i, position_y);

    force_viscous_y = calc_viscous_force(i, velocity_y);

    force_gravity_y = calc_gravity_force(i);

    // First step to initialise the scheme
    if (t == 0) {

      velocity_x[i] = scheme_init(i, velocity_x, force_pressure_x,
                                  force_viscous_x, force_gravity_x);
      position_x[i] = position_x[i] + velocity_x[i] * dt; // inlined
      velocity_y[i] = scheme_init(i, velocity_y, force_pressure_y,
                                  force_viscous_y, force_gravity_y);
      position_y[i] = position_y[i] + velocity_y[i] * dt; // inlined

    }

    // Leap frog scheme
    else {

      velocity_x[i] = velocity_integration(i, velocity_x, force_pressure_x,
                                           force_viscous_x, force_gravity_x);
      position_x[i] = position_x[i] + velocity_x[i] * dt; // inlined
      velocity_y[i] = velocity_integration(i, velocity_y, force_pressure_y,
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
