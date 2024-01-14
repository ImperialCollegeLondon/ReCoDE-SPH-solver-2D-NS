#include "particles.h"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// User defined constructor
particles::particles(const unsigned n_new) : nb_particles(n_new) {

  position_x = new double[nb_particles];
  position_y = new double[nb_particles];
  velocity_x = new double[nb_particles];
  velocity_y = new double[nb_particles];

  distance = new double[nb_particles * nb_particles];
  distance_q = new double[nb_particles * nb_particles];

  particle_density = new double[nb_particles];

  particle_pressure = new double[nb_particles];

  particle_speed_sq = new double[nb_particles];
}

// Destructor
particles::~particles() {
  delete[] particle_density;
  delete[] particle_pressure;
  delete[] distance;
  delete[] distance_q;
  delete[] particle_speed_sq;
  delete[] position_x;
  delete[] position_y;
  delete[] velocity_x;
  delete[] velocity_y;
}

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

double particles::get_position_x(int l) const { return position_x[l]; }

double particles::get_position_y(int l) const { return position_y[l]; }

double particles::get_kinetic_energy() {

  double sum = 0;
  for (int i = 0; i < nb_particles; i++) {

    particle_speed_sq[i] =
        velocity_x[i] * velocity_x[i] + velocity_y[i] * velocity_y[i];

    sum += particle_speed_sq[i];
  }

  return 0.5 * mass_assumed * sum;
}

double particles::get_potential_energy() {

  double sum = 0;
  for (int i = 0; i < nb_particles; i++) {

    sum += position_y[i] - h;
  }

  return mass_assumed * acceleration_gravity * sum;
}
