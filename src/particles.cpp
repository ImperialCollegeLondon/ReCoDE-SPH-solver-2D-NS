#include "particles.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// Default constructor
particles::particles() {}

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

// Destructor
particles::~particles() {
  delete[] position_x;
  delete[] particle_speed_sq;
  delete[] position_y;
  delete[] velocity_x;
  delete[] velocity_y;
  delete[] distance;
  delete[] distance_q;
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

// Overloading of operator=
particles &particles::operator=(const particles &particles) {
  if (this != &particles) {
    delete[] distance;
    delete[] distance_q;
    delete[] particle_speed_sq;
    delete[] position_x;
    delete[] position_y;
    delete[] velocity_x;
    delete[] velocity_y;

    nb_particles = particles.nb_particles;

    position_x = new double[nb_particles];
    position_y = new double[nb_particles];
    velocity_x = new double[nb_particles];
    velocity_y = new double[nb_particles];

    distance = new double[nb_particles * nb_particles];
    distance_q = new double[nb_particles * nb_particles];

    particle_speed_sq = new double[nb_particles];

    std::memcpy(position_x, particles.position_x,
                nb_particles * sizeof(double));
    std::memcpy(position_y, particles.position_y,
                nb_particles * sizeof(double));
    std::memcpy(velocity_x, particles.velocity_x,
                nb_particles * sizeof(double));
    std::memcpy(velocity_y, particles.velocity_y,
                nb_particles * sizeof(double));

    std::memcpy(distance, particles.distance,
                nb_particles * nb_particles * sizeof(double));
    std::memcpy(distance_q, particles.distance_q,
                nb_particles * nb_particles * sizeof(double));
    std::memcpy(particle_speed_sq, particles.particle_speed_sq,
                nb_particles * sizeof(double));
  }
  return *this;
}

// Getter functions

int particles::get_number_of_particles() { return nb_particles; }

double particles::get_position_x(int k) { return position_x[k]; }

double particles::get_position_y(int k) { return position_y[k]; }

double particles::get_distance_q(int k) { return distance_q[k]; }

// Calculation functions

void particles::calc_particle_distance() {
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
