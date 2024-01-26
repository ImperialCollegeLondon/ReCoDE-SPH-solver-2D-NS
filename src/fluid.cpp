#include "fluid.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
// User defined constructor
fluid::fluid(const unsigned n_new) : particles(n_new) {
  pressure = new double[nb_particles];
  density = new double[nb_particles];
}

fluid &fluid::operator=(const fluid &fluid) {
  if (this != &fluid) {
    delete[] distance;
    delete[] distance_q;
    delete[] particle_speed_sq;
    delete[] position_x;
    delete[] position_y;
    delete[] velocity_x;
    delete[] velocity_y;
    delete[] pressure;
    delete[] density;

    nb_particles = fluid.nb_particles;

    position_x = new double[nb_particles];
    position_y = new double[nb_particles];
    velocity_x = new double[nb_particles];
    velocity_y = new double[nb_particles];

    distance = new double[nb_particles * nb_particles];
    distance_q = new double[nb_particles * nb_particles];

    particle_speed_sq = new double[nb_particles];

    pressure = new double[nb_particles];
    density = new double[nb_particles];

    std::memcpy(position_x, fluid.position_x, nb_particles * sizeof(double));
    std::memcpy(position_y, fluid.position_y, nb_particles * sizeof(double));
    std::memcpy(velocity_x, fluid.velocity_x, nb_particles * sizeof(double));
    std::memcpy(velocity_y, fluid.velocity_y, nb_particles * sizeof(double));

    std::memcpy(distance, fluid.distance,
                nb_particles * nb_particles * sizeof(double));
    std::memcpy(distance_q, fluid.distance_q,
                nb_particles * nb_particles * sizeof(double));
    std::memcpy(particle_speed_sq, fluid.particle_speed_sq,
                nb_particles * sizeof(double));

    std::memcpy(pressure, fluid.pressure, nb_particles * sizeof(double));
    std::memcpy(density, fluid.density, nb_particles * sizeof(double));
  }
  return *this;
}

void fluid::set_gas_constant(double gas_constant) {
  this->gas_constant = gas_constant;
}

void fluid::set_density_resting(double density_resting) {
  this->density_resting = density_resting;
}

void fluid::set_rad_infl(double h) { this->h = h; }

void fluid::set_viscosity(double viscosity) { this->viscosity = viscosity; }

void fluid::set_acceleration_gravity(double acceleration_gravity) {
  this->acceleration_gravity = acceleration_gravity;
}

void fluid::calc_mass() {
  calc_particle_distance();
  calc_density();
  double sumden = 0.0;
  for (int i = 0; i < nb_particles; i++) {
    sumden += density[i];
  }

  mass_assumed = nb_particles * density_resting / sumden;
}

void fluid::calc_density() {
  double phi;
  double four_pi_h_2 =
      (4.0 / (M_PI * h * h));  // Precalculated value used to avoid multiple
                               // divisions and multiplications
  double h_inverse =
      1.0 / h;  // Precalculated value used to avoid multiple divisions

  // find φ
  for (int i = 0; i < nb_particles; i++) {
    density[i] = 0;

    for (int j = 0; j < nb_particles; j++) {
      distance_q[i * nb_particles + j] =
          std::abs(distance[i * nb_particles + j] * h_inverse);

      if (distance_q[i * nb_particles + j] < 1.0) {
        phi = four_pi_h_2 *
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

      density[i] += mass_assumed * phi;
    }
  }
}

void fluid::calc_pressure() {
  for (int i = 0; i < nb_particles; i++) {
    pressure[i] = gas_constant * (density[i] - density_resting);
  }
}

double fluid::get_pressure(int index) { return pressure[index]; }

double fluid::get_density(int index) { return density[index]; }

double fluid::get_viscosity() { return viscosity; }

double fluid::get_mass() { return mass_assumed; }

double fluid::get_acceleration_gravity() { return acceleration_gravity; }

double fluid::get_rad_infl() { return h; }

double fluid::get_kinetic_energy() {
  double sum = 0;
  for (int i = 0; i < nb_particles; i++) {
    particle_speed_sq[i] =
        velocity_x[i] * velocity_x[i] + velocity_y[i] * velocity_y[i];

    sum += particle_speed_sq[i];
  }

  return 0.5 * mass_assumed * sum;
}

double fluid::get_potential_energy() {
  double sum = 0;

  for (int i = 0; i < nb_particles; i++) {
    sum += position_y[i] - h;
  }

  return mass_assumed * acceleration_gravity * sum;
}