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

  particle_speed_sq = new double[nb_particles];
}

// Destructor
SPH::~SPH() {
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

void SPH::set_time(double t) { this->t = t; }

void SPH::set_timestep(double dt) { this->dt = dt; }

void SPH::set_rad_infl(double h) { this->h = h; }

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

void SPH::calc_density() {

  double phi;
  double four_pi_h_2 = (4.0 / (M_PI * h * h));  // Precalculated value used to avoid multiple divisions and multiplications
  double h_inverse = 1.0 / h;     // Precalculated value used to avoid multiple divisions

  // find φ
  for (int i = 0; i < nb_particles; i++) {

    particle_density[i] = 0;

    for (int j = 0; j < nb_particles; j++) {

      distance_q[i * nb_particles + j] =
          std::abs(distance[i * nb_particles + j] * h_inverse);

      if (distance_q[i * nb_particles + j] < 1) {

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

      particle_density[i] += mass_assumed * phi;
    }
  }
}

void SPH::calc_pressure() {

  for (int i = 0; i < nb_particles; i++) {

    particle_pressure[i] =
        gas_constant * (particle_density[i] - density_resting);
  }
}

double SPH::calc_pressure_force(int particle_index, double *position) {

  double sum = 0.0;                          // Initializing the sumation
  double thirty_pi_h_3 = (-30.0 / (M_PI * h * h * h)); //Precalculated value used to avoid multiple divisions and multiplications

  for (int j = 0; j < nb_particles; j++) {

    if (particle_index != j) {

      if (distance_q[particle_index * nb_particles + j] < 1) {

        sum += (mass_assumed / particle_density[j]) *
               ((particle_pressure[particle_index] + particle_pressure[j]) / 2.0) *
               (thirty_pi_h_3 * (position[particle_index] - position[j])) *
               (((1.0 - distance_q[particle_index * nb_particles + j]) *
                 (1.0 - distance_q[particle_index * nb_particles + j])) /
                distance_q[particle_index * nb_particles + j]);
      }

    }
  }

  return -sum;
}

double SPH::calc_viscous_force(int particle_index, double *v) {

  double phisq;

  double sum = 0.0;                             // Initializing the sumation
  double fourty_pi_h_4 = (40.0 / (M_PI * h * h * h * h)); //Precalculated value used to avoid multiple divisions and multiplications

  for (int j = 0; j < nb_particles; j++) {

    if (particle_index == j) {
    }

    else {

      if (distance_q[particle_index * nb_particles + j] < 1) {

        sum += (mass_assumed / particle_density[j]) * (v[particle_index] - v[j]) *
               (fourty_pi_h_4 * (1.0 - distance_q[particle_index * nb_particles + j]));
      }
    }
  }

  return -viscosity * sum;
}

double SPH::calc_gravity_force(int particle_index) {
  return -particle_density[particle_index] * acceleration_gravity;
}

double SPH::scheme_init(int particle_index, double *velocity, double &force_pressure,
                        double &force_viscous, double &force_gravity) {

  double acceleration;

  acceleration = (force_pressure + force_viscous + force_gravity) /
                 particle_density[particle_index];

  return velocity[particle_index] + acceleration * dt * 0.5;
}

double SPH::velocity_integration(int particle_index, double *velocity,
                                 double &force_pressure, double &force_viscous,
                                 double &force_gravity) {

  double acceleration;
  acceleration = (force_pressure + force_viscous + force_gravity) /
                 particle_density[particle_index];

  return velocity[particle_index] + acceleration * dt;
}

void SPH::calc_mass() {

  calc_particle_distance();
  calc_density();
  double sumden = 0.0;
  for (int i = 0; i < nb_particles; i++) {

    sumden += particle_density[i];
  }

  mass_assumed = nb_particles * density_resting / sumden;
}

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

double SPH::return_position_x(int l) { return position_x[l]; }

double SPH::return_position_y(int l) { return position_y[l]; }

double SPH::return_kinetic_energy() {

  double sum = 0;
  for (int i = 0; i < nb_particles; i++) {

    particle_speed_sq[i] =
        velocity_x[i] * velocity_x[i] + velocity_y[i] * velocity_y[i];

    sum += particle_speed_sq[i];
  }

  return 0.5 * mass_assumed * sum;
}

double SPH::return_potential_energy() {

  double sum = 0;
  for (int i = 0; i < nb_particles; i++) {

    sum += position_y[i] - h;
  }

  return mass_assumed * acceleration_gravity * sum;
}
