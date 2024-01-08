#include "sph.h"
#include "sph_calc.h"
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

void SPH::set_timestep(double dt) { this->dt = dt; }

void SPH::set_rad_infl(double h) { this->h = h; }

double SPH::scheme_init(int particle_index, double *velocity,
                        double &force_pressure, double &force_viscous,
                        double &force_gravity) {

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

void SPH::particle_iterations() {

  int i;
  for (i = 0; i < nb_particles; i++) {

    SPH_Calc::calc_pressure(*this);

    // Gathering the forces calculated by the processors
    force_pressure_x = SPH_Calc::calc_pressure_force(*this, i, position_x);

    force_viscous_x = SPH_Calc::calc_viscous_force(*this, i, velocity_x);

    force_pressure_y = SPH_Calc::calc_pressure_force(*this, i, position_y);

    force_viscous_y = SPH_Calc::calc_viscous_force(*this, i, velocity_y);

    force_gravity_y = SPH_Calc::calc_gravity_force(*this, i);

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

double SPH::get_position_x(int l) const { return position_x[l]; }

double SPH::get_position_y(int l) const { return position_y[l]; }

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
