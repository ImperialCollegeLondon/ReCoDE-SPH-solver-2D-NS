#include "sph_solver.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "fluid.h"
#include "main_prog_funcs.h"

// Setter functions
void sph_solver::set_timestep(double dt) { this->dt = dt; }

void sph_solver::set_total_iter(double total_iter) {
  this->total_iterations = total_iter;
}

void sph_solver::set_output_frequency(double f) { this->output_frequency = f; }

void sph_solver::set_coeff_restitution(double coeff_restitution) {
  this->coeff_restitution = coeff_restitution;
}

void sph_solver::set_left_wall(double left_wall) {
  this->left_wall = left_wall;
}

void sph_solver::set_right_wall(double right_wall) {
  this->right_wall = right_wall;
}

void sph_solver::set_bottom_wall(double bottom_wall) {
  this->bottom_wall = bottom_wall;
}

void sph_solver::set_top_wall(double top_wall) { this->top_wall = top_wall; }

void sph_solver::time_integration(fluid &data,
                                  std::ofstream &finalPositionsFile,
                                  std::ofstream &energiesFile) {
  std ::cout << "Time integration started -- OK"
             << "\n";

  number_of_particles = data.get_number_of_particles();

  for (int time = 0; time < total_iterations; time++) {
    t = time;
    // In each iteration the distances between the particles are recalculated,
    // as well as their densities
    data.calc_particle_distance();
    data.calc_density();
    data.calc_pressure();
    particle_iterations(data);

    if (time % output_frequency == 0) {
      storeToFile(data, "energy", energiesFile, dt, t);
    }
  }
  // Store particles' positions after integration is completed
  storeToFile(data, "position", finalPositionsFile, dt, total_iterations);

  std ::cout << "Time integration finished -- OK"
             << "\n";
}

void sph_solver::particle_iterations(fluid &data) {
  int i;

  for (i = 0; i < number_of_particles; i++) {
    // Gathering the forces calculated by the processors
    force_pressure_x = calc_pressure_force(data, i, 0);

    force_pressure_y = calc_pressure_force(data, i, 1);

    force_viscous_x = calc_viscous_force(data, i, 2);

    force_viscous_y = calc_viscous_force(data, i, 3);

    force_gravity_y = calc_gravity_force(data, i);

    // Update the position of the particle
    update_position(data, i);

    // Boundary Conditions
    boundaries(data, i);
  }
}

double sph_solver::calc_pressure_force(fluid &data, int particle_index,
                                       int dir) {
  double sum = 0.0;  // Initializing the summation
  double h = data.get_rad_infl();
  double thirty_pi_h_3 =
      (-30.0 / (M_PI * h * h * h));  // Precalculated value used to avoid
                                     // multiple divisions and multiplications

  for (int j = 0; j < number_of_particles; j++) {
    if (particle_index != j) {
      if (data.get_distance_q(particle_index * number_of_particles + j) < 1.0) {
        sum +=
            (data.get_mass() / data.get_density(j)) *
            ((data.get_pressure(particle_index) + data.get_pressure(j)) / 2.0) *
            (thirty_pi_h_3 * (data(dir, particle_index) - data(dir, j))) *
            (((1.0 -
               data.get_distance_q(particle_index * number_of_particles + j)) *
              (1.0 -
               data.get_distance_q(particle_index * number_of_particles + j))) /
             data.get_distance_q(particle_index * number_of_particles + j));
      }
    }
  }
  return -sum;
}

double sph_solver::calc_viscous_force(fluid &data, int particle_index,
                                      int dir) {
  double h = data.get_rad_infl();

  double sum = 0.0;  // Initializing the summation
  double fourty_pi_h_4 =
      (40.0 /
       (M_PI * h * h * h * h));  // Precalculated value used to avoid
                                 // multiple divisions and multiplications

  for (int j = 0; j < number_of_particles; j++) {
    if (particle_index == j) {
    }

    else {
      if (data.get_distance_q(particle_index * number_of_particles + j) < 1.0) {
        sum += (data.get_mass() / data.get_density(j)) *
               (data(dir, particle_index) - data(dir, j)) *
               (fourty_pi_h_4 *
                (1.0 - data.get_distance_q(
                           particle_index * number_of_particles + j)));
      }
    }
  }

  return -data.get_viscosity() * sum;
}

double sph_solver::calc_gravity_force(fluid &data, int particle_index) {
  return -data.get_density(particle_index) * data.get_acceleration_gravity();
}

void sph_solver::update_position(fluid &data, int particle_index) {
  // First step to initialise the scheme
  if (t == 0) {
    data(2, particle_index) =
        data(2, particle_index) + scheme_init(data, particle_index,
                                              force_pressure_x, force_viscous_x,
                                              force_gravity_x);
    data(0, particle_index) =
        data(0, particle_index) + data(2, particle_index) * dt;

    data(3, particle_index) =
        data(3, particle_index) + scheme_init(data, particle_index,
                                              force_pressure_y, force_viscous_y,
                                              force_gravity_y);
    data(1, particle_index) =
        data(1, particle_index) + data(3, particle_index) * dt;

  }

  // Leap frog scheme
  else {
    data(2, particle_index) =
        data(2, particle_index) +
        velocity_integration(data, particle_index, force_pressure_x,
                             force_viscous_x, force_gravity_x);
    data(0, particle_index) =
        data(0, particle_index) + data(2, particle_index) * dt;

    data(3, particle_index) =
        data(3, particle_index) +
        velocity_integration(data, particle_index, force_pressure_y,
                             force_viscous_y, force_gravity_y);
    data(1, particle_index) =
        data(1, particle_index) + data(3, particle_index) * dt;
  }
}

///////////MERGE THE FOLLOWING TWO AFTER TESTING///////////////////

double sph_solver::scheme_init(fluid &data, int particle_index,
                               double &force_pressure, double &force_viscous,
                               double &force_gravity) {
  double acceleration;

  acceleration = (force_pressure + force_viscous + force_gravity) /
                 data.get_density(particle_index);

  return acceleration * dt * 0.5;
}

double sph_solver::velocity_integration(fluid &data, int particle_index,
                                        double &force_pressure,
                                        double &force_viscous,
                                        double &force_gravity) {
  double acceleration;
  acceleration = (force_pressure + force_viscous + force_gravity) /
                 data.get_density(particle_index);

  return acceleration * dt;
}

void sph_solver::boundaries(fluid &data, int particle_index) {
  if (data(0, particle_index) < left_wall + data.get_rad_infl()) {
    data(0, particle_index) = left_wall + data.get_rad_infl();
    data(2, particle_index) = -coeff_restitution * data(2, particle_index);
  }

  if (data(0, particle_index) > right_wall - data.get_rad_infl()) {
    data(0, particle_index) = right_wall - data.get_rad_infl();
    data(2, particle_index) = -coeff_restitution * data(2, particle_index);
  }

  if (data(1, particle_index) < bottom_wall + data.get_rad_infl()) {
    data(1, particle_index) = bottom_wall + data.get_rad_infl();
    data(3, particle_index) = -coeff_restitution * data(3, particle_index);
  }

  if (data(1, particle_index) > top_wall - data.get_rad_infl()) {
    data(1, particle_index) = top_wall - data.get_rad_infl();
    data(3, particle_index) = -coeff_restitution * data(3, particle_index);
  }
}
