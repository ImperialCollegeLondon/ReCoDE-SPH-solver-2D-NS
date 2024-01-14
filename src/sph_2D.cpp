#include "sph_2D.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

sph_2d::sph_2d(unsigned int nb_particles) : particles(nb_particles){}

void sph_2D::set_timestep(double dt) { this->dt = dt; }

void sph_2D::set_rad_infl(double h) { this->h = h; }

void sph_2D::calc_mass() {

  calc_particle_distance();
  calc_density();
  double sumden = 0.0;
  for (int i = 0; i < nb_particles; i++) {

    sumden += particle_density[i];
  }

  mass_assumed = nb_particles * density_resting / sumden;
}

void sph_2D::calc_particle_distance() {

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

void sph_2D::calc_density() {
  
  double phi;
  double four_pi_h_2 =
      (4.0 / (M_PI * h * h)); // Precalculated value used to avoid multiple
                              // divisions and multiplications
  double h_inverse =
      1.0 / h; // Precalculated value used to avoid multiple divisions

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

void sph_2D::time_integration(int nb_particles, int total_iter, double h,
                      double dt, std::ofstream &vOut, std::ofstream &vOut2) {

  std ::cout << "Time integration started -- OK"
             << "\n";

  for (int t = 0; t < total_iter; t++) {

    // In each iteration the distances between the particles are recalculated,
    // as well as their densities
    calc_particle_distance();
    calc_density();
    calc_pressure();

    particle_iterations();

    // Write energies on the Energy-File
    vOut2 << t * dt << "  " << get_kinetic_energy() << "  "
          << get_potential_energy() << "  "
          << get_potential_energy() + get_kinetic_energy()
          << "\n";

    // Get the positions after integration is completed
    if (t == total_iter - 1) {

      for (int l = 0; l < nb_particles; l++) {

        vOut << get_position_x(l) << " " << get_position_y(l)
             << "\n";
      }
    }
  }
}

void sph_2D::calc_pressure() {

  for (int i = 0; i < nb_particles; i++) {

    particle_pressure[i] =
        gas_constant * (particle_density[i] - density_resting);
  }
}

void sph_2D::particle_iterations() {

  int i;
  for (i = 0; i < nb_particles; i++) {

    // Calculate the forces acting on the particle
    force_pressure_x = calc_pressure_force(i, position_x);

    force_viscous_x = calc_viscous_force(i, velocity_x);

    force_pressure_y = calc_pressure_force(i, position_y);

    force_viscous_y = calc_viscous_force(i, velocity_y);

    force_gravity_y = calc_gravity_force(i);

    // Update the position of the particle
    update_position(i);

    // Boundary Conditions
    boundaries(i);
  }
}

double sph_2D::calc_pressure_force(int particle_index, double *position) {

  double sum = 0.0; // Initializing the sumation
  double thirty_pi_h_3 =
      (-30.0 / (M_PI * h * h * h)); // Precalculated value used to avoid
                                    // multiple divisions and multiplications

  for (int j = 0; j < nb_particles; j++) {

    if (particle_index != j) {

      if (distance_q[particle_index * nb_particles + j] < 1) {

        sum +=
            (mass_assumed / particle_density[j]) *
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

double sph_2D::calc_viscous_force(int particle_index, double *v) {

  double phisq;

  double sum = 0.0; // Initializing the sumation
  double fourty_pi_h_4 =
      (40.0 / (M_PI * h * h * h * h)); // Precalculated value used to avoid
                                       // multiple divisions and multiplications

  for (int j = 0; j < nb_particles; j++) {

    if (particle_index == j) {
    }

    else {

      if (distance_q[particle_index * nb_particles + j] < 1) {

        sum += (mass_assumed / particle_density[j]) *
               (v[particle_index] - v[j]) *
               (fourty_pi_h_4 *
                (1.0 - distance_q[particle_index * nb_particles + j]));
      }
    }
  }

  return -viscosity * sum;
}

double sph_2D::calc_gravity_force(int particle_index) {
  return -particle_density[particle_index] * acceleration_gravity;
}

void sph_2D::update_position(int particle_index){

  // First step to initialise the scheme
    if (t == 0) {

      velocity_x[particle_index] = scheme_init(particle_index, velocity_x, force_pressure_x,
                                  force_viscous_x, force_gravity_x);
      position_x[particle_index] = position_x[particle_index] + velocity_x[particle_index] * dt; 
      velocity_y[particle_index] = scheme_init(particle_index, velocity_y, force_pressure_y,
                                  force_viscous_y, force_gravity_y);
      position_y[particle_index] = position_y[particle_index] + velocity_y[particle_index] * dt; 

    }

    // Leap frog scheme
    else {

      velocity_x[particle_index] = velocity_integration(particle_index, velocity_x, force_pressure_x,
                                           force_viscous_x, force_gravity_x);
      position_x[particle_index] = position_x[particle_index] + velocity_x[particle_index] * dt; 
      velocity_y[particle_index] = velocity_integration(particle_index, velocity_y, force_pressure_y,
                                           force_viscous_y, force_gravity_y);
      position_y[particle_index] = position_y[particle_index] + velocity_y[particle_index] * dt; 
    }

}

double sph_2D::scheme_init(int particle_index, double *velocity,
                        double &force_pressure, double &force_viscous,
                        double &force_gravity) {

  double acceleration;

  acceleration = (force_pressure + force_viscous + force_gravity) /
                 particle_density[particle_index];

  return velocity[particle_index] + acceleration * dt * 0.5;
}

double sph_2D::velocity_integration(int particle_index, double *velocity,
                                 double &force_pressure, double &force_viscous,
                                 double &force_gravity) {

  double acceleration;
  acceleration = (force_pressure + force_viscous + force_gravity) /
                 particle_density[particle_index];

  return velocity[particle_index] + acceleration * dt;
}

void sph_2D::boundaries(int particle_index) {

if (position_x[particle_index] < h) {

      position_x[particle_index] = h;
      velocity_x[particle_index] = -coeff_restitution * velocity_x[particle_index];
    }

    if (position_x[particle_index] > 1.0 - h) {

      position_x[particle_index] = 1.0 - h;
      velocity_x[particle_index] = -coeff_restitution * velocity_x[particle_index];
    }

    if (position_y[particle_index] < h) {

      position_y[particle_index] = h;
      velocity_y[particle_index] = -coeff_restitution * velocity_y[particle_index];
    }

    if (position_y[particle_index] > 1.0 - h) {

      position_y[particle_index] = 1.0 - h;
      velocity_y[particle_index] = -coeff_restitution * velocity_y[particle_index];
    }
}




