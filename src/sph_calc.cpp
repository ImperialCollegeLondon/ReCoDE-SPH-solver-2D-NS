#include "sph_calc.h"
#include <cmath>

void SPH_Calc::calc_particle_distance(SPH& data) {

  double dx;
  double dy;

  for (int i = 0; i < data.nb_particles; i++) {

    for (int j = 0; j < data.nb_particles; j++) {

      dx = data.position_x[i] - data.position_x[j];
      dy = data.position_y[i] - data.position_y[j];

      data.distance[i * data.nb_particles + j] = sqrt(dx * dx + dy * dy);
    }
  }
}

void SPH_Calc::calc_density(SPH& data) {
  
  double phi;
  double four_pi_h_2 =
      (4.0 / (M_PI * data.h * data.h)); // Precalculated value used to avoid multiple
                              // divisions and multiplications
  double h_inverse =
      1.0 / data.h; // Precalculated value used to avoid multiple divisions

  // find φ
  for (int i = 0; i < data.nb_particles; i++) {

    data.particle_density[i] = 0;

    for (int j = 0; j < data.nb_particles; j++) {

      data.distance_q[i * data.nb_particles + j] =
          std::abs(data.distance[i * data.nb_particles + j] * h_inverse);

      if (data.distance_q[i * data.nb_particles + j] < 1) {

        phi = four_pi_h_2 *
              (1.0 - data.distance_q[i * data.nb_particles + j] *
                         data.distance_q[i * data.nb_particles + j]) *
              (1.0 - data.distance_q[i * data.nb_particles + j] *
                         data.distance_q[i * data.nb_particles + j]) *
              (1.0 - data.distance_q[i * data.nb_particles + j] *
                         data.distance_q[i * data.nb_particles + j]);

      }

      else {

        phi = 0.0;
      }

      data.particle_density[i] += data.mass_assumed * phi;
    }
  }
}

void SPH_Calc::calc_pressure(SPH& data) {

  for (int i = 0; i < data.nb_particles; i++) {

    data.particle_pressure[i] =
        data.gas_constant * (data.particle_density[i] - data.density_resting);
  }
}

double SPH_Calc::calc_pressure_force(SPH& data, int particle_index, double *position) {

  double sum = 0.0; // Initializing the sumation
  double thirty_pi_h_3 =
      (-30.0 / (M_PI * data.h * data.h * data.h)); // Precalculated value used to avoid
                                    // multiple divisions and multiplications

  for (int j = 0; j < data.nb_particles; j++) {

    if (particle_index != j) {

      if (data.distance_q[particle_index * data.nb_particles + j] < 1) {

        sum +=
            (data.mass_assumed / data.particle_density[j]) *
            ((data.particle_pressure[particle_index] + data.particle_pressure[j]) / 2.0) *
            (thirty_pi_h_3 * (position[particle_index] - position[j])) *
            (((1.0 - data.distance_q[particle_index * data.nb_particles + j]) *
              (1.0 - data.distance_q[particle_index * data.nb_particles + j])) /
             data.distance_q[particle_index * data.nb_particles + j]);
      }
    }
  }

  return -sum;
}

double SPH_Calc::calc_viscous_force(SPH& data, int particle_index, double *v) {

  double phisq;

  double sum = 0.0; // Initializing the sumation
  double fourty_pi_h_4 =
      (40.0 / (M_PI * data.h * data.h * data.h * data.h)); // Precalculated value used to avoid
                                       // multiple divisions and multiplications

  for (int j = 0; j < data.nb_particles; j++) {

    if (particle_index == j) {
    }

    else {

      if (data.distance_q[particle_index * data.nb_particles + j] < 1) {

        sum += (data.mass_assumed / data.particle_density[j]) *
               (v[particle_index] - v[j]) *
               (fourty_pi_h_4 *
                (1.0 - data.distance_q[particle_index * data.nb_particles + j]));
      }
    }
  }

  return -data.viscosity * sum;
}

double SPH_Calc::calc_gravity_force(SPH& data, int particle_index) {
  return -data.particle_density[particle_index] * data.acceleration_gravity;
}

void SPH_Calc::calc_mass(SPH& data) {

  calc_particle_distance(data);
  calc_density(data);
  double sumden = 0.0;
  for (int i = 0; i < data.nb_particles; i++) {

    sumden += data.particle_density[i];
  }

  data.mass_assumed = data.nb_particles * data.density_resting / sumden;
}