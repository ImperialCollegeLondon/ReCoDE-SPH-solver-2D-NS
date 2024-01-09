#include "sph_calc.h"
#include <cmath>

void SPH_Calc::particle_iterations(SPH& data) {

  int i;
  for (i = 0; i < data.nb_particles; i++) {

    calc_pressure(data);

    // Gathering the forces
    data.force_pressure_x = calc_pressure_force(data, i, data.position_x);

    data.force_viscous_x = calc_viscous_force(data, i, data.velocity_x);

    data.force_pressure_y = calc_pressure_force(data, i, data.position_y);

    data.force_viscous_y = calc_viscous_force(data, i, data.velocity_y);

    data.force_gravity_y = calc_gravity_force(data, i);

    // Update the position of the particle
    update_position(data,i);

    // Boundary Conditions
    boundaries(data, i);
  }
}

void SPH_Calc::update_position(SPH& data,int particle_index){

  // First step to initialise the scheme
    if (data.t == 0) {

      data.velocity_x[particle_index] = scheme_init(data, particle_index, data.velocity_x, data.force_pressure_x,
                                  data.force_viscous_x, data.force_gravity_x);
      data.position_x[particle_index] = data.position_x[particle_index] + data.velocity_x[particle_index] * data.dt; 
      data.velocity_y[particle_index] = scheme_init(data,particle_index, data.velocity_y, data.force_pressure_y,
                                  data.force_viscous_y, data.force_gravity_y);
      data.position_y[particle_index] = data.position_y[particle_index] + data.velocity_y[particle_index] * data.dt; 

    }

    // Leap frog scheme
    else {

      data.velocity_x[particle_index] = velocity_integration(data, particle_index, data.velocity_x, data.force_pressure_x,
                                           data.force_viscous_x, data.force_gravity_x);
      data.position_x[particle_index] = data.position_x[particle_index] + data.velocity_x[particle_index] * data.dt; 
      data.velocity_y[particle_index] = velocity_integration(data,particle_index, data.velocity_y, data.force_pressure_y,
                                           data.force_viscous_y, data.force_gravity_y);
      data.position_y[particle_index] = data.position_y[particle_index] + data.velocity_y[particle_index] * data.dt; 
    }

}

double SPH_Calc::scheme_init(SPH& data,int particle_index, double *velocity,
                        double &force_pressure, double &force_viscous,
                        double &force_gravity) {

  double acceleration;

  acceleration = (force_pressure + force_viscous + force_gravity) /
                 data.particle_density[particle_index];

  return velocity[particle_index] + acceleration * data.dt * 0.5;
}

double SPH_Calc::velocity_integration(SPH& data,int particle_index, double *velocity,
                                 double &force_pressure, double &force_viscous,
                                 double &force_gravity) {

  double acceleration;
  acceleration = (force_pressure + force_viscous + force_gravity) /
                 data.particle_density[particle_index];

  return velocity[particle_index] + acceleration * data.dt;
}

void SPH_Calc::boundaries(SPH& data,int particle_index) {

if (data.position_x[particle_index] < data.h) {

      data.position_x[particle_index] = data.h;
     data. velocity_x[particle_index] = -data.coeff_restitution * data.velocity_x[particle_index];
    }

    if (data.position_x[particle_index] > 1.0 - data.h) {

     data. position_x[particle_index] = 1.0 - data.h;
      data.velocity_x[particle_index] = -data.coeff_restitution * data.velocity_x[particle_index];
    }

    if (data.position_y[particle_index] < data.h) {

      data.position_y[particle_index] = data.h;
      data.velocity_y[particle_index] = -data.coeff_restitution * data.velocity_y[particle_index];
    }

    if (data.position_y[particle_index] > 1.0 - data.h) {

      data.position_y[particle_index] = 1.0 - data.h;
      data.velocity_y[particle_index] = -data.coeff_restitution * data.velocity_y[particle_index];
    }
}

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