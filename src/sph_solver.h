#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include "fluid.h"

class sph_solver {
 private:
  int t;  // time

  int number_of_particles;

  int total_iterations;

  int output_frequency;

  double dt;  // time-step

  // Boundaries
  double coeff_restitution;
  double left_wall;
  double right_wall;
  double bottom_wall;
  double top_wall;

  // Forces
  double force_pressure, force_viscous, force_gravity;
  double force_pressure_x, force_pressure_y;
  double force_viscous_x, force_viscous_y;
  double force_gravity_y;
  double force_gravity_x = 0.0;

 public:
  /**********MEMBER-FUNCTIONS*********/

  // Setter Functions.

  // Assign value to dt
  void set_timestep(double dt);

  // Assign value to the total iterations
  void set_total_iter(double total_iter);

  // Assign value to the frequency
  void set_output_frequency(double output_frequency);

  // Assign value to coeff_restitution
  void set_coeff_restitution(double coeff_restitution);

  // Assign value to left_wall
  void set_left_wall(double left_wall);

  // Assign value to right_wall
  void set_right_wall(double right_wall);

  // Assign value to bottom_wall
  void set_bottom_wall(double bottom_wall);

  // Assign value to top_wall
  void set_top_wall(double top_wall);

  // Function to perform the time integration
  void time_integration(fluid &data, std::ofstream &finalPositionsFile,
                        std::ofstream &energiesFile);

  // Function to perform the particle iterations
  void particle_iterations(fluid &data);

  // Function to calculate the pressure force
  double calc_pressure_force(fluid &data, int particle_index, int dir);

  // Function to calculate the viscous force
  double calc_viscous_force(fluid &data, int particle_index, int dir);

  // Function to calculate the gravity force
  double calc_gravity_force(fluid &data, int particle_index);

  // Function to update the positions of the particles
  void update_position(fluid &data, int particle_index);

  // Function for time integration - velocity
  double velocity_integration(fluid &data, int particle_index,
                              double &force_pressure, double &force_viscous,
                              double &force_gravity);

  // Function to treat the boundaries
  void boundaries(fluid &data, int particle_index);
};

#endif