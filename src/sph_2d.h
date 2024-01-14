// Definition of the SPH-Calculation class. This should be a friend class of
// SPH main class and conduct only the calculations that are taking place.
#pragma once
#include "particles.h"
#include "fstream"

class sph_2d : public particles {
  public:

    using particles::particles;
    using particles::set_timestep;
    using particles::set_rad_infl;
    using particles::get_position_x;
    using particles::get_position_y;
    using particles::get_kinetic_energy;
    using particles::get_potential_energy;

    // Function to find the mass of the particles before the simulation starts
     void calc_mass();

        // Function to calculate the density
     void calc_density();

    // Function to calculate the pressure
     void calc_pressure();

    // Function to perform the time integration
    void time_integration(int nb_particles, int total_iter, double h,
                      double dt, std::ofstream &vOut, std::ofstream &vOut2);

    // Function to perform the particle iterations
     void particle_iterations();

    // Function to calculate the distance between the particles
     void calc_particle_distance();

     void update_position(int particle_index);

    // Function to initialise the time integration scheme - velocity
     double scheme_init(int particle_index, double *velocity,
                      double &force_pressure, double &force_viscous,
                      double &force_gravity);

    // Function for time integration - velocity
     double velocity_integration(int particle_index, double *velocity,
                                double &force_pressure, double &force_viscous,
                                double &force_gravity);

    // Function to calculate the pressure force
     double calc_pressure_force(int particle_index, double *position);

    // Function to calculate the viscous force
     double calc_viscous_force(int particle_index, double *velocity);

    // Function to calculate the gravity force
     double calc_gravity_force(int particle_index);

    // Function to treat the boundaries
     void boundaries(int particle_index);

};