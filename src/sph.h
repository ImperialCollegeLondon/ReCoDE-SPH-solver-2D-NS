// Definition of the SPH-Calculation class. This should be a friend class of
// SPH main class and conduct only the calculations that are taking place.
#pragma once
#include "particles.h"

class sph {
  public:
    // This should be a purely calculations function. All methods are
    // going to be static and no object should be created.
    sph() = delete;
    ~sph() = delete;

    // Function to perform the particle iterations
    static void particle_iterations(particles&);

    static void update_position(particles&, int particle_index);
    // Function to initialise the time integration scheme - velocity
    static double scheme_init(particles&,int particle_index, double *velocity,
                      double &force_pressure, double &force_viscous,
                      double &force_gravity);

    // Function for time integration - velocity
    static double velocity_integration(particles&,int particle_index, double *velocity,
                                double &force_pressure, double &force_viscous,
                                double &force_gravity);

    // Function to treat the boundaries
    static void boundaries(particles&, int particle_index);

    // Function to calculate the matrix with rij
    static void calc_particle_distance(particles&);

    // Function to calculate the density
    static void calc_density(particles&);

    // Function to calculate the pressure
    static void calc_pressure(particles&);

    // Function to calculate the pressure force
    static double calc_pressure_force(particles&, int particle_index, double *position);

    // Function to calculate the viscous force
    static double calc_viscous_force(particles&, int particle_index, double *velocity);

    // Function to calculate the gravity force
    static double calc_gravity_force(particles&, int particle_index);

    // Function to find the mass of the particles before the simulation starts
    static void calc_mass(particles&);
};