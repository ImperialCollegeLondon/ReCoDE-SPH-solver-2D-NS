// Definition of the SPH-Calculation class. This should be a friend class of
// SPH main class and conduct only the calculations that are taking place.
#pragma once
#include "sph.h"

class SPH_Calc {
  public:
    // This should be a purely calculations function. All methods are
    // going to be static and no object should be created.
    SPH_Calc() = delete;
    ~SPH_Calc() = delete;

    // Function to calculate the matrix with rij
    static void calc_particle_distance(SPH&);

    // Function to calculate the density
    static void calc_density(SPH&);

    // Function to calculate the pressure
    static void calc_pressure(SPH&);

    // Function to calculate the pressure force
    static double calc_pressure_force(SPH&, int particle_index, double *position);

    // Function to calculate the viscous force
    static double calc_viscous_force(SPH&, int particle_index, double *velocity);

    // Function to calculate the gravity force
    static double calc_gravity_force(SPH&, int particle_index);

    // Function to find the mass of the particles before the simulation starts
    static void calc_mass(SPH&);
};