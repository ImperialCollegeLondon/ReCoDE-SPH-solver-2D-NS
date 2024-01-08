// Definition of the SPH-Calculation class. This should be a friend class of
// SPH main class and conduct only the calculations that are taking place.
#pragma once
#include "sph.h"

class SPH_Calc {
  public:
    SPH_Calc() = default;
    ~SPH_Calc() = default;

    // Function to calculate the matrix with rij
    void calc_particle_distance(SPH&);

    // Function to calculate the density
    void calc_density(SPH&);

    // Function to calculate the pressure
    void calc_pressure(SPH&);

    // Function to calculate the pressure force
    double calc_pressure_force(SPH&, int particle_index, double *position);

    // Function to calculate the viscous force
    double calc_viscous_force(SPH&, int particle_index, double *velocity);

    // Function to calculate the gravity force
    double calc_gravity_force(SPH&, int particle_index);

    // Function to find the mass of the particles before the simulation starts
    void calc_mass(SPH&);
};