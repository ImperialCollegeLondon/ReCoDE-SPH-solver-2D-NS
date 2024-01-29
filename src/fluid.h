#pragma once
#include "fstream"
#include "particles.h"

class fluid : public particles {
 private:
  // Constants of the problem
  double gas_constant;
  double density_resting;
  double viscosity;
  double acceleration_gravity;

  double h;  // Radius of influence

  // Mass
  double mass = 1.0;

  // Density
  double *density;

  // Pressure
  double *pressure;

 public:
  fluid(const unsigned n_new);

  fluid &operator=(const fluid &fluid);

  // Assign value to gas_constant
  void set_gas_constant(double gas_constant);

  // Assign value to density_resting
  void set_density_resting(double density_resting);

  void set_rad_infl(double h);

  // Assign value to viscosity
  void set_viscosity(double viscosity);

  // Function to calculate the mass of the particles before the simulation
  // starts
  void calc_mass();

  // Function to calculate the density
  void calc_density();

  // Function to calculate the pressure
  void calc_pressure();

  // Assign value to acceleration_gravity
  void set_acceleration_gravity(double acceleration_gravity);

  // Function to get the pressure felt by a single particle
  double get_pressure(int index);

  // Function to get the density felt by a single particle
  double get_density(int index);

  double get_viscosity();

  double get_mass();

  double get_rad_infl();

  double get_acceleration_gravity();

  // Function to calculate the kinetic energy
  double get_kinetic_energy();

  // Function to calculate the potential energy
  double get_potential_energy();
};
