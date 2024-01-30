#ifndef FLUID_H
#define FLUID_H
#include "fstream"
#include "particles.h"

class fluid : public particles {
 private:
  // Constants of the problem
  double gas_constant;
  double density_resting;
  double viscosity;
  double acceleration_gravity;
  double h;

  // Mass
  double mass = 1.0;

  // Density
  double *density;

  // Pressure
  double *pressure;

 public:
  fluid();

  fluid(const unsigned n_new);

  fluid &operator=(const fluid &fluid);

  // Setter functions

  // Assign value to gas_constant
  void set_gas_constant(double gas_constant);

  // Assign value to density_resting
  void set_density_resting(double density_resting);

  // Assign value to the radius of influence
  void set_rad_infl(double h);

  // Assign value to viscosity
  void set_viscosity(double viscosity);

  // Assign value to acceleration_gravity
  void set_acceleration_gravity(double acceleration_gravity);

  // Getter functions

  // Function to get the pressure felt by a single particle
  double get_pressure(int index);

  // Function to get the density felt by a single particle
  double get_density(int index);

  // Function to get the viscosity of the fluid
  double get_viscosity();

  // Function to get the mass of the fluid
  double get_mass();

  // Function to get the radius of influence
  double get_rad_infl();

  // Function to get the gravitational acceleration
  double get_acceleration_gravity();

  // Function to calculate and get the kinetic energy
  double get_kinetic_energy();

  // Function to calculate and get the potential energy
  double get_potential_energy();

  // Calculation functions

  // Function to calculate the mass of the particles before the simulation
  // starts
  void calc_mass();

  // Function to calculate the density
  void calc_density();

  // Function to calculate the pressure
  void calc_pressure();
};
#endif