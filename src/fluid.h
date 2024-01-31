#ifndef FLUID_H
#define FLUID_H
#include "fstream"
#include "particles.h"

class Fluid : public particles {
 private:
  // Constants of the problem
  double gasConstant;
  double densityResting;
  double viscosity;
  double accelerationGravity;
  double h;

  // Mass
  double mass = 1.0;

  // Density
  double *density;

  // Pressure
  double *pressure;

 public:
  Fluid() = default;  // Default constructor

  Fluid(const unsigned nNew);

  Fluid &operator=(const Fluid &fluid);

  // Setter functions

  // Assign value to gasConstant
  void setGasConstant(double gasConstant);

  // Assign value to densityResting
  void setDensityResting(double densityResting);

  // Assign value to the radius of influence
  void setRadInfl(double h);

  // Assign value to viscosity
  void setViscosity(double viscosity);

  // Assign value to accelerationGravity
  void setAccelerationGravity(double accelerationGravity);

  // Getter functions

  // Function to get the pressure felt by a single particle
  double getPressure(int index);

  // Function to get the density felt by a single particle
  double getDensity(int index);

  // Function to get the viscosity of the fluid
  double getViscosity();

  // Function to get the mass of the fluid
  double getMass();

  // Function to get the radius of influence
  double getRadInfl();

  // Function to get the gravitational acceleration
  double getAccelerationGravity();

  // Function to calculate and get the kinetic energy
  double getKineticEnergy();

  // Function to calculate and get the potential energy
  double getPotentialEnergy();

  // Calculation functions

  // Function to calculate the mass of the particles before the simulation
  // starts
  void calculateMass();

  // Function to calculate the density
  void calculateDensity();

  // Function to calculate the pressure
  void calculatePressure();
};
#endif