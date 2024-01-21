#pragma once
#include "particles.h"
#include "fstream"

class fluid : public particles {

// Constants of the problem
  const double gas_constant ;
  const double density_resting;
  const double viscosity;
  const double acceleration_gravity;

  double h; // Radius of influence

  // Mass
  double mass_assumed = 1.0;

  // Density
  double *fluid_density;

  // Pressure
  double *fluid_pressure;


};
