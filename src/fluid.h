#ifndef FLUID_H
#define FLUID_H

#include <cmath>
#include <execution>
#include <fstream>

#include "particles.h"

class Fluid : public particles {
 public:
  explicit Fluid(const unsigned nNew);

  // Setter functions

  // Assign value to gasConstant
  inline void setGasConstant(double gasConstant) {
    this->gasConstant = gasConstant;
  }

  // Assign value to densityResting
  inline void setDensityResting(double densityResting) {
    this->densityResting = densityResting;
  }

  // Assign value to the radius of influence
  inline void setRadInfl(double radiusOfInfluence) {
    this->radiusOfInfluence = radiusOfInfluence;
    hInverse = 1.0 / radiusOfInfluence;
    fourPih2 = (4.0 / (M_PI * radiusOfInfluence * radiusOfInfluence));
  }

  // Assign value to viscosity
  inline void setViscosity(double viscosity) { this->viscosity = viscosity; }

  // Assign value to accelerationGravity
  inline void setAccelerationGravity(double accelerationGravity) {
    this->accelerationGravity = accelerationGravity;
  }

  // Getter functions

  // Function to get the pressure felt by a single particle
  inline double getPressure(int index) { return pressure[index]; }

  // Function to get the density felt by a single particle
  inline double getDensity(int index) { return density[index]; }

  // Function to get the viscosity of the fluid
  inline double getViscosity() { return viscosity; }

  // Function to get the mass of the fluid
  inline double getMass() { return mass; }

  // Function to get the gravitational acceleration
  inline double getAccelerationGravity() { return accelerationGravity; }

  // Function to get the radius of influence
  inline double getRadInfl() { return radiusOfInfluence; }

  // Function to calculate and get the kinetic energy
  inline double getKineticEnergy() {
    double sum = 0;
    for (size_t i = 0; i < nbParticles; i++) {
      particleSpeedSq[i] =
          velocityX[i] * velocityX[i] + velocityY[i] * velocityY[i];

      sum += particleSpeedSq[i];
    }

    return 0.5 * mass * sum;
  }

  // Function to calculate and get the potential energy
  inline double getPotentialEnergy() {
    double sum = 0;

    // std::for_each( std::execution::par, positionY.begin(), positionY.end(),
    // [&](const auto &Y){
    //   sum += Y - radiusOfInfluence;
    // });
    for (size_t i = 0; i < nbParticles; i++) {
      sum += positionY[i] - radiusOfInfluence;
    }

    return mass * accelerationGravity * sum;
  }

  // Calculation functions

  // Function to calculate the mass of the particles before the simulation
  // starts
  void calculateMass();

  // Function to calculate the density
  void calculateDensity();

  // Function to calculate the pressure
  void calculatePressure();

 private:
  // Constants of the problem
  double gasConstant;
  double densityResting;
  double viscosity;
  double accelerationGravity;
  double radiusOfInfluence;
  // Helper member variables
  double hInverse;
  double fourPih2;

  // Mass
  double mass = 1.0;

  // Density
  std::vector<double> density;

  // Pressure
  std::vector<double> pressure;
};
#endif