#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include <cmath>

#include "fluid.h"

class SphSolver {
 public:
  // Setter Functions

  // Assign value to dt
  inline void setTimestep(double dt) { this->dt = dt; }

  // Assign value to the total iterations
  inline void setTotalIterations(double totalIterations) {
    this->totalIterations = totalIterations;
  }

  // Assign value to the frequency
  inline void setOutputFrequency(double outputFrequency) {
    this->outputFrequency = outputFrequency;
  }

  // Assign value to coeffRestitution
  inline void setCoeffRestitution(double coeffRestitution) {
    this->coeffRestitution = coeffRestitution;
  }

  // Assign value to leftWall
  inline void setLeftWall(double leftWall) { this->leftWall = leftWall; }

  // Assign value to rightWall
  inline void setRightWall(double rightWall) { this->rightWall = rightWall; }

  // Assign value to bottomWall
  inline void setBottomWall(double bottomWall) {
    this->bottomWall = bottomWall;
  }

  // Assign value to topWall
  inline void setTopWall(double topWall) { this->topWall = topWall; }

  // Assign values to pre-calculated values
  inline void setPrecalculatedValues(double radiusOfInfluence) {
    thirtyPih3 = -30.0 / (M_PI * std::pow(radiusOfInfluence, 3.0));

    fourtyPih4 = 40.0 / (M_PI * std::pow(radiusOfInfluence, 4.0));
  }

  // Calculation functions

  // Function to perform the time integration
  void timeIntegration(Fluid &data, std::ofstream &finalPositionsFile,
                       std::ofstream &energiesFile);

  // Function to perform the particle iterations
  void particleIterations(Fluid &data);

  // Function to calculate the pressure force
  double calculatePressureForce(Fluid &data,
                                std::function<double(int)> getPosition,
                                int particleIndex);

  // Function to calculate the viscous force
  double calcViscousForce(Fluid &data, std::function<double(int)> getVelocity,
                          int particleIndex);

  // Function to calculate the gravity force
  double calcGravityForce(Fluid &data, int particleIndex);

  // Function to update the positions of the particles
  void updatePosition(Fluid &data, int particleIndex);

  // Function for time integration - velocity
  double velocityIntegration(Fluid &data, int particleIndex,
                             double forcePressure, double forceViscous,
                             double forceGravity);

  // Function to treat the boundaries
  void boundaries(Fluid &data, int particleIndex);

 private:
  int numberOfParticles;

  int t;  // Time at a specific iteration

  int totalIterations;

  int outputFrequency;

  double dt;  // Timestep

  // Boundaries
  double coeffRestitution;

  double leftWall;
  double rightWall;
  double bottomWall;
  double topWall;

  // Pre calculated values
  double thirtyPih3;
  double fourtyPih4;

  // Forces
  double forcePressure, forceViscous, forceGravity;
  double forcePressureX, forcePressureY;
  double forceViscousX, forceViscousY;
  double forceGravityY;
  double forceGravityX = 0.0;
};

#endif
