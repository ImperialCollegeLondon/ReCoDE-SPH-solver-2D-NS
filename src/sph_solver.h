#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include "fluid.h"

class SphSolver {
 private:
  int numberOfParticles;

  int t;  // Time at a specific iteration

  int totalterations;

  int outputFrequency;

  double dt;  // Timestep

  // Boundaries
  double coeffRestitution;

  double leftWall;
  double rightWall;
  double bottomWall;
  double topWall;

  // Forces
  double forcePressure, forceViscous, forceGravity;
  double forcePressureX, forcePressureY;
  double forceViscousX, forceViscousY;
  double forceGravityY;
  double forceGravityX = 0.0;

 public:
  // Setter Functions

  // Assign value to dt
  void setTimestep(double dt);

  // Assign value to the total iterations
  void setTotalIterations(double totalterations);

  // Assign value to the frequency
  void setOutputFrequency(double outputFrequency);

  // Assign value to coeffRestitution
  void setCoeffRestitution(double coeffRestitution);

  // Assign value to leftWall
  void setLeftWall(double leftWall);

  // Assign value to rightWall
  void setRightWall(double rightWall);

  // Assign value to bottomWall
  void setBottomWall(double bottomWall);

  // Assign value to topWall
  void setTopWall(double topWall);

  // Calculation functions

  // Function to perform the time integration
  void timeIntegration(Fluid &data, std::ofstream &finalPositionsFile,
                       std::ofstream &energiesFile);

  // Function to perform the particle iterations
  void particleIterations(Fluid &data);

  // Function to calculate the pressure force
  double calculatePressureForce(Fluid &data, int particleIndex, int dir);

  // Function to calculate the viscous force
  double calcViscousForce(Fluid &data, int particleIndex, int dir);

  // Function to calculate the gravity force
  double calcGravityForce(Fluid &data, int particleIndex);

  // Function to update the positions of the particles
  void updatePosition(Fluid &data, int particleIndex);

  // Function for time integration - velocity
  double velocityIntegration(Fluid &data, int particleIndex,
                             double &forcePressure, double &forceViscous,
                             double &forceGravity);

  // Function to treat the boundaries
  void boundaries(Fluid &data, int particleIndex);
};

#endif
