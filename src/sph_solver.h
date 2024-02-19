#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include "fluid.h"

class SphSolver {
 private:
  int numberOfParticles;
  std::vector<std::vector<std::pair<int, double>>> neighbourParticles;

  int numberOfCells;
  std::vector<std::vector<int>> cells;
  std::vector<std::vector<int>> neighbourCells;

  // Time related variables
  bool adaptiveTimestepBool = false;

  int t = 0;

  int outputFrequency;

  double currentIntegrationTime = 0.0;

  double totalTime;

  double dt;  // Timestep

  double coeffCfl1;  // CFL coefficients
  double coeffCfl2;  // CFL coefficients

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

  const int MAX_NEIGHBOUR_CELLS = 8;

  double memoryReservationFactor = 1.1;

  // Adaptive timestep related variables
  double maxVelocity = 0.0;      // maximum velocity
  double maxAcceleration = 0.0;  // maxumum force per unit mass

 public:
  // Setter Functions

  // Determine whether to use adaptive timestep or not
  void setAdaptiveTimestep(bool adaptiveTimestepBool);

  // Determine whether to use adaptive timestep or not
  void setCflCoefficients(double coeffCfl1, double coeffCfl2);

  // Assign value to dt
  void setTimestep(double dt);

  // Assign value to the total integration time
  void setTotalTime(double totalTime);

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

  // Assign values to pre-calculated values
  void setPrecalculatedValues(double radiusOfInfluence);

  // Getter Functions
  std::vector<std::vector<std::pair<int, double>>> getNeighbourParticles();

  // Neighbour search functions

  void createGrid(Fluid &data);

  void assignNeighbourCells(int rows, int cols);

  void placeParticlesInCells(Fluid &data);

  void neighbourParticlesSearch(Fluid &data);

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
                             double &forcePressure, double &forceViscous,
                             double &forceGravity);

  // Function to treat the boundaries
  void boundaries(Fluid &data, int particleIndex);

  // Function to update the timestep
  void adaptiveTimestep(Fluid &data);
};

#endif
