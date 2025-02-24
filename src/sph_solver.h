#ifndef SPHSOLVER_H
#define SPHSOLVER_H

#include <cmath>

#include "fluid.h"

class SphSolver {
 public:
  // Setter Functions

  // Determine whether to use adaptive timestep or not
  inline void setAdaptiveTimestep(bool adaptiveTimestepBool) {
    this->adaptiveTimestepBool = adaptiveTimestepBool;
  }

  void setCflCoefficients(double coeffCfl1, double coeffCfl2) {
    this->coeffCfl1 = coeffCfl1;
    this->coeffCfl2 = coeffCfl2;
  }

  // Assign value to dt
  inline void setTimestep(double dt) {
    if (!adaptiveTimestepBool) {
      this->dt = dt;
    } else {
      this->dt = initialTimestep;
    }
  }

  // Assign value to the total integration time
  inline void setTotalTime(double totalTime) { this->totalTime = totalTime; }

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

  // Getter Functions
  const std::vector<std::vector<std::pair<int, double>>> &
  getNeighbourParticles() const;

  // Neighbour search functions

  void createGrid(Fluid &data);

  void assignNeighbourCells(int rows, int cols);

  void placeParticlesInCells(Fluid &data);

  void neighbourParticlesSearch(Fluid &data);

  // Calculation functions

  // Function to perform the time integration
  void timeIntegration(Fluid &data, std::ofstream &simulationPositionsFile,
                       std::ofstream &finalPositionsFile,
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

  // Function to update the timestep
  void adaptiveTimestep(Fluid &data);

 private:
  constexpr static double initialTimestep = 1e-4;
  constexpr static int maxNeighbourCells = 8;
  constexpr static double memoryReservationFactor = 1.1;
  constexpr static double forceGravityX = 0.0;

  int numberOfParticles;
  std::vector<std::vector<std::pair<int, double>>> neighbourParticles;

  int numberOfCells;
  std::vector<std::vector<int>> cells;
  std::vector<std::vector<int>> neighbourCells;

  bool adaptiveTimestepBool = false;

  size_t t = 0;

  int outputFrequency;

  double currentIntegrationTime = 0.0;

  double totalTime;

  double dt;  // Timestep

  double coeffCfl1;  // CFL coefficients
  double coeffCfl2;  // CFL coefficients

  bool termination_flag = false;  // Flag to terminate the integration function
                                  // when the particles stop moving

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

  // Adaptive timestep related variables
  double maxVelocity = 0.0;      // maximum velocity
  double maxAcceleration = 0.0;  // maxumum force per unit mass
};

#endif
