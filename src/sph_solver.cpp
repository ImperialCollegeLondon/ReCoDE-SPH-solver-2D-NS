#include "sph_solver.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "fluid.h"
#include "main_prog_funcs.h"

// Setter functions
void SphSolver::setTimestep(double dt) { this->dt = dt; }

void SphSolver::setTotalIterations(double totalIterations) {
  this->totalIterations = totalIterations;
}

void SphSolver::setOutputFrequency(double f) { this->outputFrequency = f; }

void SphSolver::setCoeffRestitution(double coeffRestitution) {
  this->coeffRestitution = coeffRestitution;
}

void SphSolver::setLeftWall(double leftWall) { this->leftWall = leftWall; }

void SphSolver::setRightWall(double rightWall) { this->rightWall = rightWall; }

void SphSolver::setBottomWall(double bottomWall) {
  this->bottomWall = bottomWall;
}

void SphSolver::setTopWall(double topWall) { this->topWall = topWall; }

void SphSolver::timeIntegration(Fluid &data, std::ofstream &finalPositionsFile,
                                std::ofstream &energiesFile) {
  std ::cout << "Time integration started -- OK"
             << "\n";

  numberOfParticles = data.getNumberOfParticles();

  for (int time = 0; time < totalIterations; time++) {
    t = time;
    // In each iteration the distances between the particles are recalculated,
    // as well as their density and pressure
    data.calculateParticleDistance();
    data.calculateDensity();
    data.calculatePressure();
    particleIterations(data);

    if (time % outputFrequency == 0) {
      storeToFile(data, "energy", energiesFile, dt, t);
    }
  }
  // Store particles' positions after integration is completed
  storeToFile(data, "position", finalPositionsFile, dt, totalIterations);

  std ::cout << "Time integration finished -- OK"
             << "\n";
}

void SphSolver::particleIterations(Fluid &data) {
  int i;

  // Use std::function to store the member functions
  std::function<double(int)> ptrGetPositionX =
      std::bind(&Fluid::getPositionX, &data, std::placeholders::_1);
  std::function<double(int)> ptrGetPositionY =
      std::bind(&Fluid::getPositionY, &data, std::placeholders::_1);
  std::function<double(int)> ptrGetVelocityX =
      std::bind(&Fluid::getVelocityX, &data, std::placeholders::_1);
  std::function<double(int)> ptrGetVelocityY =
      std::bind(&Fluid::getVelocityY, &data, std::placeholders::_1);

  for (i = 0; i < numberOfParticles; i++) {
    // Gathering the forces calculated by the processors
    forcePressureX = calculatePressureForce(data, ptrGetPositionX, i);

    forcePressureY = calculatePressureForce(data, ptrGetPositionY, i);

    forceViscousX = calcViscousForce(data, ptrGetVelocityX, i);

    forceViscousY = calcViscousForce(data, ptrGetVelocityY, i);

    forceGravityY = calcGravityForce(data, i);

    // Update the position of the particle
    updatePosition(data, i);

    // Boundary Conditions
    boundaries(data, i);
  }
}

double SphSolver::calculatePressureForce(Fluid &data,
                                         std::function<double(int)> getPosition,
                                         int particleIndex) {
  double sum = 0.0;  // Initializing the summation
  double radiusOfInfluence = data.getRadInfl();
  double thirtyPih3 =
      (-30.0 / (M_PI * radiusOfInfluence * radiusOfInfluence *
                radiusOfInfluence));  // Precalculated value used to avoid
                                      // multiple divisions and multiplications

  for (int j = 0; j < numberOfParticles; j++) {
    if (particleIndex != j) {
      if (data.getDistanceQ(particleIndex * numberOfParticles + j) < 1.0) {
        sum +=
            (data.getMass() / data.getDensity(j)) *
            ((data.getPressure(particleIndex) + data.getPressure(j)) / 2.0) *
            (thirtyPih3 * (getPosition(particleIndex) - getPosition(j))) *
            (((1.0 - data.getDistanceQ(particleIndex * numberOfParticles + j)) *
              (1.0 -
               data.getDistanceQ(particleIndex * numberOfParticles + j))) /
             data.getDistanceQ(particleIndex * numberOfParticles + j));
      }
    }
  }
  return -sum;
}

double SphSolver::calcViscousForce(Fluid &data,
                                   std::function<double(int)> getVelocity,
                                   int particleIndex) {
  double radiusOfInfluence = data.getRadInfl();

  double sum = 0.0;  // Initializing the summation
  double fourtyPih4 =
      (40.0 /
       (M_PI * radiusOfInfluence * radiusOfInfluence * radiusOfInfluence *
        radiusOfInfluence));  // Precalculated value used to avoid
                              // multiple divisions and multiplications

  for (int j = 0; j < numberOfParticles; j++) {
    if (particleIndex == j) {
    }

    else {
      if (data.getDistanceQ(particleIndex * numberOfParticles + j) < 1.0) {
        sum +=
            (data.getMass() / data.getDensity(j)) *
            (getVelocity(particleIndex) - getVelocity(j)) *
            (fourtyPih4 *
             (1.0 - data.getDistanceQ(particleIndex * numberOfParticles + j)));
      }
    }
  }

  return -data.getViscosity() * sum;
}

double SphSolver::calcGravityForce(Fluid &data, int particleIndex) {
  return -data.getDensity(particleIndex) * data.getAccelerationGravity();
}

void SphSolver::updatePosition(Fluid &data, int particleIndex) {
  double newVelocity;
  double newPosition;
  double integrationCoeff = 1.0;

  // First step to initialise the scheme
  if (t == 0) {
    integrationCoeff = 0.5;
  }

  // x-direction
  newVelocity = data.getVelocityX(particleIndex) +
                integrationCoeff *
                    velocityIntegration(data, particleIndex, forcePressureX,
                                        forceViscousX, forceGravityX);
  data.setVelocityX(particleIndex, newVelocity);

  newPosition = data.getPositionX(particleIndex) + newVelocity * dt;

  data.setPositionX(particleIndex, newPosition);

  // y-direction
  newVelocity = data.getVelocityY(particleIndex) +
                integrationCoeff *
                    velocityIntegration(data, particleIndex, forcePressureY,
                                        forceViscousY, forceGravityY);

  data.setVelocityY(particleIndex, newVelocity);

  newPosition = data.getPositionY(particleIndex) + newVelocity * dt;
  data.setPositionY(particleIndex, newPosition);
}

double SphSolver::velocityIntegration(Fluid &data, int particleIndex,
                                      double &forcePressure,
                                      double &forceViscous,
                                      double &forceGravity) {
  double acceleration;
  acceleration = (forcePressure + forceViscous + forceGravity) /
                 data.getDensity(particleIndex);

  return acceleration * dt;
}

void SphSolver::boundaries(Fluid &data, int particleIndex) {
  // x-direction
  if (data.getPositionX(particleIndex) < leftWall + data.getRadInfl()) {
    data.setPositionX(particleIndex, leftWall + data.getRadInfl());
    data.setVelocityX(particleIndex,
                      -coeffRestitution * data.getVelocityX(particleIndex));
  }

  if (data.getPositionX(particleIndex) > rightWall - data.getRadInfl()) {
    data.setPositionX(particleIndex, rightWall - data.getRadInfl());
    data.setVelocityX(particleIndex,
                      -coeffRestitution * data.getVelocityX(particleIndex));
  }

  // y-direction
  if (data.getPositionY(particleIndex) < bottomWall + data.getRadInfl()) {
    data.setPositionY(particleIndex, bottomWall + data.getRadInfl());
    data.setVelocityY(particleIndex,
                      -coeffRestitution * data.getVelocityY(particleIndex));
  }

  if (data.getPositionY(particleIndex) > topWall - data.getRadInfl()) {
    data.setPositionY(particleIndex, topWall - data.getRadInfl());
    data.setVelocityY(particleIndex,
                      -coeffRestitution * data.getVelocityY(particleIndex));
  }
}
