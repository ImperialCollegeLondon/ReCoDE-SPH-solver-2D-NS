#include "sph_solver.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "fluid.h"
#include "main_prog_funcs.h"

// Setter functions
void SphSolver::setTimestep(double dt) { this->dt = dt; }

void SphSolver::setTotalIterations(double totalIter) {
  this->totalterations = totalIter;
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

  for (int time = 0; time < totalterations; time++) {
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
  storeToFile(data, "position", finalPositionsFile, dt, totalterations);

  std ::cout << "Time integration finished -- OK"
             << "\n";
}

void SphSolver::particleIterations(Fluid &data) {
  int i;

  for (i = 0; i < numberOfParticles; i++) {
    // Gathering the forces calculated by the processors
    forcePressureX = calculatePressureForce(data, i, 0);

    forcePressureY = calculatePressureForce(data, i, 1);

    forceViscousX = calcViscousForce(data, i, 2);

    forceViscousY = calcViscousForce(data, i, 3);

    forceGravityY = calcGravityForce(data, i);

    // Update the position of the particle
    updatePosition(data, i);

    // Boundary Conditions
    boundaries(data, i);
  }
}

double SphSolver::calculatePressureForce(Fluid &data, int particleIndex,
                                         int dir) {
  double sum = 0.0;  // Initializing the summation
  double h = data.getRadInfl();
  double thirtyPiH3 =
      (-30.0 / (M_PI * h * h * h));  // Precalculated value used to avoid
                                     // multiple divisions and multiplications

  for (int j = 0; j < numberOfParticles; j++) {
    if (particleIndex != j) {
      if (data.getDistanceQ(particleIndex * numberOfParticles + j) < 1.0) {
        sum +=
            (data.getMass() / data.getDensity(j)) *
            ((data.getPressure(particleIndex) + data.getPressure(j)) / 2.0) *
            (thirtyPiH3 * (data(dir, particleIndex) - data(dir, j))) *
            (((1.0 - data.getDistanceQ(particleIndex * numberOfParticles + j)) *
              (1.0 -
               data.getDistanceQ(particleIndex * numberOfParticles + j))) /
             data.getDistanceQ(particleIndex * numberOfParticles + j));
      }
    }
  }
  return -sum;
}

double SphSolver::calcViscousForce(Fluid &data, int particleIndex, int dir) {
  double h = data.getRadInfl();

  double sum = 0.0;  // Initializing the summation
  double fourtyPiH4 = (40.0 / (M_PI * h * h * h *
                               h));  // Precalculated value used to avoid
                                     // multiple divisions and multiplications

  for (int j = 0; j < numberOfParticles; j++) {
    if (particleIndex == j) {
    }

    else {
      if (data.getDistanceQ(particleIndex * numberOfParticles + j) < 1.0) {
        sum +=
            (data.getMass() / data.getDensity(j)) *
            (data(dir, particleIndex) - data(dir, j)) *
            (fourtyPiH4 *
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
  // First step to initialise the scheme
  if (t == 0) {
    // x-direction
    data(2, particleIndex) =
        data(2, particleIndex) +
        0.5 * velocityIntegration(data, particleIndex, forcePressureX,
                                  forceViscousX, forceGravityX);
    data(0, particleIndex) =
        data(0, particleIndex) + data(2, particleIndex) * dt;

    // y-direction
    data(3, particleIndex) =
        data(3, particleIndex) +
        0.5 * velocityIntegration(data, particleIndex, forcePressureY,
                                  forceViscousY, forceGravityY);
    data(1, particleIndex) =
        data(1, particleIndex) + data(3, particleIndex) * dt;

  }

  // Leap frog scheme
  else {
    // x-direction
    data(2, particleIndex) =
        data(2, particleIndex) +
        velocityIntegration(data, particleIndex, forcePressureX, forceViscousX,
                            forceGravityX);
    data(0, particleIndex) =
        data(0, particleIndex) + data(2, particleIndex) * dt;

    // y-direction
    data(3, particleIndex) =
        data(3, particleIndex) +
        velocityIntegration(data, particleIndex, forcePressureY, forceViscousY,
                            forceGravityY);
    data(1, particleIndex) =
        data(1, particleIndex) + data(3, particleIndex) * dt;
  }
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
  if (data(0, particleIndex) < leftWall + data.getRadInfl()) {
    data(0, particleIndex) = leftWall + data.getRadInfl();
    data(2, particleIndex) = -coeffRestitution * data(2, particleIndex);
  }

  if (data(0, particleIndex) > rightWall - data.getRadInfl()) {
    data(0, particleIndex) = rightWall - data.getRadInfl();
    data(2, particleIndex) = -coeffRestitution * data(2, particleIndex);
  }

  // y-direction
  if (data(1, particleIndex) < bottomWall + data.getRadInfl()) {
    data(1, particleIndex) = bottomWall + data.getRadInfl();
    data(3, particleIndex) = -coeffRestitution * data(3, particleIndex);
  }

  if (data(1, particleIndex) > topWall - data.getRadInfl()) {
    data(1, particleIndex) = topWall - data.getRadInfl();
    data(3, particleIndex) = -coeffRestitution * data(3, particleIndex);
  }
}
