#include "fluid.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// User defined constructor
Fluid::Fluid(const unsigned nNew) : particles(nNew) {
  pressure = new double[nbParticles];
  density = new double[nbParticles];
}
// Copy constructor
Fluid &Fluid::operator=(const Fluid &fluid) {
  if (this != &fluid) {
    delete[] distance;
    delete[] distanceQ;
    delete[] particleSpeedSq;
    delete[] positionX;
    delete[] positionY;
    delete[] velocityX;
    delete[] velocityY;
    delete[] pressure;
    delete[] density;

    nbParticles = fluid.nbParticles;

    positionX = new double[nbParticles];
    positionY = new double[nbParticles];
    velocityX = new double[nbParticles];
    velocityY = new double[nbParticles];

    distance = new double[nbParticles * nbParticles];
    distanceQ = new double[nbParticles * nbParticles];

    particleSpeedSq = new double[nbParticles];

    pressure = new double[nbParticles];
    density = new double[nbParticles];

    std::memcpy(positionX, fluid.positionX, nbParticles * sizeof(double));
    std::memcpy(positionY, fluid.positionY, nbParticles * sizeof(double));
    std::memcpy(velocityX, fluid.velocityX, nbParticles * sizeof(double));
    std::memcpy(velocityY, fluid.velocityY, nbParticles * sizeof(double));

    std::memcpy(distance, fluid.distance,
                nbParticles * nbParticles * sizeof(double));
    std::memcpy(distanceQ, fluid.distanceQ,
                nbParticles * nbParticles * sizeof(double));
    std::memcpy(particleSpeedSq, fluid.particleSpeedSq,
                nbParticles * sizeof(double));

    std::memcpy(pressure, fluid.pressure, nbParticles * sizeof(double));
    std::memcpy(density, fluid.density, nbParticles * sizeof(double));
  }
  return *this;
}

// Setter functions

void Fluid::setGasConstant(double gasConstant) {
  this->gasConstant = gasConstant;
}

void Fluid::setDensityResting(double densityResting) {
  this->densityResting = densityResting;
}

void Fluid::setRadInfl(double radiusOfInfluence) {
  this->radiusOfInfluence = radiusOfInfluence;
}

void Fluid::setViscosity(double viscosity) { this->viscosity = viscosity; }

void Fluid::setAccelerationGravity(double accelerationGravity) {
  this->accelerationGravity = accelerationGravity;
}

// Getter functions

double Fluid::getPressure(int index) { return pressure[index]; }

double Fluid::getDensity(int index) { return density[index]; }

double Fluid::getViscosity() { return viscosity; }

double Fluid::getMass() { return mass; }

double Fluid::getAccelerationGravity() { return accelerationGravity; }

double Fluid::getRadInfl() { return radiusOfInfluence; }

double Fluid::getKineticEnergy() {
  double sum = 0;
  for (int i = 0; i < nbParticles; i++) {
    particleSpeedSq[i] =
        velocityX[i] * velocityX[i] + velocityY[i] * velocityY[i];

    sum += particleSpeedSq[i];
  }

  return 0.5 * mass * sum;
}

double Fluid::getPotentialEnergy() {
  double sum = 0;

  for (int i = 0; i < nbParticles; i++) {
    sum += positionY[i] - radiusOfInfluence;
  }

  return mass * accelerationGravity * sum;
}

// Calculation functions

void Fluid::calculateMass() {
  calculateParticleDistance();
  calculateDensity();
  double sumDensity = 0.0;
  for (int i = 0; i < nbParticles; i++) {
    sumDensity += density[i];
  }

  mass = nbParticles * densityResting / sumDensity;
}

void Fluid::calculateDensity() {
  double phi;
  double fourPih2 =
      (4.0 / (M_PI * radiusOfInfluence *
              radiusOfInfluence));  // Precalculated value used to avoid
                                    // multiple divisions and multiplications
  double hInverse = 1.0 / radiusOfInfluence;  // Precalculated value used to
                                              // avoid multiple divisions

  // find φ
  for (int i = 0; i < nbParticles; i++) {
    density[i] = 0;

    for (int j = 0; j < nbParticles; j++) {
      distanceQ[i * nbParticles + j] =
          std::abs(distance[i * nbParticles + j] * hInverse);

      if (distanceQ[i * nbParticles + j] < 1.0) {
        phi = fourPih2 *
              (1.0 - distanceQ[i * nbParticles + j] *
                         distanceQ[i * nbParticles + j]) *
              (1.0 - distanceQ[i * nbParticles + j] *
                         distanceQ[i * nbParticles + j]) *
              (1.0 -
               distanceQ[i * nbParticles + j] * distanceQ[i * nbParticles + j]);

      }

      else {
        phi = 0.0;
      }

      density[i] += mass * phi;
    }
  }
}

void Fluid::calculatePressure() {
  for (int i = 0; i < nbParticles; i++) {
    pressure[i] = gasConstant * (density[i] - densityResting);
  }
}
