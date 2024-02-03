#include "particles.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// User defined constructor
particles::particles(const unsigned nNew) : nbParticles(nNew) {
  positionX.reserve(nbParticles);
  positionY.reserve(nbParticles);
  velocityX.reserve(nbParticles);
  velocityY.reserve(nbParticles);

  distance.reserve(nbParticles * nbParticles);
  distanceQ.reserve(nbParticles * nbParticles);

  particleSpeedSq.reserve(nbParticles);
}

// Overloading of ()
double &particles::operator()(unsigned row, unsigned col) {
  switch (row) {
    case 0:
      return this->positionX[col];
      break;
    case 1:
      return this->positionY[col];
      break;
    case 2:
      return this->velocityX[col];
      break;
    case 3:
      return this->velocityY[col];
      break;
    default:
      std::cerr << "ERROR: Out of bounds on row selection" << std::endl;
      abort();
  }
}

// Overloading of operator=
particles &particles::operator=(const particles &particles) {
  if (this != &particles) {

    nbParticles = particles.nbParticles;

    positionX.reserve(nbParticles);
    positionY.reserve(nbParticles);
    velocityX.reserve(nbParticles);
    velocityY.reserve(nbParticles);

    distance.reserve(nbParticles * nbParticles);
    distanceQ.reserve(nbParticles * nbParticles);

    particleSpeedSq.reserve(nbParticles);

    positionX = particles.positionX;
    positionY = particles.positionY;
    velocityX = particles.velocityX;
    velocityY = particles.velocityY;

    distance = particles.distance;
    distanceQ = particles.distanceQ;
    particleSpeedSq = particles.particleSpeedSq;
  }
  return *this;
}

// Getter functions

int particles::getNumberOfParticles() { return nbParticles; }

double particles::getPositionX(int k) { return positionX[k]; }

double particles::getPositionY(int k) { return positionY[k]; }

double particles::getDistanceQ(int k) { return distanceQ[k]; }

// Calculation functions

void particles::calculateParticleDistance() {
  double dx;
  double dy;

  for (int i = 0; i < nbParticles; i++) {
    for (int j = 0; j < nbParticles; j++) {
      dx = positionX[i] - positionX[j];
      dy = positionY[i] - positionY[j];

      distance[i * nbParticles + j] = sqrt(dx * dx + dy * dy);
    }
  }
}
