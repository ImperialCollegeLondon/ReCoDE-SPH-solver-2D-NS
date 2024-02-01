#include "particles.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

// User defined constructor
particles::particles(const unsigned nNew) : nbParticles(nNew) {
  positionX = new double[nbParticles];
  positionY = new double[nbParticles];
  velocityX = new double[nbParticles];
  velocityY = new double[nbParticles];

  distance = new double[nbParticles * nbParticles];
  distanceQ = new double[nbParticles * nbParticles];

  particleSpeedSq = new double[nbParticles];
}

// Destructor
particles::~particles() {
  delete[] positionX;
  delete[] particleSpeedSq;
  delete[] positionY;
  delete[] velocityX;
  delete[] velocityY;
  delete[] distance;
  delete[] distanceQ;
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
    delete[] distance;
    delete[] distanceQ;
    delete[] particleSpeedSq;
    delete[] positionX;
    delete[] positionY;
    delete[] velocityX;
    delete[] velocityY;

    nbParticles = particles.nbParticles;

    positionX = new double[nbParticles];
    positionY = new double[nbParticles];
    velocityX = new double[nbParticles];
    velocityY = new double[nbParticles];

    distance = new double[nbParticles * nbParticles];
    distanceQ = new double[nbParticles * nbParticles];

    particleSpeedSq = new double[nbParticles];

    std::memcpy(positionX, particles.positionX, nbParticles * sizeof(double));
    std::memcpy(positionY, particles.positionY, nbParticles * sizeof(double));
    std::memcpy(velocityX, particles.velocityX, nbParticles * sizeof(double));
    std::memcpy(velocityY, particles.velocityY, nbParticles * sizeof(double));

    std::memcpy(distance, particles.distance,
                nbParticles * nbParticles * sizeof(double));
    std::memcpy(distanceQ, particles.distanceQ,
                nbParticles * nbParticles * sizeof(double));
    std::memcpy(particleSpeedSq, particles.particleSpeedSq,
                nbParticles * sizeof(double));
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
