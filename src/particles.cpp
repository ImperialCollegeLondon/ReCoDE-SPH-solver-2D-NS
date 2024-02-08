#include "particles.h"

#include <cmath>
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

// Getter functions

int particles::getNumberOfParticles() { return nbParticles; }

double particles::getPositionX(int k) { return positionX[k]; }

double particles::getPositionY(int k) { return positionY[k]; }

double particles::getVelocityX(int k) { return velocityX[k]; }

double particles::getVelocityY(int k) { return velocityY[k]; }

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

// Setter functions

void particles::setPositionX(int k, double newPositionX) {
  positionX[k] = newPositionX;
}

void particles::setPositionY(int k, double newPositionY) {
  positionY[k] = newPositionY;
}

void particles::setVelocityX(int k, double newVelocityX) {
  velocityX[k] = newVelocityX;
}

void particles::setVelocityY(int k, double newVelocityY) {
  velocityY[k] = newVelocityY;
}
