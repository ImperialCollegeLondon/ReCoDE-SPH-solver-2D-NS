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