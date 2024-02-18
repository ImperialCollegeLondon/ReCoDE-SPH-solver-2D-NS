#include "particles.h"

#include <cmath>
#include <iostream>

// User defined constructor
particles::particles(const unsigned nNew)
    : nbParticles(nNew),
      positionX(nNew, 0.0),
      positionY(nNew, 0.0),
      velocityX(nNew, 0.0),
      velocityY(nNew, 0.0),
      particleSpeedSq(nNew, 0.0),
      distance(nNew * nNew, 0.0),
      distanceQ(nNew * nNew, 0.0) {}

// Calculation functions

void particles::calculateParticleDistance() {
  double dx;
  double dy;

  for (size_t i = 0; i < nbParticles; i++) {
    for (size_t j = 0; j < nbParticles; j++) {
      dx = positionX[i] - positionX[j];
      dy = positionY[i] - positionY[j];

      distance[i * nbParticles + j] = sqrt(dx * dx + dy * dy);
    }
  }
}