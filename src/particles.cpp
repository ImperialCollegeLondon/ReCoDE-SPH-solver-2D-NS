#include "particles.h"

#include <cmath>
#include <iostream>

// User defined constructor
particles::particles(const unsigned nNew) : nbParticles(nNew) {
  positionX.reserve(nbParticles);
  positionY.reserve(nbParticles);
  velocityX.reserve(nbParticles);
  velocityY.reserve(nbParticles);

  particleSpeedSq.reserve(nbParticles);
}