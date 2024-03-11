#include "particles.h"

#include <cmath>
#include <iostream>

// User defined constructor
Particles::Particles(const unsigned nNew)
    : nbParticles(nNew),
      positionX(nNew, 0.0),
      positionY(nNew, 0.0),
      velocityX(nNew, 0.0),
      velocityY(nNew, 0.0),
      particleSpeedSq(nNew, 0.0) {}