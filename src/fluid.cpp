#include "fluid.h"

#include <cmath>

// User defined constructor
Fluid::Fluid(const unsigned nNew) : particles(nNew) {
  pressure.reserve(nbParticles);
  density.reserve(nbParticles);
}
// Copy constructor
Fluid &Fluid::operator=(const Fluid &fluid) {
  if (this != &fluid) {
    nbParticles = fluid.nbParticles;

    positionX.reserve(nbParticles);
    positionY.reserve(nbParticles);
    velocityX.reserve(nbParticles);
    velocityY.reserve(nbParticles);

    distance.reserve(nbParticles * nbParticles);
    distanceQ.reserve(nbParticles * nbParticles);

    particleSpeedSq.reserve(nbParticles);

    pressure.reserve(nbParticles);
    density.reserve(nbParticles);

    positionX = fluid.positionX;
    positionY = fluid.positionY;
    velocityX = fluid.velocityX;
    velocityY = fluid.velocityY;

    distance = fluid.distance;
    distanceQ = fluid.distanceQ;
    particleSpeedSq = fluid.particleSpeedSq;

    pressure = fluid.pressure;
    density = fluid.density;
  }
  return *this;
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
