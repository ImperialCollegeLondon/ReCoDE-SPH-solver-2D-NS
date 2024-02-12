#include "fluid.h"

#include <cmath>
#include <iostream>

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

    particleSpeedSq.reserve(nbParticles);

    pressure.reserve(nbParticles);
    density.reserve(nbParticles);

    positionX = fluid.positionX;
    positionY = fluid.positionY;
    velocityX = fluid.velocityX;
    velocityY = fluid.velocityY;

    particleSpeedSq = fluid.particleSpeedSq;

    pressure = fluid.pressure;
    density = fluid.density;
  }
  return *this;
}

// Calculation functions

void Fluid::calculateMass(
    std::vector<std::vector<std::pair<int, double>>> neighbours) {
  calculateDensity(neighbours);
  double sumDensity = 0.0;
  for (int i = 0; i < nbParticles; i++) {
    sumDensity += density[i];
  }

  mass = nbParticles * densityResting / sumDensity;
  std::cout << "Mass = " << mass << std::endl;
}

void Fluid::calculateDensity(
    std::vector<std::vector<std::pair<int, double>>> neighbours) {
  double phi;
  double fourPih2 =
      (4.0 / (M_PI * radiusOfInfluence *
              radiusOfInfluence));  // Precalculated value used to avoid
                                    // multiple divisions and multiplications
  double hInverse = 1.0 / radiusOfInfluence;  // Precalculated value used to
                                              // avoid multiple divisions
  double normalisedDistance;
  // find φ
  for (int i = 0; i < nbParticles; i++) {

  density[i] = mass * fourPih2;
    for (int j = 0; j < neighbours[i].size(); j++) {
      normalisedDistance = neighbours[i][j].second * hInverse;
      phi = fourPih2 * (1.0 - normalisedDistance * normalisedDistance) *
            (1.0 - normalisedDistance * normalisedDistance) *
            (1.0 - normalisedDistance * normalisedDistance);
      density[i] += mass * phi;
    }
  }
}

void Fluid::calculatePressure() {
  for (int i = 0; i < nbParticles; i++) {
    pressure[i] = gasConstant * (density[i] - densityResting);
  }
}
