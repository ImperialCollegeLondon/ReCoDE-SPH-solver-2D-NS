#include "fluid.h"

#include <numeric>

// User defined constructor
Fluid::Fluid(const unsigned nNew)
    : particles(nNew), density(nNew, 0.0), pressure(nNew, 0.0) {}

// Calculation functions

void Fluid::calculateMass() {
  calculateParticleDistance();
  calculateDensity();
  double sumDensity = std::accumulate(density.begin(), density.end(), 0.0);

  mass = nbParticles * densityResting / sumDensity;
}

void Fluid::calculateDensity() {
  double phi, normalisedDistance, normalisedDistanceSqr;

  // find φ
  for (size_t i = 0; i < nbParticles; i++) {
    density[i] = 0.0;

    for (size_t j = 0; j < nbParticles; j++) {
      normalisedDistance = distanceQ[i * nbParticles + j] =
          std::abs(distance[i * nbParticles + j] * hInverse);

      normalisedDistanceSqr = (1.0 - normalisedDistance * normalisedDistance);
      if (normalisedDistance < 1.0) {
        phi = fourPih2 * normalisedDistanceSqr * normalisedDistanceSqr *
              normalisedDistanceSqr;

      } else {
        phi = 0.0;
      }

      density[i] += mass * phi;
    }
  }
}

void Fluid::calculatePressure() {
  for (size_t i = 0; i < nbParticles; i++) {
    pressure[i] = gasConstant * (density[i] - densityResting);
  }
}
