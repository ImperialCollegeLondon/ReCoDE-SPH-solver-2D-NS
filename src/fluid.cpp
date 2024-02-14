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
  double phi, distanceQcurr, distanceQcurr3;

  // find φ
  for (size_t i = 0; i < nbParticles; i++) {
    density[i] = 0.0;

    for (size_t j = 0; j < nbParticles; j++) {
      distanceQcurr = distanceQ[i * nbParticles + j] =
          std::abs(distance[i * nbParticles + j] * hInverse);

      distanceQcurr3 = (1.0 - distanceQcurr * distanceQcurr);
      if (distanceQcurr < 1.0) {
        phi = fourPih2 * distanceQcurr3 * distanceQcurr3 * distanceQcurr3;

      }
      else {
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
