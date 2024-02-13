#include "fluid.h"

#include <numeric>
#include <iostream>

// User defined constructor
Fluid::Fluid(const unsigned nNew)
    : particles(nNew), density(nNew, 0.0), pressure(nNew, 0.0) {}

// Calculation functions

void Fluid::calculateMass(
    std::vector<std::vector<std::pair<int, double>>> neighbours) {
  calculateDensity(neighbours);
  double sumDensity = std::accumulate(density.begin(), density.end(), 0.0);

  mass = nbParticles * densityResting / sumDensity;
}

void Fluid::calculateDensity(
    std::vector<std::vector<std::pair<int, double>>> neighbours) {
  double phi, normalisedDistance, normalisedDistanceSqr;
  // find φ
  for (size_t i = 0; i < nbParticles; i++) {
    density[i] = mass * fourPih2;
    for (size_t j = 0; j < neighbours[i].size(); j++) {
      normalisedDistance = neighbours[i][j].second * hInverse;
      normalisedDistanceSqr = (1.0 - normalisedDistance * normalisedDistance);
      phi = fourPih2 * normalisedDistanceSqr * normalisedDistanceSqr *
              normalisedDistanceSqr;
      density[i] += mass * phi;
    }
  }
}

void Fluid::calculatePressure() {
  for (size_t i = 0; i < nbParticles; i++) {
    pressure[i] = gasConstant * (density[i] - densityResting);
  }
}
