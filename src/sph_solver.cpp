#include "sph_solver.h"

#include <iomanip>
#include <iostream>

#include "fluid.h"
#include "main_prog_funcs.h"


// Getter functions
std::vector<std::vector<std::pair<int, double>>>
SphSolver::getNeighbourParticles() {
  return neighbourParticles;
}

// Neighbour search functions

// Split the domain into cells
void SphSolver::createGrid(Fluid &data) {
  numberOfParticles = data.getNumberOfParticles();
  neighbourParticles.resize(numberOfParticles);

  double radiusOfInfluence = data.getRadInfl();
  int cellsRows =
      static_cast<int>(std::ceil((topWall - bottomWall) / radiusOfInfluence));
  int cellsCols =
      static_cast<int>(std::ceil((rightWall - leftWall) / radiusOfInfluence));
  numberOfCells = cellsRows * cellsCols;
  cells.reserve(numberOfCells);
  neighbourCells.reserve(numberOfCells);

  assignNeighbourCells(cellsRows, cellsCols);
}

void SphSolver::assignNeighbourCells(int cellsRows, int cellsCols) {
  // Each cell could have at most 8 neighbours (and most of them do), so reserve
  // the memory
  for (int i = 0; i < numberOfCells; i++) {
    neighbourCells[i].reserve(MAX_NEIGHBOUR_CELLS);
  }
  // Flags to check if the cell is on the edge or in the middle
  bool top = false;
  bool left = false;
  bool right = false;
  bool bottom = false;
  bool middle = false;

  for (int i = 0; i < numberOfCells; i++) {
    // Cell has a bottom neighbour
    if (i >= cellsCols) {
      neighbourCells[i].push_back(i - cellsCols);
      bottom = true;
    }
    // Cell has a top neighbour
    if (i < numberOfCells - cellsCols) {
      neighbourCells[i].push_back(i + cellsCols);
      top = true;
    }
    // Cell has a left neighbour
    if (i % cellsCols != 0) {
      neighbourCells[i].push_back(i - 1);
      left = true;
    }
    // Cell has a right neighbour
    if ((i + 1) % cellsCols != 0) {
      neighbourCells[i].push_back(i + 1);
      right = true;
    }
    // If the cell is not on the edge, add the diagonal neighbours
    if (bottom && top && left && right) {
      neighbourCells[i].push_back(i - 1 - cellsCols);
      neighbourCells[i].push_back(i - 1 + cellsCols);
      neighbourCells[i].push_back(i + 1 - cellsCols);
      neighbourCells[i].push_back(i + 1 + cellsCols);
      middle = true;
    }
    // If the cell is on the edge, add only specific neighbours
    if (!middle) {
      // Add bottom-left diagonal neighbour
      if (bottom && left) {
        neighbourCells[i].push_back(i - 1 - cellsCols);
      }
      // Add bottom-right diagonal neighbour
      if (bottom && right) {
        neighbourCells[i].push_back(i + 1 - cellsCols);
      }
      // Add top-left diagonal neighbour
      if (top && left) {
        neighbourCells[i].push_back(i - 1 + cellsCols);
      }
      // Add top-right diagonal neighbour
      if (top && right) {
        neighbourCells[i].push_back(i + 1 + cellsCols);
      }
    }
    // Reset the flags
    top = false;
    left = false;
    right = false;
    bottom = false;
    middle = false;
  }
}

void SphSolver::neighbourParticlesSearch(Fluid &data) {
  int currentNumberOfNeighbours;
  for (int i = 0; i < numberOfParticles; i++) {
    currentNumberOfNeighbours = neighbourParticles[i].size();
    neighbourParticles[i].clear();
    neighbourParticles[i].reserve(
        static_cast<int>(memoryReservationFactor * currentNumberOfNeighbours));
  }

  placeParticlesInCells(data);

  // For each cell, for each particle in the cell, find neighbour particles in
  // the cell
  for (int i = 0; i < numberOfCells; i++) {
    for (int j = 0; j < cells[i].size(); j++) {
      for (int k = 0; k < cells[i].size(); k++) {
        if (cells[i][j] != cells[i][k]) {
          double distance = sqrt(pow(data.getPositionX(cells[i][j]) -
                                         data.getPositionX(cells[i][k]),
                                     2) +
                                 pow(data.getPositionY(cells[i][j]) -
                                         data.getPositionY(cells[i][k]),
                                     2));

          if (distance <= data.getRadInfl()) {
            neighbourParticles[cells[i][j]].push_back({cells[i][k], distance});
          }
        }
      }
    }
  }

  // For each cell, for each particle in the cell, find neighbour particles in
  // all neighbour cells
  for (int i = 0; i < numberOfCells; i++) {
    for (int j = 0; j < cells[i].size(); j++) {
      for (int k = 0; k < neighbourCells[i].size(); k++) {
        for (int q = 0; q < cells[neighbourCells[i][k]].size(); q++) {
          double distance =
              sqrt(pow(data.getPositionX(cells[i][j]) -
                           data.getPositionX(cells[neighbourCells[i][k]][q]),
                       2) +
                   pow(data.getPositionY(cells[i][j]) -
                           data.getPositionY(cells[neighbourCells[i][k]][q]),
                       2));
          if (distance <= data.getRadInfl()) {
            neighbourParticles[cells[i][j]].push_back(
                {cells[neighbourCells[i][k]][q], distance});
          }
        }
      }
    }
  }
}

void SphSolver::placeParticlesInCells(Fluid &data) {
  int currentCellSize;
  for (int i = 0; i < numberOfCells; i++) {
    currentCellSize = cells[i].size();
    cells[i].clear();
    cells[i].reserve(
        static_cast<int>(memoryReservationFactor * currentCellSize));
  }

  double radiusOfInfluence = data.getRadInfl();
  int cellsCols =
      static_cast<int>(std::ceil((rightWall - leftWall) / radiusOfInfluence));
  for (int i = 0; i < numberOfParticles; i++) {
    double positionX = data.getPositionX(i);
    double positionY = data.getPositionY(i);
    int j = static_cast<int>(positionX / radiusOfInfluence) +
            static_cast<int>(positionY / radiusOfInfluence) * cellsCols;
    cells[j].push_back(i);
  }
}

// Time integration
void SphSolver::timeIntegration(Fluid &data, std::ofstream &finalPositionsFile,
                                std::ofstream &energiesFile) {
  std ::cout << "Time integration started -- OK"
             << "\n";

  while (timeInteg < totalTime) {
    if (adaptiveTimestepBool) {
      vmax = 0.0;
      amax = 0.0;
      if (t == 0) {
        dt = 1e-8;
      }
    }

    // In each iteration the distances between the particles are recalculated,
    // as well as their density and pressure
    neighbourParticlesSearch(data);
    data.calculateDensity(neighbourParticles);
    data.calculatePressure();
    particleIterations(data);

    if (t % outputFrequency == 0) {
      storeToFile(data, "energy", energiesFile, dt, timeInteg);
    }

    timeInteg += dt;
    t++;

    if (adaptiveTimestepBool) {
      adaptiveTimestep(data);
    }
  }
  // Store particles' positions after integration is completed
  storeToFile(data, "position", finalPositionsFile, dt, timeInteg);

  std ::cout << "Time integration finished -- OK"
             << "\n";
}

void SphSolver::particleIterations(Fluid &data) {
  int i;

  // Use std::function to store the member functions
  std::function<double(int)> ptrGetPositionX =
      std::bind(&Fluid::getPositionX, &data, std::placeholders::_1);
  std::function<double(int)> ptrGetPositionY =
      std::bind(&Fluid::getPositionY, &data, std::placeholders::_1);
  std::function<double(int)> ptrGetVelocityX =
      std::bind(&Fluid::getVelocityX, &data, std::placeholders::_1);
  std::function<double(int)> ptrGetVelocityY =
      std::bind(&Fluid::getVelocityY, &data, std::placeholders::_1);

  for (i = 0; i < numberOfParticles; i++) {
    // Gathering the forces calculated by the processors
    forcePressureX = calculatePressureForce(data, ptrGetPositionX, i);

    forcePressureY = calculatePressureForce(data, ptrGetPositionY, i);

    forceViscousX = calcViscousForce(data, ptrGetVelocityX, i);

    forceViscousY = calcViscousForce(data, ptrGetVelocityY, i);

    forceGravityY = calcGravityForce(data, i);

    // Update the position of the particle
    updatePosition(data, i);

    // Boundary Conditions
    boundaries(data, i);
  }
}

double SphSolver::calculatePressureForce(Fluid &data,
                                         std::function<double(int)> getPosition,
                                         int particleIndex) {
  double sum = 0.0;  // Initializing the summation
  double normalisedDistance;
  double position = getPosition(particleIndex);
  double pressure = data.getPressure(particleIndex);
  double mass = data.getMass();
  double radiusOfInfluence = data.getRadInfl();

  for (int j = 0; j < neighbourParticles[particleIndex].size(); j++) {
    if (particleIndex != neighbourParticles[particleIndex][j].first) {
      normalisedDistance =
          neighbourParticles[particleIndex][j].second / radiusOfInfluence;

      sum +=
          (mass / data.getDensity(neighbourParticles[particleIndex][j].first)) *
          ((pressure +
            data.getPressure(neighbourParticles[particleIndex][j].first)) /
           2.0) *
          (thirtyPih3 *
           (position -
            getPosition(neighbourParticles[particleIndex][j].first))) *
          (((1.0 - normalisedDistance) * (1.0 - normalisedDistance)) /
           normalisedDistance);
    }
  }
  return -sum;
}

double SphSolver::calcViscousForce(Fluid &data,
                                   std::function<double(int)> getVelocity,
                                   int particleIndex) {
  double sum = 0.0;  // Initializing the summation
  double velocity = getVelocity(particleIndex);
  double mass = data.getMass();
  double radiusOfInfluence = data.getRadInfl();

  for (int j = 0; j < neighbourParticles[particleIndex].size(); j++) {
    if (particleIndex != neighbourParticles[particleIndex][j].first) {
      sum +=
          (mass / data.getDensity(neighbourParticles[particleIndex][j].first)) *
          (velocity - getVelocity(neighbourParticles[particleIndex][j].first)) *
          (fourtyPih4 * (1.0 - neighbourParticles[particleIndex][j].second /
                                   radiusOfInfluence));
    }
  }

  return -data.getViscosity() * sum;
}

double SphSolver::calcGravityForce(Fluid &data, int particleIndex) {
  return -data.getDensity(particleIndex) * data.getAccelerationGravity();
}

void SphSolver::updatePosition(Fluid &data, int particleIndex) {
  double newVelocity;
  double newPosition;
  double integrationCoeff = 1.0;

  // First step to initialise the scheme
  if (t == 0) {
    integrationCoeff = 0.5;
  }

  // x-direction
  newVelocity = data.getVelocityX(particleIndex) +
                integrationCoeff *
                    velocityIntegration(data, particleIndex, forcePressureX,
                                        forceViscousX, forceGravityX);
  data.setVelocityX(particleIndex, newVelocity);

  newPosition = data.getPositionX(particleIndex) + newVelocity * dt;

  data.setPositionX(particleIndex, newPosition);

  if (adaptiveTimestepBool) {
    vmax = std::max(vmax, std::abs(newVelocity));
  }

  // y-direction
  newVelocity = data.getVelocityY(particleIndex) +
                integrationCoeff *
                    velocityIntegration(data, particleIndex, forcePressureY,
                                        forceViscousY, forceGravityY);

  data.setVelocityY(particleIndex, newVelocity);

  newPosition = data.getPositionY(particleIndex) + newVelocity * dt;
  data.setPositionY(particleIndex, newPosition);

  if (adaptiveTimestepBool) {
    vmax = std::max(vmax, std::abs(newVelocity));
  }
}

double SphSolver::velocityIntegration(Fluid &data, int particleIndex,
                                      double forcePressure, double forceViscous,
                                      double forceGravity) {
  double acceleration = (forcePressure + forceViscous + forceGravity) /
                        data.getDensity(particleIndex);

  if (adaptiveTimestepBool) {
    amax = std::max(amax, std::abs(acceleration));
  }

  return acceleration * dt;
}

void SphSolver::boundaries(Fluid &data, int particleIndex) {
  // x-direction
  if (data.getPositionX(particleIndex) < leftWall + data.getRadInfl()) {
    data.setPositionX(particleIndex, leftWall + data.getRadInfl());

    data.setVelocityX(particleIndex,
                      -coeffRestitution * data.getVelocityX(particleIndex));

  } else if (data.getPositionX(particleIndex) > rightWall - data.getRadInfl()) {
    data.setPositionX(particleIndex, rightWall - data.getRadInfl());

    data.setVelocityX(particleIndex,
                      -coeffRestitution * data.getVelocityX(particleIndex));
  }

  // y-direction
  if (data.getPositionY(particleIndex) < bottomWall + data.getRadInfl()) {
    data.setPositionY(particleIndex, bottomWall + data.getRadInfl());

    data.setVelocityY(particleIndex,
                      -coeffRestitution * data.getVelocityY(particleIndex));

  } else if (data.getPositionY(particleIndex) > topWall - data.getRadInfl()) {
    data.setPositionY(particleIndex, topWall - data.getRadInfl());

    data.setVelocityY(particleIndex,
                      -coeffRestitution * data.getVelocityY(particleIndex));
  }
}

void SphSolver::adaptiveTimestep(Fluid &data) {
  double h = data.getRadInfl();

  dt = std::min(0.025 * h / vmax, 0.05 * pow(h / amax, 0.5));
}
