#include "sph_solver.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "fluid.h"
#include "main_prog_funcs.h"

// Setter functions
void SphSolver::setTimestep(double dt) { this->dt = dt; }

void SphSolver::setTotalIterations(double totalIterations) {
  this->totalIterations = totalIterations;
}

void SphSolver::setOutputFrequency(double f) { this->outputFrequency = f; }

void SphSolver::setCoeffRestitution(double coeffRestitution) {
  this->coeffRestitution = coeffRestitution;
}

void SphSolver::setLeftWall(double leftWall) { this->leftWall = leftWall; }

void SphSolver::setRightWall(double rightWall) { this->rightWall = rightWall; }

void SphSolver::setBottomWall(double bottomWall) {
  this->bottomWall = bottomWall;
}

void SphSolver::setTopWall(double topWall) { this->topWall = topWall; }

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
  cells.resize(numberOfCells);
  neighbourCells.resize(numberOfCells);

  assignNeighbourCells(cellsRows, cellsCols);
}

void SphSolver::assignNeighbourCells(int cellsRows, int cellsCols) {
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

void SphSolver::placeParticlesInCells(Fluid &data) {
  for (int i = 0; i < numberOfCells; i++) {
    cells[i].clear();
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

void SphSolver::neighbourParticlesSearch(Fluid &data) {
  for (int i = 0; i < numberOfParticles; i++) {
    neighbourParticles[i].clear();
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

          if (distance <= data.getRadInfl())
            neighbourParticles[cells[i][j]].push_back({cells[i][k], distance});
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
          if (distance <= data.getRadInfl())
            neighbourParticles[cells[i][j]].push_back(
                {cells[neighbourCells[i][k]][q], distance});
        }
      }
    }
  }
}

// Time integration
void SphSolver::timeIntegration(Fluid &data, std::ofstream &finalPositionsFile,
                                std::ofstream &energiesFile) {
  std ::cout << "Time integration started -- OK"
             << "\n";

  for (int time = 0; time < totalIterations; time++) {
    t = time;
    // In each iteration the distances between the particles are recalculated,
    // as well as their density and pressure
    neighbourParticlesSearch(data);
    data.calculateDensity(neighbourParticles);
    data.calculatePressure();
    particleIterations(data);

    if (time % outputFrequency == 0) {
      storeToFile(data, "energy", energiesFile, dt, t);
    }
  }
  // Store particles' positions after integration is completed
  storeToFile(data, "position", finalPositionsFile, dt, totalIterations);

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
  double radiusOfInfluence = data.getRadInfl();
  double thirtyPih3 =
      (-30.0 / (M_PI * radiusOfInfluence * radiusOfInfluence *
                radiusOfInfluence));  // Precalculated value used to avoid
                                      // multiple divisions and multiplications
  // For each neighbour particle, calculate the pressure force
  for (int j = 0; j < neighbourParticles[particleIndex].size(); j++) {
    normalisedDistance =
        neighbourParticles[particleIndex][j].second / radiusOfInfluence;
    sum += (data.getMass() /
            data.getDensity(neighbourParticles[particleIndex][j].first)) *
           ((data.getPressure(particleIndex) +
             data.getPressure(neighbourParticles[particleIndex][j].first)) /
            2.0) *
           (thirtyPih3 *
            (getPosition(particleIndex) -
             getPosition(neighbourParticles[particleIndex][j].first))) *
           (((1.0 - normalisedDistance) * (1.0 - normalisedDistance)) /
            normalisedDistance);
  }
  return -sum;
}

double SphSolver::calcViscousForce(Fluid &data,
                                   std::function<double(int)> getVelocity,
                                   int particleIndex) {
  double radiusOfInfluence = data.getRadInfl();

  double sum = 0.0;  // Initializing the summation
  double fourtyPih4 =
      (40.0 /
       (M_PI * radiusOfInfluence * radiusOfInfluence * radiusOfInfluence *
        radiusOfInfluence));  // Precalculated value used to avoid
                              // multiple divisions and multiplications
  // For each neighbour particle, calculate the viscous force
  for (int j = 0; j < neighbourParticles[particleIndex].size(); j++) {
    sum += (data.getMass() /
            data.getDensity(neighbourParticles[particleIndex][j].first)) *
           (getVelocity(particleIndex) -
            getVelocity(neighbourParticles[particleIndex][j].first)) *
           (fourtyPih4 * (1.0 - neighbourParticles[particleIndex][j].second /
                                    radiusOfInfluence));
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

  // y-direction
  newVelocity = data.getVelocityY(particleIndex) +
                integrationCoeff *
                    velocityIntegration(data, particleIndex, forcePressureY,
                                        forceViscousY, forceGravityY);

  data.setVelocityY(particleIndex, newVelocity);

  newPosition = data.getPositionY(particleIndex) + newVelocity * dt;
  data.setPositionY(particleIndex, newPosition);
}

double SphSolver::velocityIntegration(Fluid &data, int particleIndex,
                                      double &forcePressure,
                                      double &forceViscous,
                                      double &forceGravity) {
  double acceleration;
  acceleration = (forcePressure + forceViscous + forceGravity) /
                 data.getDensity(particleIndex);

  return acceleration * dt;
}

void SphSolver::boundaries(Fluid &data, int particleIndex) {
  // x-direction
  if (data.getPositionX(particleIndex) < leftWall + data.getRadInfl()) {
    data.setPositionX(particleIndex, leftWall + data.getRadInfl());
    data.setVelocityX(particleIndex,
                      -coeffRestitution * data.getVelocityX(particleIndex));
  }

  if (data.getPositionX(particleIndex) > rightWall - data.getRadInfl()) {
    data.setPositionX(particleIndex, rightWall - data.getRadInfl());
    data.setVelocityX(particleIndex,
                      -coeffRestitution * data.getVelocityX(particleIndex));
  }

  // y-direction
  if (data.getPositionY(particleIndex) < bottomWall + data.getRadInfl()) {
    data.setPositionY(particleIndex, bottomWall + data.getRadInfl());
    data.setVelocityY(particleIndex,
                      -coeffRestitution * data.getVelocityY(particleIndex));
  }

  if (data.getPositionY(particleIndex) > topWall - data.getRadInfl()) {
    data.setPositionY(particleIndex, topWall - data.getRadInfl());
    data.setVelocityY(particleIndex,
                      -coeffRestitution * data.getVelocityY(particleIndex));
  }
}
