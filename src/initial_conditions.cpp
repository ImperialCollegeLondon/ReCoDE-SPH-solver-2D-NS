#include "initial_conditions.h"

#include <cmath>
#include <iostream>

#include "fluid.h"

// ========== Initial Conditions ==========

void icBasic(std::unique_ptr<Fluid> &fluidPtr, unsigned int nbParticles,
             std::vector<double> &positionX, std::vector<double> &positionY) {
  // Instead of using `new` operator, we're creating a
  // std::unique_ptr that can manage its resources (which
  // means that we don't need to care about
  // deleting the memory.). std::make_unique is using the
  // constructor of Fluid to create an std::unique_ptr.
  // See also RAII for an explanation on smart pointers.
  fluidPtr = std::make_unique<Fluid>(nbParticles);

  for (size_t i = 0; i < nbParticles; i++) {
    fluidPtr->setPositionX(i, positionX[i]);
    fluidPtr->setPositionY(i, positionY[i]);
    fluidPtr->setVelocityX(i, 0.0);
    fluidPtr->setVelocityY(i, 0.0);
  }
}

// Block drop
void icBlockDrop(std::unique_ptr<Fluid> &fluidPtr, unsigned int nbParticles,
                 double length, double width, double centerX, double centerY) {
  int n1, n2;

  nbParticles = rectangleN(nbParticles, length, width, n1, n2);

  fluidPtr = std::make_unique<Fluid>(nbParticles);

  // Distance between neighboring particles in x and y
  double dx = length / static_cast<double>((n1 - 1));
  double dy = width / static_cast<double>((n2 - 1));

  // Starting position in x
  double positionX = centerX - length / 2.0;
  double positionY;
  size_t particleIndex;  // indices

  // Assign the values in x for all particles
  for (size_t i = 0; i < n1; i++) {
    for (size_t j = 0; j < n2; j++) {
      particleIndex = i * n2 + j;

      fluidPtr->setPositionX(
          particleIndex,
          positionX + static_cast<double>(rand()) / RAND_MAX / 100000.0);
      fluidPtr->setVelocityX(particleIndex, 0.0);
    }
    positionX += dx;
  }

  // Assign the values in y for all particles
  for (size_t i = 0; i < n1; i++) {
    positionY = centerY - width / 2.0;
    for (size_t j = 0; j < n2; j++) {
      particleIndex = i * n2 + j;
      fluidPtr->setPositionY(
          particleIndex,
          positionY + static_cast<double>(rand()) / RAND_MAX / 100000.0);
      fluidPtr->setVelocityY(particleIndex, 0.0);
      positionY += dy;
    }
  }
}

// Droplet
void icDroplet(std::unique_ptr<Fluid> &fluidPtr, unsigned int nbParticles,
               double radius, double centerX, double centerY) {
  nbParticles = closestIntegerSqrt(nbParticles);

  std::vector<double> positionXStore(nbParticles, 0.0);
  std::vector<double> positionYStore(nbParticles, 0.0);

  unsigned int el = std::sqrt(nbParticles);

  // For uniform distribution the step in y has to be equal to the step in x
  double step = 2 * radius / (el - 1);
  double positionX = centerX - radius;  // Starting position in x
  double positionY;
  size_t positionIndex;

  for (size_t i = 0; i < el; i++) {
    positionY = centerY - radius;
    for (size_t j = 0; j < el; j++) {
      positionIndex = i * el + j;
      positionXStore[positionIndex] = positionX;
      positionYStore[positionIndex] = positionY;
      positionY += step;
    }
    positionX += step;
  }
  // After the initial square is created, the number of particles that are in
  // that square and from a distance from the centre less or equal to the radius
  // of the circle is calculated
  size_t count = 0;
  for (size_t i = 0; i < el; i++) {
    for (size_t j = 0; j < el; j++) {
      positionIndex = i * el + j;
      if (std::hypot(positionYStore[positionIndex] - centerY,
                     positionXStore[positionIndex] - centerX) <= radius) {
        count++;
      }
    }
  }

  nbParticles = count;

  fluidPtr = std::make_unique<Fluid>(nbParticles);

  size_t particleIndex = 0;
  for (size_t i = 0; i < el; i++) {
    for (size_t j = 0; j < el; j++) {
      positionIndex = i * el + j;
      if (std::hypot(positionYStore[positionIndex] - centerY,
                     positionXStore[positionIndex] - centerX) <= radius) {
        fluidPtr->setPositionX(particleIndex, positionXStore[positionIndex] +
                                                  static_cast<double>(rand()) /
                                                      RAND_MAX / 100000.0);
        fluidPtr->setPositionY(particleIndex, positionYStore[positionIndex] +
                                                  static_cast<double>(rand()) /
                                                      RAND_MAX / 100000.0);
        fluidPtr->setVelocityX(particleIndex, 0.0);
        fluidPtr->setVelocityY(particleIndex, 0.0);
        particleIndex++;
      }
    }
  }
}

unsigned int rectangleN(unsigned int num, double length, double width, int &n1,
                        int &n2) {
  // Function that transforms the user's particle related input to the closest
  // values that can be use to create a rectangle block
  double division = length / width;
  n2 = std::sqrt(num / division);

  n1 = ceil(division * n2);
  n2 = ceil(n2);
  // New number of particles
  return n1 * n2;
}

unsigned int closestIntegerSqrt(unsigned int num) {
  // Function that returns the closest number that has an integer square root
  double root = std::sqrt(num);

  auto integerRoot = static_cast<unsigned int>(root + 0.5);

  return integerRoot * integerRoot;
}