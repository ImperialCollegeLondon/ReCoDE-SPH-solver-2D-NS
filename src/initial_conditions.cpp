#include "initial_conditions.h"

#include <cmath>
#include <iostream>

#include "sph.h"

// ========== Initial Conditions ==========

SPH ic_basic(int nb_particles, double *position_x, double *position_y) {
  SPH sph(nb_particles);

  for (int i = 0; i < nb_particles; i++) {
    sph(0, i) = position_x[i];
    sph(1, i) = position_y[i];
    sph(2, i) = 0.0;
    sph(3, i) = 0.0;
  }

  return sph;
}

SPH ic_block_drop(int nb_particles, double length, double width,
                  double center_x, double center_y) {
  int n1, n2;
  nb_particles = rectangle_n(nb_particles, length, width, n1, n2);
  SPH sph(nb_particles);

  // Distance between neighboring particles in x and y
  double dx = length / double((n1 - 1));
  double dy = width / double((n2 - 1));

  // Starting position in x
  double position_x = center_x - length / 2.0;
  double position_y;
  int kx, ky;  // indices

  // Assign the values in x for all particles
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      kx = i * n2 + j;
      sph(0, kx) = position_x + double(rand()) / RAND_MAX / 100000;
      sph(2, kx) = 0.0;
    }
    position_x += dx;
  }

  // Assign the values in y for all particles
  for (int i = 0; i < n1; i++) {
    position_y = center_y - width / 2.0;
    for (int j = 0; j < n2; j++) {
      ky = i * n2 + j;
      sph(1, ky) = position_y + double(rand()) / RAND_MAX / 100000;
      sph(3, ky) = 0.0;
      position_y += dy;
    }
  }

  return sph;
}

// Droplet
SPH ic_droplet(int nb_particles, double radius, double center_x,
               double center_y) {
  nb_particles = closest_power_of_two(nb_particles);
  std::cout << "Number of particles: " << nb_particles << std::endl;
  double *position_x_store = new double[nb_particles];
  double *position_y_store = new double[nb_particles];
  int el = std::sqrt(nb_particles);
  std::cout << "Number of particles per side: " << el << std::endl;
  int kx;

  // For uniform distribution the step in y has to be equal to the step in x
  double step = 2 * radius / (el - 1);
  double position_x = center_x - radius;  // Starting position in x
  double position_y;

  for (int i = 0; i < el; i++) {
    position_y = center_y - radius;
    for (int j = 0; j < el; j++) {
      position_x_store[i * el + j] = position_x;
      position_y_store[i * el + j] = position_y;
      position_y += step;
    }
    position_x += step;
  }
  // After the initial square is created, the number of particles that are in
  // that square and from a distance from the centre less or equal to the radius
  // of the circle is calculated
  int count = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (std::hypot(position_y_store[i * el + j] - center_y,
                     position_x_store[i * el + j] - center_x) <= radius) {
        count++;
      }
    }
  }
  std::cout << "Number of particles in the circle: " << count << std::endl;

  SPH sph(count);
  kx = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (std::hypot(position_y_store[i * el + j] - center_y,
                     position_x_store[i * el + j] - center_x) <= radius) {
        sph(0, kx) =
            position_x_store[i * el + j] + double(rand()) / RAND_MAX / 100000;
        sph(1, kx) =
            position_y_store[i * el + j] + double(rand()) / RAND_MAX / 100000;
        sph(2, kx) = 0.0;
        sph(3, kx) = 0.0;
        kx++;
      }
    }
  }

  delete[] position_x_store;
  delete[] position_y_store;

  return sph;
}

int rectangle_n(int nb_particles, double length, double width, int &n1,
                int &n2) {
  double division = length / width;
  n2 = std::sqrt(nb_particles / division);

  n1 = ceil(division * n2);
  n2 = ceil(n2);

  return n1 * n2;
}

int closest_power_of_two(int nb_particles) {
  // Ensure n is not negative
  if (nb_particles <= 0) {
    return 0;  // or handle the case as needed
  }

  int result = 1;
  while (result < nb_particles) {
    result <<=
        1;  // Left shift by 1 to double the value (same as multiplying by 2)
  }

  // Check which power of 2 is closest to n
  int lowerPower =
      result >> 1;  // Right shift by 1 to get the previous power of 2

  return (result - nb_particles) < (nb_particles - lowerPower) ? result
                                                               : lowerPower;
}