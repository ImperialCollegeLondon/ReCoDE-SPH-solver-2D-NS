#include "initial_conditions.h"

#include <cmath>
#include <iostream>

#include "fluid.h"

// ========== Initial Conditions ==========

void ic_basic(fluid *&fluid_ptr, int nb_particles, double *position_x,
              double *position_y) {
  // Allocate memory for the fluid object and call the constructor
  // This needs to be deleted by the caller.
  fluid_ptr = new fluid(nb_particles);

  fluid &fluid = *fluid_ptr;  // Use a reference to the object

  for (int i = 0; i < nb_particles; i++) {
    fluid(0, i) = position_x[i];
    fluid(1, i) = position_y[i];
    fluid(2, i) = 0.0;
    fluid(3, i) = 0.0;
  }

  return;
}

// Block drop
void ic_block_drop(fluid *&fluid_ptr, int &nb_particles, double length,
                   double width, double center_x, double center_y) {
  int n1, n2;
  nb_particles = rectangle_n(nb_particles, length, width, n1, n2);

  // Allocate memory for the fluid object and call the constructor
  // This needs to be deleted by the caller.
  fluid_ptr = new fluid(nb_particles);

  fluid &fluid = *fluid_ptr;  // Use a reference to the object

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
      fluid(0, kx) = position_x + double(rand()) / RAND_MAX / 100000;
      fluid(2, kx) = 0.0;
    }
    position_x += dx;
  }

  // Assign the values in y for all particles
  for (int i = 0; i < n1; i++) {
    position_y = center_y - width / 2.0;
    for (int j = 0; j < n2; j++) {
      ky = i * n2 + j;
      fluid(1, ky) = position_y + double(rand()) / RAND_MAX / 100000;
      fluid(3, ky) = 0.0;
      position_y += dy;
    }
  }

  return;
}

// Droplet
void ic_droplet(fluid *&fluid_ptr, int &nb_particles, double radius,
                double center_x, double center_y) {
  nb_particles = closest_integer_sqrt(nb_particles);

  double *position_x_store = new double[nb_particles];
  double *position_y_store = new double[nb_particles];
  int el = std::sqrt(nb_particles);

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

  nb_particles = count;
  // Allocate memory for the fluid object and call the constructor
  // This needs to be deleted by the caller.
  fluid_ptr = new fluid(nb_particles);

  fluid &fluid = *fluid_ptr;  // Use a reference to the object

  kx = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (std::hypot(position_y_store[i * el + j] - center_y,
                     position_x_store[i * el + j] - center_x) <= radius) {
        fluid(0, kx) =
            position_x_store[i * el + j] + double(rand()) / RAND_MAX / 100000;
        fluid(1, kx) =
            position_y_store[i * el + j] + double(rand()) / RAND_MAX / 100000;
        fluid(2, kx) = 0.0;
        fluid(3, kx) = 0.0;
        kx++;
      }
    }
  }

  delete[] position_x_store;
  delete[] position_y_store;

  return;
}

int rectangle_n(int nb_particles, double length, double width, int &n1,
                int &n2) {
  // Function that transforms the user's particle related input to the closest
  // values that can be use to create a rectangle block
  double division = length / width;
  n2 = std::sqrt(nb_particles / division);

  n1 = ceil(division * n2);
  n2 = ceil(n2);
  // New number of particles
  return n1 * n2;
}

int closest_integer_sqrt(int num) {
  // Function that returns the closest number that has an integer square root
  double root = std::sqrt(num);

  int integerRoot = static_cast<int>(root + 0.5);

  return integerRoot * integerRoot;
}