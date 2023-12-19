#include "initial_conditions.h"
#include "sph.h"
#include <cmath>

// ========== Initial Conditions ==========

void ic_one_particle(int nb_particles, SPH &sph) {

  sph(0, 0) = 0.5;
  sph(1, 0) = 0.5;
  sph(2, 0) = 0.0;
  sph(3, 0) = 0.0;
}

void ic_two_particles(int nb_particles, SPH &sph) {

  sph(0, 0) = 0.5;
  sph(0, 1) = 0.5;

  sph(1, 0) = 0.5;
  sph(1, 1) = 0.01;

  sph(2, 0) = 0.0;
  sph(2, 1) = 0.0;

  sph(3, 0) = 0.0;
  sph(3, 1) = 0.0;
}

void ic_three_particles(int nb_particles, SPH &sph) {

  sph(0, 0) = 0.5;
  sph(0, 1) = 0.495;
  sph(0, 2) = 0.505;

  sph(1, 0) = 0.5;
  sph(1, 1) = 0.01;
  sph(1, 2) = 0.01;

  sph(2, 0) = 0.0;
  sph(2, 1) = 0.0;
  sph(2, 2) = 0.0;

  sph(3, 0) = 0.0;
  sph(3, 1) = 0.0;
  sph(3, 2) = 0.0;
}

void ic_four_particles(int nb_particles, SPH &sph) {

  sph(0, 0) = 0.505;
  sph(0, 1) = 0.515;
  sph(0, 2) = 0.51;
  sph(0, 3) = 0.5;

  sph(1, 0) = 0.5;
  sph(1, 1) = 0.5;
  sph(1, 2) = 0.45;
  sph(1, 3) = 0.45;

  sph(2, 0) = 0.0;
  sph(2, 1) = 0.0;
  sph(2, 2) = 0.0;
  sph(2, 3) = 0.0;

  sph(3, 0) = 0.0;
  sph(3, 1) = 0.0;
  sph(3, 2) = 0.0;
  sph(3, 3) = 0.0;
}

void ic_dam_break(int nb_particles, SPH &sph) {

  int el = pow(nb_particles, 0.5);
  // Initial distance between the particles in both directions
  double step = 0.19 / (el - 1);
  // Starting position in x
  double position_x = 0.01;
  double position_y;
  // Assing the values in x for all particles
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      sph(0, i * el + j) = position_x + double(rand()) / RAND_MAX / 100000;
      sph(2, i * el + j) = 0.0;
    }
    position_x += step;
  }

  // For uniform distribution the step in y has to be equal to the step in x
  step = 0.19 / (el - 1);

  // Assing values in y for all particles
  for (int i = 0; i < el; i++) {
    position_y = 0.01;
    for (int j = 0; j < el; j++) {
      sph(1, i * el + j) = position_y + double(rand()) / RAND_MAX / 100000;
      sph(3, i * el + j) = 0.0;
      position_y += step;
    }
  }
}

void ic_block_drop(int nb_particles, int n1, int n2, SPH &sph) {

  // Distance between neighbouring particles in x and y
  // 0.2 is the total distance in x and 0.3 in y
  double dx = 0.2 / double((n1 - 1));
  double dy = 0.3 / double((n2 - 1));

  // Starting position in x
  double position_x = 0.1;
  double position_y;
  int kx, ky;

  // Assing the values in x for all particles
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      kx = i * n2 + j;
      sph(0, kx) = position_x + double(rand()) / RAND_MAX / 100000;
      sph(2, kx) = 0.0;
    }
    position_x += dx;
  }

  // Assing the values in y for all particles
  for (int i = 0; i < n1; i++) {
    position_y = 0.3;
    for (int j = 0; j < n2; j++) {
      ky = i * n2 + j;
      sph(1, ky) = position_y + double(rand()) / RAND_MAX / 100000;
      sph(3, ky) = 0.0;
      position_y += dy;
    }
  }
}

// Droplet
void ic_droplet(int nb_particles, SPH &sph) {

  double *position_x_i = new double[nb_particles];
  double *position_y_i = new double[nb_particles];
  int el = pow(nb_particles, 0.5);
  int kx;

  // For uniform distribution the step in y has to be equal to the step in x
  double step = 0.2 / (el - 1);
  double position_x = 0.4; // Starting position in x
  double position_y;       // Starting position in y

  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      position_x_i[i * el + j] = position_x;
    }
    position_x += step;
  }

  step = 0.2 / (el - 1);

  for (int i = 0; i < el; i++) {
    position_y = 0.6;
    for (int j = 0; j < el; j++) {
      position_y_i[i * el + j] = position_y;
      position_y += step;
    }
  }
  kx = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (sqrt(pow((position_y_i[i * el + j] - 0.7), 2) +
               pow((position_x_i[i * el + j] - 0.5), 2)) <= 0.1) {
        sph(0, kx) =
            position_x_i[i * el + j] + double(rand()) / RAND_MAX / 100000;
        sph(1, kx) =
            position_y_i[i * el + j] + double(rand()) / RAND_MAX / 100000;
        sph(2, kx) = 0;
        sph(3, kx) = 0;
        kx++;
      }
    }
  }
  delete[] position_x_i;
  delete[] position_y_i;
}

// Defines the number of particles that will be in the circular region
int dropletn(int nb_particles) {

  // Process similar to dam break. Creates an initial square
  double *position_x_i = new double[nb_particles];
  double *position_y_i = new double[nb_particles];
  int el = pow(nb_particles, 0.5);
  double step = 0.2 / (el - 1);
  double position_x = 0.4;
  double position_y;

  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      position_x_i[i * el + j] = position_x;
    }
    position_x += step;
  }

  step = 0.2 / (el - 1);

  for (int i = 0; i < el; i++) {
    position_y = 0.6;
    for (int j = 0; j < el; j++) {
      position_y_i[i * el + j] = position_y;
      position_y += step;
    }
  }

  // After the initial square is created, the number of particles that are in
  // that square and from a distance from the centre less or equal to the radius
  // of the circle is calculated
  int count = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (sqrt(pow((position_y_i[i * el + j] - 0.7), 2) +
               pow((position_x_i[i * el + j] - 0.5), 2)) <= 0.1) {
        count++;
      } else {
        count += 0;
      }
    }
  }
  delete[] position_x_i;
  delete[] position_y_i;

  return count;
}