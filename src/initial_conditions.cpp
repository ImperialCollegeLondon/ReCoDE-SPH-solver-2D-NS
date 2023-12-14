#include "class.h"
#include "initial_conditions.h"
#include <cmath>

// ========== Initial Conditions ==========

void ic_one_particle(int n, SPH &sph) {

  sph(0, 0) = 0.5;
  sph(1, 0) = 0.5;
  sph(2, 0) = 0.0;
  sph(3, 0) = 0.0;
}

void ic_two_particles(int n, SPH &sph) {

  sph(0, 0) = 0.5;
  sph(0, 1) = 0.5;

  sph(1, 0) = 0.5;
  sph(1, 1) = 0.01;

  sph(2, 0) = 0.0;
  sph(2, 1) = 0.0;

  sph(3, 0) = 0.0;
  sph(3, 1) = 0.0;
}

void ic_three_particles(int n, SPH &sph) {

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

void ic_four_particles(int n, SPH &sph) {

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

void ic_dam_break(int n, SPH &sph) {

  int el = pow(n, 0.5);
  // Initial distance between the particles in both directions
  double step = 0.19 / (el - 1);
  // Starting position in x
  double posx = 0.01;
  double posy;
  // Assing the values in x for all particles
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      sph(0, i * el + j) = posx + double(rand()) / RAND_MAX / 100000;
      sph(2, i * el + j) = 0.0;
    }
    posx += step;
  }

  // For uniform distribution the step in y has to be equal to the step in x
  step = 0.19 / (el - 1);

  // Assing values in y for all particles
  for (int i = 0; i < el; i++) {
    posy = 0.01;
    for (int j = 0; j < el; j++) {
      sph(1, i * el + j) = posy + double(rand()) / RAND_MAX / 100000;
      sph(3, i * el + j) = 0.0;
      posy += step;
    }
  }
}

void ic_block_drop(int n, int n1, int n2, SPH &sph) {

  // Distance between neighbouring particles in x and y
  // 0.2 is the total distance in x and 0.3 in y
  double dx = 0.2 / double((n1 - 1));
  double dy = 0.3 / double((n2 - 1));

  // Starting position in x
  double posx = 0.1;
  double posy;
  int kx, ky;

  // Assing the values in x for all particles
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      kx = i * n2 + j;
      sph(0, kx) = posx + double(rand()) / RAND_MAX / 100000;
      sph(2, kx) = 0.0;
    }
    posx += dx;
  }

  // Assing the values in y for all particles
  for (int i = 0; i < n1; i++) {
    posy = 0.3;
    for (int j = 0; j < n2; j++) {
      ky = i * n2 + j;
      sph(1, ky) = posy + double(rand()) / RAND_MAX / 100000;
      sph(3, ky) = 0.0;
      posy += dy;
    }
  }
}

// Droplet
void ic_droplet(int n, SPH &sph) {

  double *xi = new double[n];
  double *yi = new double[n];
  int el = pow(n, 0.5);
  int kx;

  // For uniform distribution the step in y has to be equal to the step in x
  double step = 0.2 / (el - 1);
  double posx = 0.4; // Starting position in x
  double posy;       // Starting position in y

  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      xi[i * el + j] = posx;
    }
    posx += step;
  }

  step = 0.2 / (el - 1);

  for (int i = 0; i < el; i++) {
    posy = 0.6;
    for (int j = 0; j < el; j++) {
      yi[i * el + j] = posy;
      posy += step;
    }
  }
  kx = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (sqrt(pow((yi[i * el + j] - 0.7), 2) +
               pow((xi[i * el + j] - 0.5), 2)) <= 0.1) {
        sph(0, kx) = xi[i * el + j] + double(rand()) / RAND_MAX / 100000;
        sph(1, kx) = yi[i * el + j] + double(rand()) / RAND_MAX / 100000;
        sph(2, kx) = 0;
        sph(3, kx) = 0;
        kx++;
      }
    }
  }
  delete[] xi;
  delete[] yi;
}

// Defines the number of particles that will be in the circular region
int dropletn(int n) {

  // Process similar to dam break. Creates an initial square
  double *xi = new double[n];
  double *yi = new double[n];
  int el = pow(n, 0.5);
  double step = 0.2 / (el - 1);
  double posx = 0.4;
  double posy;

  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      xi[i * el + j] = posx;
    }
    posx += step;
  }

  step = 0.2 / (el - 1);

  for (int i = 0; i < el; i++) {
    posy = 0.6;
    for (int j = 0; j < el; j++) {
      yi[i * el + j] = posy;
      posy += step;
    }
  }

  // After the initial square is created, the number of particles that are in
  // that square and from a distance from the centre less or equal to the radius
  // of the circle is calculated
  int count = 0;
  for (int i = 0; i < el; i++) {
    for (int j = 0; j < el; j++) {
      if (sqrt(pow((yi[i * el + j] - 0.7), 2) +
               pow((xi[i * el + j] - 0.5), 2)) <= 0.1) {
        count++;
      } else {
        count += 0;
      }
    }
  }
  delete[] xi;
  delete[] yi;

  return count;
}