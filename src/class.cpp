
#include "class.h"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
using namespace std;

// Defining the user defined constructor
SPH::SPH(const unsigned n_new) {

  n = n_new;

  for (int i = 0; i < 4; i++) {

    vMatrix[i] = new double[n]; // Initializing an array that will be used to
                                // pass and stor the values of the initial
                                // positions and velocities inside the class
  }

  x = new double[n];  // Array to store the positions in x of the particles
  y = new double[n];  // Array to store the positions in y of the particles
  vx = new double[n]; // Array to store the values of the velocity in x of the
                      // particles
  vy = new double[n]; // Array to store the values of the velocity in y of the
                      // particles
  r = new double[n * n]; // Array to store the values of the distances between
                         // the particles
  q = new double[n * n]; // Array to store the values of q
  ro = new double[n];    // Array to store the densitites of the particles
  p = new double[n];     // Array to store the pressure on each particle
  vi = new double[n];    // Array to store the norm of the velocity, used in the
                      // calculation of the kinetic energy

  // Array that will be broadcasted to the all the processors and will contain
  // the above values
  xyv = new double[4 * n];
}

// Destructor
SPH::~SPH() {
  delete[] vMatrix;
  delete[] ro;
  delete[] p;
  delete[] r;
  delete[] vi;
  delete[] x;
  delete[] y;
  delete[] vx;
  delete[] vy;
  delete[] Freduce;
  delete[] Ftake;
  delete[] xyv;
}

// Defining the Overloading of ()
double &SPH::operator()(unsigned row, unsigned col) {

  return vMatrix[row][col];
}

// Defining the Overloading of > to pass the number of the current time
// iteratiion to the class
int &SPH::operator>(unsigned t_new) {

  t = t_new;
  return t;
}

// Defining the Overloading of >> to pass the time-step inside tha class
double &SPH::operator>>(double dt_new) {

  dt = dt_new;
  return dt;
}

// Defining the Overloading of < to pass the radius of indfluence inside the
// class
double &SPH::operator<(double h_new) {

  h = h_new;
  return h;
}

/**Assigning values to x: The values of the positions on the
 * x-axis are stored in the first row of the vMatrix, and they are now assigned
 * to a single vector **/
void SPH::x0() {

  for (int i = 0; i < n; i++) {

    x[i] = vMatrix[0][i];
  }
}

/**Assigning values to y: The values of the positions on the
 * y-axis are stored int he second row of the vMatrix, and they are now assigned
 * to a single vector **/
void SPH::y0() {

  for (int i = 0; i < n; i++) {

    y[i] = vMatrix[1][i];
  }
}

/**Assigning values to vx: The values of the initial velocities in the x
 * direction are stored int the third row of the vMatrix, and they are now
 * assigned to a single vector **/
void SPH::vx0() {

  for (int i = 0; i < n; i++) {

    vx[i] = vMatrix[2][i];
  }
}

/**Assigning values to vy: The values of the initial velocities in the x
 * direction are stored int the fourth row of the vMatrix, and they are now
 * assigned to a single vector **/
void SPH::vy0() {

  for (int i = 0; i < n; i++) {

    vy[i] = vMatrix[3][i];
  }
}

// Creating a funcation to calculate the matrix with rij
void SPH::rVec() {

  double dx;
  double dy;

  for (int i = 0; i < n; i++) {

    for (int j = 0; j < n; j++) {

      dx = x[i] - x[j];
      dy = y[i] - y[j];

      r[i * n + j] = sqrt(dx * dx + dy * dy);
    }
  }
}

// creating the function for density

void SPH::den() {

  double phi;
  double pre = (4.0 / (M_PI * h * h)); // Precalculated value
  double hinv = 1.0 / h;               // This is to avoid many divisions

  // finding φ
  for (int i = 0; i < n; i++) {

    ro[i] = 0;

    for (int j = 0; j < n; j++) {

      q[i * n + j] = abs(r[i * n + j] * hinv);

      if (q[i * n + j] < 1) {

        phi = pre * (1.0 - q[i * n + j] * q[i * n + j]) *
              (1.0 - q[i * n + j] * q[i * n + j]) *
              (1.0 - q[i * n + j] * q[i * n + j]);

      }

      else {

        phi = 0.0;
      }

      ro[i] += m * phi;
    }
  }
}

// creating the function for the pressure

void SPH::pressure() {

  for (int i = 0; i < n; i++) {

    p[i] = k * (ro[i] - roo);
  }
}

double SPH::FiP(double *x_y) {

  double sum = 0.0;                          // Initializing the sumation
  double pre = (-30.0 / (M_PI * h * h * h)); // precalculated value

  for (int j = 0; j < n; j++) {

    if (i == j) {
    } else {

      if (q[i * n + j] < 1) {

        // Equation (5)
        sum += (m / ro[j]) * ((p[i] + p[j]) / 2.0) * (pre * (x_y[i] - x_y[j])) *
               (((1.0 - q[i * n + j]) * (1.0 - q[i * n + j])) / q[i * n + j]);
      }

      else {
      }
    }
  }

  // Equation (4)

  return -sum;
}

double SPH::FiV(double *v) {

  double phisq;

  double sum = 0.0;                             // Initializing the sumation
  double pre = (40.0 / (M_PI * h * h * h * h)); // precalculated value

  for (int j = 0; j < n; j++) {

    if (i == j) {
    }

    else {

      if (q[i * n + j] < 1) {

        // Equation (7)
        sum += (m / ro[j]) * (v[i] - v[j]) * (pre * (1.0 - q[i * n + j]));
      }
    }
  }

  // Equation (6)

  return -miou * sum;
}

// Creating the function for the gravity force
double SPH::FiG() { return -ro[i] * g; }

// Function for initialisation of the time integration scheme - velocity
double SPH::vInit(double *v, double &Fip, double &Fiv, double &Fig) {

  double a;

  a = (Fip + Fiv + Fig) / ro[i];

  return v[i] + a * dt * 0.5;
}

// Function for time integration - velocity
double SPH::v_x_y(double *v, double &Fip, double &Fiv, double &Fig) {

  double a;
  a = (Fip + Fiv + Fig) / ro[i];

  return v[i] + a * dt;
}

// Function to find the mass of the particles before the simulation starts
void SPH::mass() {

  rVec();
  den();
  double sumden = 0.0;
  for (int i = 0; i < n; i++) {

    sumden += ro[i];
  }

  m = n * roo / sumden;
}

// Function to perform the spatial iterations
void SPH::spatial() {

  for (i = 0; i < n; i++) {

    pressure();

    // Gathering the forces calculated by the processors

    Fipx = FiP(x);

    Fivx = FiV(vx);

    Fipy = FiP(y);

    Fivy = FiV(vy);

    Figy = FiG();

    // First step to initialise the scheme
    if (t == 0) {

      vx[i] = vInit(vx, Fipx, Fivx, Figx);
      x[i] = x[i] + vx[i] * dt; // inlined
      vy[i] = vInit(vy, Fipy, Fivy, Figy);
      y[i] = y[i] + vy[i] * dt; // inlined

    }

    // Leap frog scheme
    else {

      vx[i] = v_x_y(vx, Fipx, Fivx, Figx);
      x[i] = x[i] + vx[i] * dt; // inlined
      vy[i] = v_x_y(vy, Fipy, Fivy, Figy);
      y[i] = y[i] + vy[i] * dt; // inlined
    }

    // Boundary Conditions
    if (x[i] < h) {

      x[i] = h;
      vx[i] = -e * vx[i];
    }

    if (x[i] > 1.0 - h) {

      x[i] = 1.0 - h;
      vx[i] = -e * vx[i];
    }

    if (y[i] < h) {

      y[i] = h;
      vy[i] = -e * vy[i];
    }

    if (y[i] > 1.0 - h) {

      y[i] = 1.0 - h;
      vy[i] = -e * vy[i];
    }
  }
}

// Function to return the position x
double SPH::retx(int l) { return x[l]; }

// Function to return the position y
double SPH::rety(int l) { return y[l]; }

// Function to calculate the kinetic energy
double SPH::Ek() {

  double sum = 0;
  for (int i = 0; i < n; i++) {

    vi[i] = vx[i] * vx[i] + vy[i] * vy[i];

    sum += vi[i];
  }

  return 0.5 * m * sum;
}

// Function to calculate the potential energy
double SPH::Ep() {

  double sum = 0;
  for (int i = 0; i < n; i++) {

    sum += y[i];
  }

  return m * g * sum;
}

// Getting the new data calculated at each time step
void SPH::getdata() {

  for (int k = 0; k < n; k++) {

    x[k] = xyv[k];
    y[k] = xyv[n + k];
    vx[k] = xyv[2 * n + k];
    vy[k] = xyv[3 * n + k];
  }
}
