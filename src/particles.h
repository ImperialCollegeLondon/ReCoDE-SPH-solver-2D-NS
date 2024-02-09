#ifndef PARTICLES_H
#define PARTICLES_H

#include <functional>
#include <vector>

class particles {
 public:
  /******** CONSTRUCTORS/DESTRUCTOR********/

  explicit particles(
      const unsigned
          nNew);  // User defined constructor for allocating
                  // the dimensions of the arrays. One parameter
                  // constructors should declared `explicit` to avoid
                  // implicit conversions. Check
                  // https://google.github.io/styleguide/cppguide.html#Implicit_Conversions

  // Getter functions

  // Function to get the number of particles
  inline unsigned int getNumberOfParticles() { return nbParticles; }

  // Function to get the position x
  inline double getPositionX(int k) { return positionX[k]; }

  // Function to get the position y
  inline double getPositionY(int k) { return positionY[k]; }

  // Function to get the velocity x
  inline double getVelocityX(int k) { return velocityX[k]; }

  // Function to get the velocity y
  inline double getVelocityY(int k) { return velocityY[k]; }

  // Setter functions

  // Function to set the position x
  inline void setPositionX(int k, double newPositionX) {
    positionX[k] = newPositionX;
  }

  // Function to set the position y
  inline void setPositionY(int k, double newPositionY) {
    positionY[k] = newPositionY;
  }

  // Function to set the velocity x
  inline void setVelocityX(int k, double newVelocityX) {
    velocityX[k] = newVelocityX;
  }

  // Function to set the velocity y
  inline void setVelocityY(int k, double newVelocityY) {
    velocityY[k] = newVelocityY;
  }

  // Calculation functions

  // Function to calculate the matrix with rij
  void calculateParticleDistance();

 protected:
  unsigned int nbParticles;  // number of particles and characteristic size of
                             // the class arrays

  // Positions
  std::vector<double> positionX;
  std::vector<double> positionY;

  // Velocities
  std::vector<double> velocityX;
  std::vector<double> velocityY;

  std::vector<double> particleSpeedSq;  // u(i)^2+v(i)^2
};

#endif