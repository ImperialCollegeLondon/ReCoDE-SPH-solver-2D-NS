#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
class particles {
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

  // Distances
  std::vector<double>
      distance;  // Array to store the distances between the particles
  std::vector<double>
      distanceQ;  // Array to store the values of the normalised distance q

 public:
  /******** CONSTRUCTORS/DESTRUCTOR********/

  explicit particles(
      const unsigned
          nNew);  // User defined constructor for allocating
                  // the dimensions of the Matrix. One parameter
                  // constructors should declared `explicit` to avoid
                  // implicit conversions. Check
                  // https://google.github.io/styleguide/cppguide.html#Implicit_Conversions

  /**********OVERLOADINGS**********/

  double &operator()(unsigned row, unsigned col);

  // Getter functions

  // Function to get the number of particles
  int getNumberOfParticles();

  // Function to get the position x
  double getPositionX(int k);

  // Function to get the position y
  double getPositionY(int k);

  // Function to get the normalised distance between the particles
  double getDistanceQ(int k);

  // Calculation functions

  // Function to calculate the matrix with rij
  void calculateParticleDistance();
};

#endif