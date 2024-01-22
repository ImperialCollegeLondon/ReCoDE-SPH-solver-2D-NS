// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include "sph.h"

SPH ic_basic(unsigned int n, double *position_x, double *position_y);

SPH ic_block_drop(unsigned int &n, double length, double width, double center_x,
                  double center_y);

SPH ic_droplet(unsigned int &n, double radius, double center_x,
               double center_y);

unsigned int closest_integer_sqrt(unsigned int n);

unsigned int rectangle_n(unsigned int n, double length, double width,
                         unsigned int &n1, unsigned int &n2);

#endif