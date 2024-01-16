// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include "sph.h"

SPH ic_basic(int n, double *position_x, double *position_y);

SPH ic_block_drop(int n, double length, double width, double center_x,
                  double center_y);

SPH ic_droplet(int n, double radius, double center_x, double center_y);

int closest_power_of_two(int n);

int rectangle_n(int n, double length, double width, int &n1, int &n2);

#endif