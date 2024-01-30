// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include "fluid.h"

void ic_basic(fluid **fluid_ptr, int n, double *position_x, double *position_y);

void ic_block_drop(fluid **fluid_ptr, int &n, double length, double width,
                   double center_x, double center_y);

void ic_droplet(fluid **fluid_ptr, int &n, double radius, double center_x,
                double center_y);

int closest_integer_sqrt(int n);

int rectangle_n(int n, double length, double width, int &n1, int &n2);

#endif