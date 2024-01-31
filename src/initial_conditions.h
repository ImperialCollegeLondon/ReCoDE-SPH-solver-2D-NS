// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include "fluid.h"

void icBasic(Fluid *&fluidPtr, int n, double *positionX, double *positionY);

void icBlockDrop(Fluid *&fluidPtr, int &n, double length, double width,
                 double centerX, double centerY);

void icDroplet(Fluid *&fluidPtr, int &n, double radius, double centerX,
               double centerY);

int closestIntegerSqrt(int n);

int rectangleN(int n, double length, double width, int &n1, int &n2);

#endif