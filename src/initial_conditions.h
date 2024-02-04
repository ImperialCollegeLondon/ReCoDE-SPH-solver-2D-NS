// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include <memory>

#include "fluid.h"

void icBasic(std::unique_ptr<Fluid> &fluidPtr, int n,
             std::vector<double> &positionX, std::vector<double> &positionY);

void icBlockDrop(std::unique_ptr<Fluid> &fluidPtr, int &n, double length,
                 double width, double centerX, double centerY);

void icDroplet(std::unique_ptr<Fluid> &fluidPtr, int &n, double radius,
               double centerX, double centerY);

int closestIntegerSqrt(int n);

int rectangleN(int n, double length, double width, int &n1, int &n2);

#endif