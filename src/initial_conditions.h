// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include <memory>

#include "fluid.h"

void icBasic(std::unique_ptr<Fluid> &fluidPtr, unsigned int n,
             std::vector<double> &positionX, std::vector<double> &positionY);

void icBlockDrop(std::unique_ptr<Fluid> &fluidPtr, unsigned int n,
                 double length, double width, double centerX, double centerY);

void icDroplet(std::unique_ptr<Fluid> &fluidPtr, unsigned int n, double radius,
               double centerX, double centerY);

unsigned int closestIntegerSqrt(unsigned int n);

unsigned int rectangleN(unsigned int n, double length, double width, int &n1,
                        int &n2);

#endif