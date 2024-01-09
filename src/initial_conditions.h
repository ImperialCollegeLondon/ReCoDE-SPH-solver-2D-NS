// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include "particles.h"

void ic_one_particle(int n, particles &fluid);

void ic_two_particles(int n, particles &fluid);

void ic_three_particles(int n, particles &fluid);

void ic_four_particles(int n, particles &fluid);

void ic_dam_break(int n, particles &fluid);

void ic_block_drop(int n, int n1, int n2, particles &fluid);

void ic_droplet(int n, particles &fluid);

int dropletn(int n);

#endif