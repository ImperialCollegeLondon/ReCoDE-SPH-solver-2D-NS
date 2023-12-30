// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include "sph.h"

void ic_one_particle(int n, SPH &sph);

void ic_two_particles(int n, SPH &sph);

void ic_three_particles(int n, SPH &sph);

void ic_four_particles(int n, SPH &sph);

void ic_dam_break(int n, SPH &sph);

void ic_block_drop(int n, int n1, int n2, SPH &sph);

void ic_droplet(int n, SPH &sph);

int dropletn(int n);

#endif