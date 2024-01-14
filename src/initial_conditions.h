// Functions to create the initial conditions
#ifndef IC_H
#define IC_H
#include "sph_2D.h"

void ic_one_particle(int n, sph_2D &solver);

void ic_two_particles(int n, sph_2D &solver);

void ic_three_particles(int n, sph_2D &solver);

void ic_four_particles(int n, sph_2D &solver);

void ic_dam_break(int n, sph_2D &solver);

void ic_block_drop(int n, int n1, int n2, sph_2D &solver);

void ic_droplet(int n, sph_2D &solver);

int dropletn(int n);

#endif