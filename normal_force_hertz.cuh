//
// Created by bluet on 06/04/2021.
//

#include <stdio.h>
#include <math.h>

#ifndef CUDAUBUNTU_NORMAL_FORCE_HERTZ_CUH
#define CUDAUBUNTU_NORMAL_FORCE_HERTZ_CUH

__device__ double normal_force_hertz(double deltan,double deltandot, double rij, double Eij, double Aij);

#endif //CUDAUBUNTU_NORMAL_FORCE_HERTZ_CUH
