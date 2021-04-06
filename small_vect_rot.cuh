//
// Created by bluet on 06/04/2021.
//
#include <stdio.h>
#include <math.h>

#ifndef CUDAUBUNTU_SMALL_VECT_ROT_CUH
#define CUDAUBUNTU_SMALL_VECT_ROT_CUH

__device__ vector small_vect_rot(vector u,vector nnew,vector nold,vector wi,vector wj,double deltat);

#endif //CUDAUBUNTU_SMALL_VECT_ROT_CUH

