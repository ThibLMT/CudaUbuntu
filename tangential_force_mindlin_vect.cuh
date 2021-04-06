//
// Created by bluet on 06/04/2021.
//
#include <stdio.h>
#include <math.h>

#ifndef CUDAUBUNTU_TANGENTIAL_FORCE_MINDLIN_VECT_CUH
#define CUDAUBUNTU_TANGENTIAL_FORCE_MINDLIN_VECT_CUH

__device__ vector tangential_force_mindlin_vect(double deltanl,double fnijl,double mul, double Rij, double Gijl,vector ftan,vector vtanrel,double deltat);

#endif //CUDAUBUNTU_TANGENTIAL_FORCE_MINDLIN_VECT_CUH
