//
// Created by bluet on 06/04/2021.
//
#include <stdio.h>
#include <math.h>
#include "def_const.h"

#ifndef CUDAUBUNTU_DLVO_FORCES_CUH
#define CUDAUBUNTU_DLVO_FORCES_CUH

__device__ double van_der_waals_force(double hij,double ri,double rj,double rij,int model_vdw);
__device__ double NormalElectrostaticForce(double hij,double ri, double rj,double rij);

#endif //CUDAUBUNTU_DLVO_FORCES_CUH
