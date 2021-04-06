//
// Created by bluet on 06/04/2021.
//

#ifndef CUDAUBUNTU_ROLLING_RESISTANCE_TORQUE_CUH
#define CUDAUBUNTU_ROLLING_RESISTANCE_TORQUE_CUH

#include "def_types.h"
__device__ vector rolling_resistance_torque(double fnijl,double murol, double Rij, vector resis_tor,vector wi,vector wj);

#endif //CUDAUBUNTU_ROLLING_RESISTANCE_TORQUE_CUH
