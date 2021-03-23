/**
*\fn vector rolling_resistance_torque(...)
*\brief function used to compute the rolling friction torque
*\param vector
*\return vector
*/




#ifndef __ROLLING_RESISTANCE_TORQUE_H
#define __ROLLING_RESISTANCE_TORQUE_H

#include "def_types.h"
vector rolling_resistance_torque(double fnijl,double murol, double Rij, vector resis_tor,vector wi,vector wj);

#endif

