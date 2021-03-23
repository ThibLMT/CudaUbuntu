/**
*\file fluid_forces.h
*\brief header file of fluid_force and fluid_torque function
*\fn fluid_force(vector vrif)
*\brief 
*\return vector 
*\fn fluid_torque(vector wrif);
*\brief 
*/

#include <stdio.h>
#include <math.h>
#include "def_types.h"

#ifndef __FLUID_FORCES_H
#define __FLUID_FORCES_H
vector fluid_force(double radius,vector vpart,vector vflu);
vector fluid_torque(double radius,vector wpart,vector wflu);
#endif

