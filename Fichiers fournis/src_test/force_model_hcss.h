/**
*\fn normal_force_hcss (...)
*\fn tangential_force_hcss (...)
*\brief function used to compute the rolling friction torque
*\param void
*\return vector
*/




#ifndef __FORCE_MODEL_HCSS_H
#define __FORCE_MODEL_HCSS_H

#include "def_types.h"
double normal_force_hcss(double deltan);
vector tangential_force_hcss(double deltan,vector nij,double r1,double r1eq,vector v1,vector w1,double r2,double r2eq,vector v2,vector w2);

#endif
