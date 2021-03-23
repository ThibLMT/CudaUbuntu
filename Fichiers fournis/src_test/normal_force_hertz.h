/**
*\file normal_force_hertz.h
*\brief header file of normal_force_hertz function
*\fn double normal_force_hertz(double deltan,double deltandot, double rij, double Eij, double Aij)
*\brief This function computes the elastic (Hertz model) and dissipative components of the normal contact forces
*\param deltan indentation
*\param deltandot Normal
*\param rij
*\param Eij
*\param Aij
*\return double
*/

#include <stdio.h>
#include <math.h>


#ifndef __NORMAL_FORCE_HERTZ_H
#define __NORMAL_FORCE_HERTZ_H

double normal_force_hertz(double deltan,double deltandot, double rij, double Eij, double Aij);

#endif

