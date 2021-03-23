/**
*\file dlvo_forces.h
*\brief header file of van_der_waals_force function
*\fn van_der_walls_force(double hij,double ri,double rj,double rij,int model_vdw)
*\brief This function computes the van der Waals function. Two models are implemented.
*\param hij distance beteween two centers of the particles i and j
*\param ri radius of particle i
*\param rj radius of particle j
*\param rij Equivalent radius
*\return double normal of the force of van der Walls
*\fn NormalElectrostaticForce(double hij,double ri, double rj,double rij)
*\brief This function computes the electrostatic force
*/

#include <stdio.h>
#include <math.h>
#include "def_const.h"

#ifndef __DLVO_FORCES_H
#define __DLVO_FORCES_H

double van_der_waals_force(double hij,double ri,double rj,double rij,int model_vdw);
double NormalElectrostaticForce(double hij,double ri, double rj,double rij);
#endif

