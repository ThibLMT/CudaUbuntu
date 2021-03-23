/**
*\fn double compacity_z(double zmin, double zmax,double sxy)
*\brief This function computes the compacity in each layer defined by zmin & zmax
*\author K. Kimbonguila
*\param lmin represents the lower altitude of the layer
*\param lmax represents the upper altitude of the layer
*\param sxy surface of the volume's base
*\return compacity
*\fn void compacity_profil_z(const char *Nprofilfile)
*\brief This function computes the gradient of compacity folloving the z axis
*\author K. El-Cheikh
**/


#include <stdio.h>
#include <stdlib.h>
#include "def_const.h"
#ifndef __CAL_COMPACITY_H
#define __CAL_COMPACITY_H

double compacity_z(double zmin, double zmax,double sxy);
void compacity_profil(long int npoints,int ifile);
void compacity_profil_z(const char *Nprofilfile,double sxy);



#endif
