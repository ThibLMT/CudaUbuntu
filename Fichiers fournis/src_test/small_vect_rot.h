/**
*\file tangential_force_mindlin.h
*\brief header file of tangential_force_mindlin function
*\fn double tangential_force_mindlin(double deltan,double deltatan,double fnij,double mu, double Eij, double Gij)
*\brief this function computes the tangential force between two grains. The norme of tangential force is limited by Coulomb criteria
*\param double deltan,double deltatan,double fnij,double mu, double Eij, double Gij
*\return double
*
*/
#include <stdio.h>
#include <math.h>


#ifndef __SMALL_VECT_ROT_H
#define __SMALL_VECT_ROT_H

vector small_vect_rot(vector u,vector nnew,vector nold,vector wi,vector wj,double deltat);

#endif

