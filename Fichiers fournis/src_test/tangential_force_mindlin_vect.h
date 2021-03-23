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


#ifndef __TANGENTIAL_FORCE_MINDLIN_VECT_H
#define __TANGENTIAL_FORCE_MINDLIN_VECT_H

vector tangential_force_mindlin_vect(double deltanl,double fnijl,double mul, double Rij, double Gijl,vector ftan,vector vtanrel,double deltat);

#endif

