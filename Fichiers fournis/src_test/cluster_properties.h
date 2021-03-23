/**
*\file cluster_properties.h
*\brief header file of cluster_properties.c
*\fn void cluster_properties(void)
*\brief Subroutine to characterize clusters
*\param void
*\return void
*
*/
#include <stdio.h>
#include <stdlib.h>
 
#ifndef __CLUSTER_PROPERTIES_H
#define __CLUSTER_PROPERTIES_H

double *cluster_properties(double *freturn);

#endif

struct agglomerat      //!< data structure used for the description of clusters of agglomerated particles
	{
	int first;								//!< id of the first particle composing the cluster
	int np;									//!< number of particles composing the cluster
	double Rg;								//!< giration radius
	double Xg;								//!< X coordinate of the cluster's barycentre
	double Yg;								//!< Y coordinate of the cluster's barycentre
	double Zg;								//!< Z coordinate of the cluster's barycentre
	double Df;								//!< Z fractal dimension
	double Ks;								//!< structure factor
	double Comp;							//!< packing density of the cluster
	};
