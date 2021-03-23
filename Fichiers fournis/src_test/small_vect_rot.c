/**
*\file tangential_force_mindlin.c
*\brief body of the function tangential_force_mindlin
*/
#include "def_types.h"
//#include "def_global_variables.h"
#include "small_vect_rot.h"
//extern double deltat; //!< Time step

vector small_vect_rot(vector u,vector nnew,vector nold,vector wi,vector wj,double deltat)
	{
	vector Wij,fu;
	double comp12,comp31,comp23; //!< Scalar components of the first cross product
	double wijmean; //!<
	//double dot_prod;

	// Firstly make a rotation due to the moving of the local base (nold -> nnews)
	// This "rolling part" is orthogonal to the plane containing nold and nnews
	
	comp12=nold.x*nnew.y-nold.y*nnew.x;
	comp31=nold.z*nnew.x-nold.x*nnew.z;
	comp23=nold.y*nnew.z-nold.z*nold.y;

	fu.x=u.x-u.y*comp12+u.z*comp31;
	fu.y=u.y+u.x*comp12-u.z*comp23;
	fu.z=u.z-u.x*comp31+u.y*comp23;

	// Then, Pivoting around the new normal direction. The rotation is equal to the average rotation of the two particles in contact
	
	wijmean=0.5*((wi.x+wj.x)*nnew.x+(wi.y+wj.y)*nnew.y+(wi.z+wj.z)*nnew.z);
	Wij.x=wijmean*nnew.x;
	Wij.y=wijmean*nnew.y;
	Wij.z=wijmean*nnew.z;

	fu.x=fu.x-fu.y*Wij.z*deltat+fu.z*Wij.y*deltat;
	fu.y=fu.y+fu.x*Wij.z*deltat-fu.z*Wij.x*deltat;
	fu.z=fu.z-fu.x*Wij.y*deltat+fu.y*Wij.x*deltat;

	// Just test, specified in the Fazekas thesis Projection on the rotational plane to ensure that vector u to be perpendicular to nnews
	/*
	dot_prod=fu.x*nnew.x+fu.y*nnew.y+fu.z*nnew.z;

	fu.x=fu.x-nnew.x*dot_prod;
	fu.y=fu.y-nnew.y*dot_prod;
	fu.z=fu.z-nnew.z*dot_prod;*/

	return(fu);
	}

