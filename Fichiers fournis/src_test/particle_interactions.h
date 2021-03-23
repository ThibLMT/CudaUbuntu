/**
*\fn void choc(unsigned int ni,unsigend int nj,int k_ij,double *fcontloc)
*\brief Routine to start the computation of contacts, forces and displacements
*\param ni id of the particle
*\param nj id of the particle in contact with ni
*\param k_ij	number of the contact conerning the particle i. This number is the position of the contact i and j in the history arrays (linked to the structure)
*\param fcontloc Contact force in the local basis
*\return void
*/


#ifndef __PARTICLE_INTERACTIONS_H
#define __PARTICLE_INTERACTIONS_H

#include "def_types.h"


void particle_interactions(unsigned int ni,unsigned int nj,int k_ij,vector* forceji_pointer,vector* torqueji_pointer,discrete_elt *particle,geom_struct geom);

#endif
