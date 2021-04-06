//
// Created by ThibLMT on 01/04/2021.
//

#ifndef CUDAUBUNTU_PARTICLE_INTERACTIONS_CUH
#define CUDAUBUNTU_PARTICLE_INTERACTIONS_CUH

#include "def_types.h"


__device__ void particle_interactions(unsigned int ni,unsigned int nj,int k_ij,vector* forceji_pointer,vector* torqueji_pointer,discrete_elt *particle,geom_struct *geom);
__device__ vector tangential_force_mindlin_vect(double deltanl,double fnijl,double mul, double Rij, double Gijl,vector ftan,vector vtanrel,double deltat);

#endif //CUDAUBUNTU_PARTICLE_INTERACTIONS_CUH
