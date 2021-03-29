//
// Created by ThibLMT on 24/03/2021.
//
#include <stdio.h>
#include <stdlib.h>
#include "def_types.h"

#ifndef CUDAUBUNTU_SUB_DOMAIN_CUH
#define CUDAUBUNTU_SUB_DOMAIN_CUH

__global__ void initialize_backgrid(unsigned int *backgrid,int *backgrid_insert,geom_struct *geom);
void set_id_backgrid(int xv,int yv,int zv,unsigned int idp,unsigned int *backgrid, int *backgrid_insert,geom_struct *geom);
__global__ void insert_sph_backgrid(discrete_elt *particle, unsigned int *backgrid,int *backgrid_insert,geom_struct *geom);

#endif //CUDAUBUNTU_SUB_DOMAIN_CUH
