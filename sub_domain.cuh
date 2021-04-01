//
// Created by ThibLMT on 24/03/2021.
//
#include <stdio.h>
#include <stdlib.h>
#include "def_types.h"

#ifndef CUDAUBUNTU_SUB_DOMAIN_CUH
#define CUDAUBUNTU_SUB_DOMAIN_CUH

__global__ void initialize_backgrid(unsigned int *backgrid,int *backgrid_insert,geom_struct *geom);
__host__ void set_id_backgrid(int xv,int yv,int zv,unsigned int idp,unsigned int *backgrid, int *backgrid_insert,geom_struct *geom);
__global__ void insert_sph_backgrid(discrete_elt *particle, unsigned int *backgrid,int *backgrid_insert,geom_struct *geom);
__device__ unsigned int get_id_backgrid(int xv,int yv,int zv,int lv,unsigned int *backgrid,geom_struct *geom);
__device__ int detect_contact_sph_backgrid(discrete_elt *particle, double rad, geom_struct *geom,unsigned int idpart,unsigned int* list_part,unsigned int *backgrid);
__device__ int list_part_pot_sorting(unsigned int* array, int nelement);

#endif //CUDAUBUNTU_SUB_DOMAIN_CUH
