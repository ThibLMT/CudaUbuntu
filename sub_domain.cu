//
// Created by ThibLMT on 24/03/2021.
//

#include "sub_domain.cuh"
#include <stdio.h>
#include <math.h>

unsigned int * allocation_backgrid(geom_struct *geom)
{
    int nb_elements;
    unsigned int *backgrid;
    nb_elements=geom->sizex*geom->sizey*geom->sizez*geom->sizel;
    cudaMallocManaged(&backgrid,nb_elements * sizeof(unsigned int));
    return backgrid;
}

int * allocation_backgrid_insert(geom_struct *geom)
{
    int nb_elements;
    int *backgrid_insert;
    nb_elements=geom->sizex*geom->sizey*geom->sizez;
    cudaMallocManaged(&backgrid_insert,nb_elements * sizeof(int));
    return backgrid_insert;
}