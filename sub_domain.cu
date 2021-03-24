//
// Created by ThibLMT on 24/03/2021.
//

#include "sub_domain.cuh"
#include <stdio.h>
#include <math.h>

unsigned int * allocation_backgrid(geom_struct *geom)
{
    int nb_elements;
    nb_elements=geom->sizex*geom->sizey*geom->sizez*geom->sizel;
    return static_cast<unsigned int *>(malloc(nb_elements * sizeof(unsigned int)));
}

int * allocation_backgrid_insert(geom_struct *geom)
{
    int nb_elements;
    nb_elements=geom->sizex*geom->sizey*geom->sizez;
    return static_cast<int *>(malloc(nb_elements * sizeof(int)));
}