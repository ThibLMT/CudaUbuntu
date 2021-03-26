//
// Created by ThibLMT on 24/03/2021.
//

#include "sub_domain.cuh"
#include <stdio.h>
#include <math.h>

__global__ void initialize_backgrid(unsigned int *backgrid,int *backgrid_insert,geom_struct *geom)
{
    int size_backgrid,size_backgrid_insert;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    size_backgrid_insert=geom->sizex*geom->sizey*geom->sizez;
    size_backgrid=size_backgrid_insert*geom->sizel;

    for (int i = index; i < size_backgrid_insert; i+= stride)
    {
        backgrid_insert[i] = 0;
    }

    for (int i = index; i < size_backgrid; i+= stride)
    {
        backgrid[i] = 0;
    }
}