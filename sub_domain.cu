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

void set_id_backgrid(int xv,int yv,int zv,unsigned int idp,unsigned int *backgrid, int *backgrid_insert,geom_struct *geom)
{
    int index_backgrid,index_backgrid_insert,lv;
    // Keep the last position lv
    index_backgrid_insert=zv+yv*geom->sizez+xv*geom->sizez*geom->sizey;
    lv=backgrid_insert[index_backgrid_insert];
    index_backgrid=lv+zv*geom->sizel+yv*geom->sizel*geom->sizez+xv*geom->sizel*geom->sizez*geom->sizey;
    backgrid[index_backgrid]=idp;
    backgrid_insert[index_backgrid_insert]++;
}
