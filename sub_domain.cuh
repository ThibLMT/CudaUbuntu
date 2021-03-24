//
// Created by ThibLMT on 24/03/2021.
//
#include <stdio.h>
#include <stdlib.h>
#include "def_types.h"

#ifndef CUDAUBUNTU_SUB_DOMAIN_CUH
#define CUDAUBUNTU_SUB_DOMAIN_CUH

unsigned int * allocation_backgrid(geom_struct *geom);
int * allocation_backgrid_insert(geom_struct *geom);

#endif //CUDAUBUNTU_SUB_DOMAIN_CUH
