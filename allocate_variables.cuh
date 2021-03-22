//
// Created by ThibLMT on 22/03/2021.
//

#include "def_types.h"

#ifndef CUDAUBUNTU_ALLOCATE_VARIABLES_CUH
#define CUDAUBUNTU_ALLOCATE_VARIABLES_CUH

__global__ void initialize_particle(discrete_elt *particle,geom_struct *geom);

#endif //CUDAUBUNTU_ALLOCATE_VARIABLES_CUH
