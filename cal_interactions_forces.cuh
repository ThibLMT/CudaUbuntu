//
// Created by ThibLMT on 01/04/2021.
//
#include <stdio.h>
#include <stdlib.h>

#ifndef CUDAUBUNTU_CAL_INTERACTIONS_FORCES_CUH
#define CUDAUBUNTU_CAL_INTERACTIONS_FORCES_CUH
#include "def_types.h"
__global__ void cal_interaction_forces(int idparti,discrete_elt *particle,geom_struct geom,unsigned int *backgrid);
#endif //CUDAUBUNTU_CAL_INTERACTIONS_FORCES_CUH
