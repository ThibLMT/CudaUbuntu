/**
*\file read_micro.h
*\brief contain function to read data in microfile
*/
#include <stdio.h>
#include <stdlib.h>
#include "def_types.h"
#ifndef CUDAUBUNTU_READ_MICRO_H
#define CUDAUBUNTU_READ_MICRO_H

void microfile_read_particle(const char *Nfile,discrete_elt *particle,geom_struct *geom);
unsigned int microfile_read_npart(const char *Nfile);
void microcontfile_read_contact(const char *Nfile,discrete_elt *particle,geom_struct *geom);

#endif //CUDAUBUNTU_READ_MICRO_H
