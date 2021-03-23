/**
*\file add_obj_spheres.h
*\brief header file of read_micro function
*\fn int add_obj_spheres(const char *Nobjfile, int iinit,Flag check_overlap)
*\brief function to read a file describing the position and radius of objects. These objects are spheres witch can overlap together.
*\param *Nobjfile name of the object file
*\param iinit id if the object (equal to iobj + npart)
*\param check_overlap
*\return int numver of the objects defined in the object file
*
*/
#include <stdio.h>
#include <stdlib.h>
 
#ifndef __ADD_OBJ_SPHERES_H
#define __ADD_OBJ_SPHERES_H

int add_obj_spheres(const char *Nobjfile, unsigned int iinit,Flag check_overlap);

#endif

