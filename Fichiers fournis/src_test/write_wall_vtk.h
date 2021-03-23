/**
*\file write_wall_vtk.h
*\brief header file of write_wall_vtk function
*\fn 
*\brief this function writes vtk file to visualize the rigid walls on paraview software
*\param *Nvtkfile
*\return void
*
*/

#include <stdio.h>
#include <stdlib.h>
 
#ifndef __WRITE_WALL_VTK_H
#define __WRITE_WALL_VTK_H

void write_wall_vtk(const char *Nvtkfile);

#endif

