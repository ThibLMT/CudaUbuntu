/**
*\file write_vtk.h
*\brief header file of write_vtk function
*\fn 
*\brief this function writes vtk file to visualize the microstrure on paraview software
*\param *Nvtkfile
*\return void
*
*/

#include <stdio.h>
#include <stdlib.h>
 
#ifndef __WRITE_VTK_H
#define __WRITE_VTK_H

void write_vtk(const char *Nvtkfile);

#endif

