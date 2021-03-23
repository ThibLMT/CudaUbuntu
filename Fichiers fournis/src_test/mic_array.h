/**
*\file mic array.h
*\brief header file of write_cyl_vtk function
*\fn 
*\brief this function writes vtk file to visualize the cylindre on paraview software
*\param *Nvtkfile
*\return void
*
*/

#include <stdio.h>
#include <stdlib.h>
 
#ifndef __MIC_ARRAY_H
#define __MIC_ARRAY_H
#include "def_types.h"

int insert_sphere_backgrid(double xpart, double ypart, double zpart, double radpart, unsigned int idpart,unsigned int *backgrid,int *backgrid_insert);
void alloc_backgrid(); 
unsigned int get_idp_backgrid(int xv,int yv,int zv,int lv,unsigned int *backgrid);
void set_idp_backgrid(int xv,int yv,int zv,unsigned int idp,unsigned int *backgrid, int *backgrid_insert);
int detect_contact_sphere_backgrid(double xin, double yin, double zin, double rad, unsigned int idpart,unsigned int* list_part,unsigned int *backgrid);
double dist2_min_point_periodic_voxel(double xpoint,double ypoint,double zpoint,int xv,int yv,int zv);
int list_part_pot_sorting_back(unsigned int* array, int nelement);
int detect_contact_sphere_backgrid_ini(double xin, double yin, double zin, double rad, unsigned int idpart,unsigned int* list_part,unsigned int *backgrid);

#endif

