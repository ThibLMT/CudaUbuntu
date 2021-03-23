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
#include "def_types.h"
 
#ifndef __SUB_DOMAIN_H
#define __SUB_DOMAIN_H

unsigned int * allocation_backgrid(geom_struct geom);
int * allocation_backgrid_insert(geom_struct geom);
void initialize_backgrid(unsigned int *backgrid,int *backgrid_insert,geom_struct geom);
unsigned int get_id_backgrid(int xv,int yv,int zv,int lv,unsigned int *backgrid,geom_struct geom);
void set_id_backgrid(int xv,int yv,int zv,unsigned int idp,unsigned int *backgrid, int *backgrid_insert,geom_struct geom);
int insert_sph_backgrid(double xpart, double ypart, double zpart, double radpart, unsigned int idpart,unsigned int *backgrid,int *backgrid_insert,geom_struct geom);
int detect_contact_sph_backgrid(discrete_elt *particle, double rad, geom_struct geom,unsigned int idpart,unsigned int* list_part,unsigned int *backgrid);
int list_part_pot_sorting(unsigned int* array, int nelement);

#endif

