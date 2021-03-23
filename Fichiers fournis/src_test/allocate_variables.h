/**
*\file allocate_variables.h
*\brief header file of allocate_variables function
*\fn void allocate_variables(void)
*\brief function to allocate mic array for the moment
*\param void	
*\return void
*
*/
#include <stdio.h>
#include <stdlib.h>
#include "def_types.h"
 
#ifndef __ALLOCATE_VARIABLES_H
#define __ALLOCATE_VARIABLES_H

discrete_elt* allocate_particle(int nb_elements,size_t size);
void initialize_particle(discrete_elt *particle,geom_struct geom);
void give_properties_particle(discrete_elt *particle,double unity,material_data properties);
void set_vect_0(vector *vec);
void update_particle(discrete_elt *particle,geom_struct geom);

#endif
