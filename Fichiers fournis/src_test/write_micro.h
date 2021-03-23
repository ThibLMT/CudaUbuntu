/**
*\file write_micro.h
*\brief header file of write_micro function
*\fn void write_micro(char )
*\brief function to write typical micro and contact files
*\param void	
*\return void
*
*/
#include <stdio.h>
#include <stdlib.h>
 
#ifndef __WRITE_MICRO_H
#define __WRITE_MICRO_H
#include "def_types.h"
void microfile_write(const char *Nfile,discrete_elt *particle,geom_struct geom);

#endif

