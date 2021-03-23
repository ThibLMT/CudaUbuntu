/**
*\file add_cone.h
*\brief header file of add_cone function
*\fn struct sphere add_cone(int nwall,double height,double rbot,double rtop)
*\brief function to add a conic boudary condition
*\param struct sphere wall_i,int dist,vector vectnorm	
*\return void
*
*/



#ifndef __ADD_CONE_H
#define __ADD_CONE_H

#include "def_types.h"

void add_cone(int nwall);

#endif

