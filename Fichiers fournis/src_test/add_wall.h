/**
*\file add_wall.h
*\brief header file of add_wall function
*\fn struct sphere add_wall(int nwall,vector Rp,vector vnorm)
*\brief function to read typical microtructure and contact files
*\param struct sphere wall_i,int dist,vector vectnorm	
*\return void
*
*/
#include <stdio.h>
#include <stdlib.h>


#ifndef __ADD_WALL_H
#define __ADD_WALL_H

struct sphere add_wall(int nwall,vector Rp,vector vnorm);

#endif

