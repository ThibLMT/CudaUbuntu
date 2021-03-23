/**
*\file add_wall.c
*\brief body of the function add_wall
*/
#include "def_global_variables.h"
#include "add_wall.h"

struct sphere add_wall(int nwall,vector Rp,vector vnorm)
	{
	int i,j,k,lc;
	int xp,yp,zp;  //!< Coordinate of the wall in pixel belonging x, y and z axes.
	struct sphere part_wall;
	int ltab;
	
	zp= (int) Rp.z ;
	xp= (int) Rp.x ;
	yp= (int) Rp.y ;
	
	
	if(Rp.x-xp>0.5){xp+=1;}
	else if(Rp.x-xp<-0.5){xp-=1;}
	
	if(Rp.y-yp>0.5){yp+=1;}
	else if(Rp.y-yp<-0.5){yp-=1;}
	
	if(Rp.z-zp>0.5){zp+=1;}
	else if(Rp.z-zp<-0.5){zp-=1;}
		
	//Wall is modeled by flat plane of the system considered as a sphere of infinite radius
	//Wall z =cte
	if (vnorm.z != 0 )
		{
		if ((vnorm.z > 0)&&(Rp.z - (double) zp) > 0.5 ) {zp++;}	
		
		for(i=0;i<syssizex;i++)
			{
			for(j=0;j<syssizey;j++)
				{
				lc=0;
				while (mic_boundary[i][j][zp][lc]!=0) {lc++;}
				mic_boundary[i][j][zp][lc]=nwall;
				mic_insert_boundary[i][j][zp]++;
				}
			}
		}
	//Wall x =cte	
	if (vnorm.x != 0 )
		{		
		if ((vnorm.x > 0)&&(Rp.x - (double) xp) > 0.5 ) {xp++;}						
		for(j=0;j<syssizey;j++)
			{
			for(k=0;k<syssizez;k++)
				{
				lc=0;
				while (mic_boundary[xp][j][k][lc]!=0) {lc++;}
				mic_boundary[xp][j][k][lc]=nwall;
				mic_insert_boundary[xp][j][k]++;
				}
			}
		}
		
	//Wall y =cte	
	if (vnorm.y != 0 )
		{
		if ((vnorm.y > 0)&&(Rp.y - (double) yp) > 0.5 ) {yp++;}					
		for(i=0;i<syssizex;i++)
			{
			for(k=0;k<syssizez;k++)
				{
				lc=0;
				while (mic_boundary[i][yp][k][lc]!=0) {lc++;}
				mic_boundary[i][yp][k][lc]=nwall;
				mic_insert_boundary[i][yp][k]++;
				}
			}
		}
	
	
	// Geometrie 
	part_wall.radius=1000000.0;
	part_wall.Ri.x=Rp.x + vnorm.x * part_wall.radius;
	part_wall.Ri.y=Rp.y + vnorm.y * part_wall.radius;
	part_wall.Ri.z=Rp.z + vnorm.z * part_wall.radius;
				
	// Initialisation vitesses 
	part_wall.Vi=vect0;
	part_wall.Wi=vect0;
				
	// Initialisation forces et moments
	part_wall.Fi=vect0;
	part_wall.Mi=vect0;

	// printing some information for logfile and screen
	fprintf(stdout,"wall : id %d - Position %d %d %d - Normal %1f %1f %1f \n",nwall,xp,yp,zp,vnorm.x,vnorm.y,vnorm.z);	
	fprintf(flogfile," Characteritics of the wall- id:%d \n",nwall);
	fprintf(flogfile," =======================================\n");

	fprintf(flogfile,"Wall position : %d %d %d\n",xp,yp,zp);
	fprintf(flogfile,"(center of the particle : %3.2lf; %3.2lf; %3.2lf)\n",part_wall.Ri.x,part_wall.Ri.y,part_wall.Ri.z );

	return part_wall;
	} 


