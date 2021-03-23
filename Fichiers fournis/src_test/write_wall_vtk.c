
/**
*\file write_wall_vtk.c
*\brief body of the function write_wall_vtk
*/

#include "def_global_variables.h"
#include "write_wall_vtk.h"


void write_wall_vtk(const char *Nvtkfile)
	{	
	FILE *fparaview;
	int nwall;  //!< Number of rigid walls
	int i;

	fparaview = fopen(Nvtkfile,"w");
	if (fparaview != NULL) 
		{
		// Writing typical header for a vtk file
		fprintf(fparaview,"# vtk DataFile Version 3.4.0 \n");
		fprintf(fparaview,"Box file\n");
		fprintf(fparaview,"ASCII\n");
		
		// Writing of the wall objects
		// Using polydata dataset and wrinting of the 8 points defining the simulation box
		fprintf(fparaview,"DATASET POLYDATA\n");
		fprintf(fparaview,"POINTS 8   float\n");
		fprintf(fparaview,"%e %e %e\n",0.0,0.0,0.0);
		fprintf(fparaview,"%e %e %e\n",box_size.x,0.0,0.0);
		fprintf(fparaview,"%e %e %e\n",box_size.x,box_size.y,0.0);
		fprintf(fparaview,"%e %e %e\n",0.0,box_size.y,0.0);
		fprintf(fparaview,"%e %e %e\n",0.0,0.0,box_size.z);
		fprintf(fparaview,"%e %e %e\n",box_size.x,0.0,box_size.z);
		fprintf(fparaview,"%e %e %e\n",box_size.x,box_size.y,box_size.z);
		fprintf(fparaview,"%e %e %e\n",0.0,box_size.y,box_size.z);

		// Firstly, compute the number of planes and then write it
		nwall=0;
		for(i=1;i<=6;i++)
			{
			if(state_obj[i]!=0)
				{
				nwall++;
				}
			}


		fprintf(fparaview,"POLYGONS %d %d \n",nwall,5*nwall);
		if(state_obj[1]!=0)
			{
			fprintf(fparaview,"4 0 3 7 4 \n");
			}
		if(state_obj[2]!=0)
			{
			fprintf(fparaview,"4 1 2 6 5 \n");
			}
		if(state_obj[3]!=0)
			{
			fprintf(fparaview,"4 0 1 5 4 \n");
			}
		if(state_obj[4]!=0)
			{
			fprintf(fparaview,"4 3 2 6 7 \n");
			}
		if(state_obj[5]!=0)
			{
			fprintf(fparaview,"4 0 1 2 3 \n");
			}

		if(state_obj[6]!=0)
			{
			fprintf(fparaview,"4 4 5 6 7 \n");
			}
		

		fprintf(fparaview,"\n");
		
		fclose(fparaview);
		} 
	else 
		{
		perror(Nvtkfile);				
		fprintf(stderr," Check if the vtk file is used by another program \n");		
		exit(ierror=errno);
		}
	
	}


