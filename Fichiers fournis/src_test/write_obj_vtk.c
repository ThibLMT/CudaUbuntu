/**
*\file write_obj_vtk.c
*\brief body of the function write_obj_vtk
*/

#include "write_obj_vtk.h"
#include "def_global_variables.h"

void write_obj_vtk(const char *Nvtkfile)
	{	
	FILE *fparaview;
	int nppcell=1; //!< Number of points by cell
	long int ip;  //!< integer for loop on the particles
	int nobjmin,nobjmax;  //!< Range of identities of objects except planes and cylinder
	int nobjsphere;  //!< Number of sphere objects

	fparaview = fopen(Nvtkfile,"w");
	nobjmin=npart+nsurf+1;
	nobjmax=npart+nobj;
	nobjsphere=nobjmax-nobjmin+1;

	if (fparaview != NULL) 
		{
		// Writing typical header for a vtk file
		fprintf(fparaview,"# vtk DataFile Version 3.4.0 \n");
		fprintf(fparaview,"demGCE code\n");
		fprintf(fparaview,"ASCII\n");
		
		// Writing of the sphere objects
		// Using unstructured grid dataset because particle packing consists of arbitrary dispositions. 
		// Each point corresponds to one cell compsed of only one particle coordinate 
		fprintf(fparaview,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(fparaview,"POINTS %d   float\n",nobjsphere);
		for (ip=nobjmin;ip<=nobjmax;ip++)
			{
			fprintf(fparaview,"%e %e %e\n",particle[ip].Ri.x,particle[ip].Ri.y,particle[ip].Ri.z);	
			}

		fprintf(fparaview,"CELLS %d %d \n",nobjsphere,2*nobjsphere);
		//for(ip=nobjmin;ip<=nobjmax;ip++)
		for(ip=1;ip<=nobjsphere;ip++)
			{
			fprintf(fparaview,"%d %ld \n",nppcell,ip-1);
			}
		// Cell type is one by default (value of VTK_VERTEX which corresponds to 3D point)		
		fprintf(fparaview,"CELL_TYPES %d \n",nobjsphere);
		for(ip=nobjmin;ip<=nobjmax;ip++)
			{
			fprintf(fparaview,"%d \n",1);	
			}		   
		fprintf(fparaview,"POINT_DATA %d \n",nobjsphere);
		
		// Writing scalar information 
		// Radius of the particles 
		fprintf(fparaview,"\n");		 
		fprintf(fparaview,"SCALARS Ri float \n");	 
		fprintf(fparaview,"LOOKUP_TABLE default \n");	
		
		for(ip=nobjmin;ip<=nobjmax;ip++)
			{
			fprintf(fparaview,"%lf\n",particle[ip].radius);
			}
		
		//Id of the particles
		fprintf(fparaview,"\n");
		fprintf(fparaview,"SCALARS IDP(-) int \n");	 
		fprintf(fparaview,"LOOKUP_TABLE default \n");			
		for(ip=nobjmin;ip<=nobjmax;ip++)
			{
			fprintf(fparaview,"%d\n",(int)ip);
			}
				
		// Linear velocity of the particles		
		fprintf(fparaview,"\n");
		fprintf(fparaview,"VECTORS Vi(u/s) float \n");
		for (ip=nobjmin;ip<=nobjmax;ip++)
			{
			fprintf(fparaview,"%e %e %e\n",particle[ip].Vi.x,particle[ip].Vi.y,particle[ip].Vi.z);
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


