/**
*\file write_vtk.c
*\brief body of the function write_vtk
*/

#include "write_vtk.h"
#include "def_global_variables.h"

void write_vtk(const char *Nvtkfile)
	{	
	FILE *fparaview;
	int nppcell=1; //!< Number of points by cell
	unsigned int ip;  //!< integer for loop on the particles
    
	fparaview = fopen(Nvtkfile,"w");
	
	if (fparaview != NULL) 
		{
		// Writing typical header for a vtk file
		fprintf(fparaview,"# vtk DataFile Version 3.4.0 \n");
		fprintf(fparaview,"demGCE code\n");
		fprintf(fparaview,"ASCII\n");
		
		// Using unstructured grid dataset because particle packing consists of arbitrary dispositions. 
		// Each point corresponds to one cell compsed of only one particle coordinate 
		fprintf(fparaview,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(fparaview,"POINTS %d   float\n",npart);		
		for (ip=1;ip<=npart;ip++)
			{
			fprintf(fparaview,"%e %e %e\n",particle[ip].Ri.x,particle[ip].Ri.y,particle[ip].Ri.z);
			}

		fprintf(fparaview,"CELLS %d %d \n",npart,2*npart);
		for(ip=1;ip<npart+1;ip++)
			{
			fprintf(fparaview,"%d %d \n",nppcell,ip-1);
			}
		// Cell type is one by default (value of VTK_VERTEX which corresponds to 3D point)		
		fprintf(fparaview,"CELL_TYPES %d \n",npart);
		for(ip=1;ip<npart+1;ip++)
			{
			fprintf(fparaview,"%d \n",1);
			}		   
		fprintf(fparaview,"POINT_DATA %d \n",npart);
		
		// Writing scalar information 
		// Radius of the particles 
		fprintf(fparaview,"\n");		 
		fprintf(fparaview,"SCALARS Ri float \n");	 
		fprintf(fparaview,"LOOKUP_TABLE default \n");	
		
		for(ip=1;ip<=npart;ip++)
			{
			// hmax is equal to zero for solid interaction and negligible for the DLVO forces
			if(forceHCSS){fprintf(fparaview,"%lf\n",particle[ip].radius+hmax/2);}
			else{fprintf(fparaview,"%lf\n",particle[ip].radius);}
			}
		
		//Id of the particles
		fprintf(fparaview,"\n");
		fprintf(fparaview,"SCALARS IDP(-) int \n");
		fprintf(fparaview,"LOOKUP_TABLE default \n");			
		for(ip=1;ip<npart+1;ip++)
			{
			fprintf(fparaview,"%d\n",ip);
			}

		//Id of the cluster to which particles are belonging
		if(Fclust)
			{
			fprintf(fparaview,"\n");
			fprintf(fparaview,"SCALARS IDClust(-) int \n");
			fprintf(fparaview,"LOOKUP_TABLE default \n");
			for(ip=1;ip<npart+1;ip++)
				{
				fprintf(fparaview,"%d\n",particle[ip].clust);
				}
			}
				
		// Linear velocity of the particles		
		fprintf(fparaview,"\n");
		fprintf(fparaview,"VECTORS Vi(u/s) float \n");
		for (ip=1;ip<=npart;ip++)
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


