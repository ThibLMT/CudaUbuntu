/**
 *\file add_obj_spheres.c
 *\brief body of the function add_obj_spheres
 */

#include "def_global_variables.h"
#include "add_obj_spheres.h"


int check(double dxmin, double dymin, double dzmin, double rayon, int idgrain,int drapeau);


int add_obj_spheres(const char *Nobjfile,unsigned int iinit,Flag check_overlap)
{
  // \var iinit number  of readed line
	
  FILE *fobj;
  unsigned int id_sphere,iobj,rnobj;
  int fit_check;
	
  // iinit = 0, read only first line
  // iinit != 0, read coordinate of the objects and set the particule array from i=iinit
  fobj=fopen(Nobjfile,"r");
  if (fobj != NULL)
    {
      fscanf(fobj,"%d",&rnobj);
		
      if(iinit==0)
	{
	  // Just read the nimber of object
	  fclose(fobj);
	  // check if object or not in object file
	  if (rnobj==0)
	    {
	      fprintf(stdout,"  Warning: no object in %s file - rnobj %d \n",Nobjfile,rnobj);
	      exit(EXIT_FAILURE);
	    }			
	  return rnobj;
	}
				
      for(id_sphere=iinit;id_sphere<=(nobj+npart);id_sphere++)
	{
	  fscanf(fobj,"\n");
	  fscanf(fobj,"%d",&iobj);
	  fscanf(fobj,"%lf",&particle[id_sphere].Ri.x);
	  fscanf(fobj,"%lf",&particle[id_sphere].Ri.y);
	  fscanf(fobj,"%lf",&particle[id_sphere].Ri.z);

	  fscanf(fobj,"%lf",&particle[id_sphere].Vi.x);
	  fscanf(fobj,"%lf",&particle[id_sphere].Vi.y);
	  fscanf(fobj,"%lf",&particle[id_sphere].Vi.z);
	  fscanf(fobj,"%lf",&particle[id_sphere].Wi.x);
	  fscanf(fobj,"%lf",&particle[id_sphere].Wi.y);
	  fscanf(fobj,"%lf",&particle[id_sphere].Wi.z);
	  fscanf(fobj,"%lf",&particle[id_sphere].radius);

	 
	  // Placement of the sphere object in the mic array (set the identity in the voxel)
	  printf("error cleanopenmp : link addsphere to backgrid_n");
	  exit(1);
	  //check(particle[id_sphere].Ri.x,particle[id_sphere].Ri.y,particle[id_sphere].Ri.z,particle[id_sphere].radius+hmax/2.0,id_sphere,1);

	  // Initialization of the forces and torques
	  particle[id_sphere].Fi=vect0;
	  particle[id_sphere].Mi=vect0;

	}
      fclose(fobj);
    }
  else
    {
      perror(Nobjfile);
      exit(ierror=errno);
    }
	
  fprintf(stdout,"Object : id %d to %d - number of objjects %d - file : %s \n",iinit,nobj+npart,rnobj,Nobjfile);
  fprintf(flogfile,"  ** Object file : %s \n",Nobjfile);
  fprintf(flogfile,"  Number of the objects %d \n",nobj);
  return 0;
}


