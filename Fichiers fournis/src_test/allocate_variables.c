/**
 *\file allocate_variables.c
 *\brief body of the function allocate_variables
 */


#include "allocate_variables.h"
#include <omp.h>

extern void set_vect_0(vector *vec);

discrete_elt* allocate_particle(int nb_elements,size_t size)
  {
  return malloc( nb_elements * size );
  }

void initialize_particle(discrete_elt *particle,geom_struct geom)
  {
   unsigned int nelement; //!< Number of elements of particle arrays
   unsigned int i,k; //!< Indice for loop
   nelement=geom.nb_part+geom.nb_bc;
  
   fflush(stdout);
   #pragma omp parallel for private(i,k) shared(particle,nelement)
   for(i=0;i<=nelement;i++)
	  {
	  particle[i].Ri.x=0.0;	// vectoru
	  particle[i].Ri.y=0.0;
	  particle[i].Ri.z=0.0;
	  particle[i].Vi.x=0.0;	// vector
	  particle[i].Vi.y=0.0;
	  particle[i].Vi.z=0.0;
	  particle[i].Ai.x=0.0;	// vector
	  particle[i].Ai.y=0.0;	// vecto
	  particle[i].Ai.z=0.0;	// vecto
	  particle[i].Wi.x=0.0;		// vector
	  particle[i].Wi.y=0.0;		// vector
	  particle[i].Wi.z=0.0;		// vector
	  particle[i].Aroti.x=0.0; 	// vector
	  particle[i].Aroti.y=0.0;
	  particle[i].Aroti.z=0.0;
	  particle[i].Fi.x=0.0;   // vector
	  particle[i].Fi.y=0.0;
	  particle[i].Fi.z=0.0;
	  particle[i].Mi.x=0.0;		// vector
	  particle[i].Mi.x=0.0;
	  particle[i].Mi.x=0.0;
	  particle[i].radius=0.0;		// double
	  particle[i].mass=0.0;		// double
	  particle[i].inertia=0.0;	// double
	  particle[i].Yn=0.0;	//double
	  particle[i].Nu=0.0;		// double
	  particle[i].Ndamp=0.0;		// double
	  particle[i].Mu=0.0;	// double
	  particle[i].Mur=0.0;	// double
	  particle[i].next=0;	// int
	  particle[i].clust=0;    // int
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].type[k]=0;}	// I
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].contact[k]=0;}	// Int array
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].ut[k]=0.0;}			// double array
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].ftanold[k].x=0.0;}	//vector array
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].ftanold[k].y=0.0;}	//vector array
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].ftanold[k].z=0.0;}	//vector array
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].nijold[k].x=0.0;}	//vector array
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].nijold[k].y=0.0;}
	  for(k=0;k<geom.max_cont_per_part;k++){particle[i].nijold[k].z=0.0;}
	  
	  }
	  
  }
  
void give_properties_particle(discrete_elt *particle,double unity,material_data properties)
  {
	 double radius;
	 double mass;
	 discrete_elt particle_i;
  // Mass of particle (Kg)
  
  particle_i=*particle;
  radius=particle_i.radius;
  mass=4.0/3.0*PI*radius*radius*radius*properties.density;
  particle_i.mass=mass;
  // Internia of particle (Kg.m2)
  particle_i.inertia=2.0/5.0*mass*radius*radius*unity*unity;
  // Elastic properties	
    particle_i.Yn=properties.E;
    particle_i.Nu=properties.nu;
    particle_i.Ndamp=properties.cn;
    particle_i.Mu=properties.mu_gg;
    particle_i.Mur=properties.mu_roll_gg;
    *particle=particle_i;
    }

void update_particle(discrete_elt *particle,geom_struct geom) 
    {
		double deltat;
		deltat=geom.deltat;
	//**********
      //  Update particle position for the next ime step - loop on the particle
      // Semi-implicite euler scheme
      // ai(t)=1/mi*sum(Fext->i)
      // vi(t+dt)=vi(t)+ai(t)*dt
      // xi(t+dt)=xi(t)+vi(t+dt)*dt
      // awi(t)=1/Ii*sum(Text->i)
      // wi(t+dt)=wi(t)+awi(t)*dt
	
     // ai(t)=1/mi*sum(Fext->i)
	  particle->Ai.x=particle->Fi.x/geom.unity/particle->mass;
	  particle->Ai.y=particle->Fi.y/geom.unity/particle->mass;
	  particle->Ai.z=particle->Fi.z/geom.unity/particle->mass;

	  // vi(t+dt)=vi(t)+ai(t)*dt
	  particle->Vi.x=(particle->Vi.x+particle->Ai.x*deltat);
	  particle->Vi.y=(particle->Vi.y+particle->Ai.y*deltat);
	  particle->Vi.z=(particle->Vi.z+particle->Ai.z*deltat);
			
	  // Computing of angular acceleration -> awi(t)=1/Ii*sum(Text->i)
	  particle->Aroti.x=particle->Mi.x*geom.unity/particle->inertia;
	  particle->Aroti.y=particle->Mi.y*geom.unity/particle->inertia;
	  particle->Aroti.z=particle->Mi.z*geom.unity/particle->inertia;
							
	  // actualization of angular velocity -> wi(t+dt)=wi(t)+awi(t)*dt
	  particle->Wi.x=particle->Wi.x+particle->Aroti.x*deltat;
	  particle->Wi.y=particle->Wi.y+particle->Aroti.y*deltat;
	  particle->Wi.z=particle->Wi.z+particle->Aroti.z*deltat;

	  // actualization of position by checking the periodic condtition 
	  // xi(t+dt)=xi(t)+vi(t+dt)*dt
	  particle->Ri.x+=particle->Vi.x*deltat;		
	  if(particle->Ri.x > geom.sizex){particle->Ri.x-=geom.sizex;}
	  else if(particle->Ri.x < 0.0){particle->Ri.x+=geom.sizex;}
	  
	  particle->Ri.y+=particle->Vi.y*deltat;
	  if(particle->Ri.y > geom.sizey){particle->Ri.y-=geom.sizey;}
	  else if(particle->Ri.y < 0.0){particle->Ri.y+=geom.sizey;}

	  particle->Ri.z+=particle->Vi.z*deltat;
	  if(particle->Ri.z > geom.sizez){particle->Ri.z-=geom.sizez;}
	  else if(particle->Ri.z < 0.0){particle->Ri.z+=geom.sizez;}
}



/*
void allocate_variables(void)
{
  int ia,ja,ka;
  int dimpartarray;

  // allocation of mic array   
  mic = calloc(syssizex, sizeof(unsigned int***));
  for(ia= 0; ia < syssizex; ia++)
    {
      mic[ia]= calloc(syssizey, sizeof(unsigned int**));
      if(mic[ia] == NULL)
	{
	  fprintf(stderr,"Problem during mic array allocation loop ia\n");
	  exit(ierror=EXIT_FAILURE);
	}
      for(ja= 0; ja < syssizey; ja++)
	{
	  mic[ia][ja]= calloc(syssizez, sizeof(unsigned int*));
	  if(mic[ia][ja] == NULL)
	    {
	      fprintf(stderr,"Problem during mic array allocation loop ja\n");
	      exit(EXIT_FAILURE);
	    }
	  for(ka= 0; ka < syssizez; ka++)
	    {
	      mic[ia][ja][ka]= calloc(ncont, sizeof(unsigned int));
	      if(mic[ia][ja][ka] == NULL)
		{
		  fprintf(stderr,"Problem during mic array allocation loop Ka\n");
		  exit(EXIT_FAILURE);
		}
	    }
	}
    }
		
  // allocation of mic_boundary array   
  mic_boundary = calloc(syssizex, sizeof(unsigned int***));
  for(ia= 0; ia < syssizex; ia++)
    {
      mic_boundary[ia]= calloc(syssizey, sizeof(unsigned int**));
      if(mic_boundary[ia] == NULL)
	{
	  fprintf(stderr,"Problem during mic_boundary array allocation loop ia\n");
	  exit(ierror=EXIT_FAILURE);
	}
      for(ja= 0; ja < syssizey; ja++)
	{
	  mic_boundary[ia][ja]= calloc(syssizez, sizeof(unsigned int*));
	  if(mic_boundary[ia][ja] == NULL)
	    {
	      fprintf(stderr,"Problem during mic_boundary array allocation loop ja\n");
	      exit(EXIT_FAILURE);
	    }
	  for(ka= 0; ka < syssizez; ka++)
	    {
	      mic_boundary[ia][ja][ka]= calloc(ncont, sizeof(unsigned int));
	      if(mic_boundary[ia][ja][ka] == NULL)
		{
		  fprintf(stderr,"Problem during mic_boundary array allocation loop Ka\n");
		  exit(EXIT_FAILURE);
		}
	    }
	}
    }
    

  // allocation of mic_insert array   
  mic_insert = calloc(syssizex, sizeof(unsigned int**));
  for(ia= 0; ia < syssizex; ia++)
    {
      mic_insert[ia]= calloc(syssizey, sizeof(unsigned int*));
      if(mic_insert[ia] == NULL)
	{
	  fprintf(stderr,"Problem during mic_insert array allocation loop ia\n");
	  exit(ierror=EXIT_FAILURE);
	}
      for(ja= 0; ja < syssizey; ja++)
	{
	  mic_insert[ia][ja]= calloc(syssizez, sizeof(unsigned int));
	  if(mic_insert[ia][ja] == NULL)
	    {
	      fprintf(stderr,"Problem during mic_insert array allocation loop ja\n");
	      exit(EXIT_FAILURE);
	    }
	}
    }
		
  // allocation of mic_insert_boundary array   
  mic_insert_boundary = calloc(syssizex, sizeof(unsigned int**));
  for(ia= 0; ia < syssizex; ia++)
    {
      mic_insert_boundary[ia]= calloc(syssizey, sizeof(unsigned int*));
      if(mic_insert_boundary[ia] == NULL)
	{
	  fprintf(stderr,"Problem during mic_insert_boundary array allocation loop ia\n");
	  exit(ierror=EXIT_FAILURE);
	}
      for(ja= 0; ja < syssizey; ja++)
	{
	  mic_insert_boundary[ia][ja]= calloc(syssizez, sizeof(unsigned int));
	  if(mic_insert_boundary[ia][ja] == NULL)
	    {
	      fprintf(stderr,"Problem during mic_insert_boundary array allocation loop ja\n");
	      exit(EXIT_FAILURE);
	    }
	}
    }
	
  // allocation of the particle array 
  // npart and nobj are set up un read_param function. We increment the size of array because our loop on particle started from 1. 	
  dimpartarray=npart+nobj+1;
  particle=malloc(dimpartarray*sizeof(struct sphere));	
}
* */
