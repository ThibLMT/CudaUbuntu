/**
 *\file main.c          \brief main function                               
 *\version $Revision: 1.65 $  $Date: 2013/07/16 18:07:23 $ 
 */
/**************************************************************************
 *  demGCE in-house code					  	   *
 *  Discrete element modeling of Granular materials for Civil Engineering  *
 *				   			                   *
 * Developped by the civil and Environmental Engineering Department        *
 * Research centre of Ecole Mines Douai, France                            *         	                   
 * IDDN.FR.010.0118996.000.R.A.2013.035.30000        	                   *
 * 
 * Parallelisation du code - Maison de la Simulation *
 * version du 25/09/2013 *
 ***************************************************************************/ 

// General definitions
//********************
// Standard library inclusions

#include <time.h>
#include <omp.h>


// Project inclusion
#include "def_types.h"
#include "def_global_variables.h"

// Declaration of the functions
//*****************************
// Insertion of the header files linking extern functions with corresponding source files (*.c) 
#include "list_header_files.h"

// Functions used to generate sphere packing and boundary conditions
void boundary_conditions(void);
void particles_properties(void);
struct sphere wall_properties(struct sphere part_wall);

// Functions of the DEM computation -> integrate contains the DEM loop and check used to detect the contacts

double num_coor();
void output_data(int istep,double timet);


// Function to run program
void end_program(void);
static void print_about(void);

/**
 *\fn int main(int argc, char *argv[])
 *\brief main function or the program
 *\param argc number of arguments supplied to the program
 *\param *argv  Arguments supplied to the program	
 *\return EXIT_SUCCESS  - Normal stop of the program
 */

int main(int argc, char *argv[])
{
  char Nlogfile[50]="logfile";
  char Nmicrofile[50]="init_micro";
  // Initialization of some global variables
  //** Flag variables
  Nparafile="params.in";	
  int i,j,k,k1;
  int nb_taches;
  int iter,niter,imicro;
  char filename[50];
  discrete_elt *particle; 
	
	
  // Initialization of Ierror
  ierror=EXIT_SUCCESS;
	
  // Function called by exit() function for both failure or sucess cases
  atexit(end_program);
	
  // Treatement of the main() arguments
  fprintf(stdout,"***********************************\n");
  fprintf(stdout,"* demGCE code \n");
  if (argc>1)
    {
      if (!strcmp (argv[1],"-about"))
	{
	  print_about();
	  return 0;
	}
      else if (!strcmp (argv[1],"-file"))
	{
	  Nparafile=argv[2];
	}
      else 
	{
	  fprintf(stdout,"%s is not a valid option \n",argv[1]);
	  fprintf(stdout,"list of available options : \n");
	  fprintf(stdout,"  -file Name_of_parameter_file \n");
	  fprintf(stdout,"  -about \n");
	  return 0;
	}
    }
	
  // Open first the log file (history of simulation)	
  flogfile=fopen(Nlogfile,"w");
  // Test the logfile
  if(flogfile != NULL)
    {
    fprintf(flogfile,"***********************************\n");
    fprintf(flogfile,"* demGCE code : log file \n");
    fprintf(flogfile,"  Compilation "__DATE__" a "__TIME__"\n");
    }
  else	
	{
    perror(Nlogfile);
    exit(ierror=errno);
    }

  // Test openmp and print the number of threads
  #ifdef _OPENMP
    #pragma omp parallel
    {
    #pragma omp master
      {
         nb_taches = omp_get_num_threads() ;
      }
     }
    fprintf (flogfile,"** parallel execution : %d thread(s) **", nb_taches) ;
    fprintf (stdout,"** parallel execution : %d thread(s) **", nb_taches) ;
   #endif
   
  // Define the size of the particle systems 
  // * Get   the number of particle
  geom.nb_part=microfile_read_npart(Nmicrofile);
  // * Get the number of boundary conditions
  geom.nb_bc=20;
  geom.max_cont_per_part=(int)MAXCONT;
  read_geom_param(&geom);
  
  
  // Allocate particle array
  int nb_elements=geom.nb_part + geom.nb_bc+1;
  
  particle=allocate_particle(nb_elements,sizeof(discrete_elt));
  
  // Initialize particle array via openmp
  initialize_particle(particle,geom);
  
  // Read "micro" particle, this file contains the coordinate of the particles
  microfile_read_particle(Nmicrofile,particle,geom);
  microcontfile_read_contact(Nmicrofile,particle,geom);
  
  //********************
  // Initialize parameters
  // Bulk parameters
  read_table_mat(&prop_mat_part);
  
  // Friction parameters
  prop_mat_part.mu_gg=0.3; //!< Friction coefficient grain-grain
  prop_mat_part.mu_gw=0.3; //!< Friction coefficient grain-wall
  // -- Rolling resistant parameters
  prop_mat_part.mu_roll_gg=0.01; //!< Rolling friction coefficient grain-grain
  prop_mat_part.mu_roll_gw=0.01;
  // Set bulk forces
  gravity.x=0.0;
  gravity.y=0.0;
  gravity.z=-9.81;  // m.s-2
  
  //adimention of length 
  adi_params(&prop_mat_part,geom);
	
  // give properties for particle
  for(i=1;i<=geom.nb_part;i++)
    {
	give_properties_particle(&particle[i],geom.unity,prop_mat_part);
	}
   
  // set parameter for the plane
  i=geom.nb_part+1;
  give_properties_particle(&particle[i],geom.unity,prop_mat_part);
  particle[i].radius=1000000.0;
	particle[i].Ri.x=geom.sizex/2.0;
	particle[i].Ri.y=geom.sizey/2.0;
	particle[i].Ri.z=0.0 + -1.0 * particle[i].radius;
  
	
  // Start DEM computation
  // Allocation of subdomain backgrid
  backgrid=allocation_backgrid(geom);
  backgrid_insert=allocation_backgrid_insert(geom);
  
  // Initialization via openmp loop of backgrid
  initialize_backgrid(backgrid,backgrid_insert,geom);
  
  // Add boundary condition
  for(i=0;i<geom.sizex;i++)
	{
	  for(j=0;j<geom.sizey;j++)
	    {
	      set_id_backgrid(i,j,0,geom.nb_part+1,backgrid,backgrid_insert,geom);
	    }
			
	}  
  geom.deltat=0.000001;
  niter=1000000;
  imicro=0;
  iter=0;
  microfile_write("micro_ini",particle,geom);
  do{
	 
	//***********************************************************
    //reset forces for all the particles and the boundary conditions (objects and surfaces) 
    #pragma omp parallel for shared(particle,geom) private(i) schedule(static)
     for(i=1;i<=geom.nb_part+geom.nb_bc;i++)
	    {
	    set_vect_0(&particle[i].Fi);
	    set_vect_0(&particle[i].Mi);
	    }

	 // Loop to compute the volumic forces 
    #pragma omp parallel for shared (particle,geom) private(i) firstprivate(gravity) default(none)
      for(i=1;i<=geom.nb_part;i++)
	    {	
	  // Gravity force
	  particle[i].Fi.x+=gravity.x*particle[i].mass;
	  particle[i].Fi.y+=gravity.y*particle[i].mass;
	  particle[i].Fi.z+=gravity.z*particle[i].mass;
			
	  // Hydrodynamic forces
	  //... 
        }
        
     // add id in subdomain backgrid
     #pragma omp  parallel for shared(backgrid,backgrid_insert,geom) private(i) schedule(static) 
      for(i=1;i<=geom.nb_part;i++)
	{
	  
	  insert_sph_backgrid(particle[i].Ri.x,particle[i].Ri.y,particle[i].Ri.z,particle[i].radius,i,backgrid,backgrid_insert,geom);
	}    
      
    // Loop to compute the contact forces
    #pragma omp parallel for shared (particle,geom,backgrid) private(i) schedule(static)
      for(i=1;i<=geom.nb_part;i++)
	{
	cal_interaction_forces(i,particle,geom,backgrid);
	  
	}
	
	
    if(iter%100000==0)
    {
        sprintf(filename,"micro_%04d",imicro);
        printf("micro iter %d %d \n",iter,imicro);
        microfile_write(filename,particle,geom);
        imicro++;
    }
    
    //... 
    //**********
      //  Update particle position for the next ime step - loop on the particle
      // Semi-implicite euler scheme
      // ai(t)=1/mi*sum(Fext->i)
      // vi(t+dt)=vi(t)+ai(t)*dt
      // xi(t+dt)=xi(t)+vi(t+dt)*dt
      // awi(t)=1/Ii*sum(Text->i)
      // wi(t+dt)=wi(t)+awi(t)*dt  
        #pragma omp parallel for shared(particle,geom,unity) private(i) schedule(static)
          for(i=1;i<=geom.nb_part;i++)
          {
          update_particle(&particle[i],geom);
          }

        // Reset backgrid
        memset(backgrid,0,geom.sizex*geom.sizey*geom.sizez*geom.sizel*sizeof(unsigned int));
        memset(backgrid_insert,0,geom.sizex*geom.sizey*geom.sizez*sizeof(int));

          // Ajout plan z=0 dans backgrid
          for(i=0;i<geom.sizex;i++)
        {
          for(j=0;j<geom.sizey;j++)
            {
              set_id_backgrid(i,j,0,geom.nb_part+1,backgrid,backgrid_insert,geom);
            }

        }


      iter++;
  }
  while(iter<=niter);
  microfile_write("micro_fin",particle,geom);


  
		

   exit(errno);

    
  return EXIT_SUCCESS;
}


































//*****************	
// functions
//************
// Functions used to generate sphere packing and boundary conditions
//******************************************************************

/**
 *\fn void boundary_conditions(void)
 *\brief function to set up the boundary conditions
 *\param void	
 *\return void
 */


void boundary_conditions(void)
{
  vector vnwall,Xp;
  unsigned int i, isurf;
  fprintf(stdout,"\n  List of the boundary conditions\n");
  fprintf(stdout,"=================================\n");
  fprintf(flogfile,"\n  List of the boundary conditions\n");
  fprintf(flogfile,"=================================\n");
	
  box_size.x=(double)syssizex;
  box_size.y=(double)syssizey;
  box_size.z=(double)syssizez;
	
  if(state_obj[1]==1)
    {
      // add wall x=0	 isurf=2
      isurf=2;
      Xp.x=0;
      Xp.y=syssizey/2.0;
      Xp.z=syssizez/2.0;

      vnwall.x=-1.0;
      vnwall.y=0.0;
      vnwall.z=0.0;
      particle[npart+isurf]=add_wall(npart+isurf,Xp,vnwall);
      particle[npart+isurf]=wall_properties(particle[npart+isurf]);
      box_size.x-=0.5;
    }	
	
  if(state_obj[2]==1)
    {
      // add wall x=syssizex-1 isurf=3
      isurf=3;
      Xp.x=syssizex-1;
      Xp.y=syssizey/2.0;
      Xp.z=syssizez/2.0;
						
      vnwall.x=1.0;
      vnwall.y=0.0;
      vnwall.z=0.0;
      particle[npart+isurf]=add_wall(npart+isurf,Xp,vnwall);
      particle[npart+isurf]=wall_properties(particle[npart+isurf]);
      box_size.x-=0.5;
    }	
	
  if(state_obj[3]==1)
    {			
      // add wall y=0  isurf=4
      isurf=4;
      Xp.x=syssizex/2.0;
      Xp.y=0.0;
      Xp.z=syssizez/2.0;
      vnwall.x=0.0;
      vnwall.y=-1.0;
      vnwall.z=0.0;
      particle[npart+isurf]=add_wall(npart+isurf,Xp,vnwall);
      particle[npart+isurf]=wall_properties(particle[npart+isurf]);
      box_size.y-=0.5;
    }
	
  if(state_obj[4]==1)
    {
      // add wall y=syssizey-1 isurf=5
      isurf=5;
      Xp.x=syssizex/2.0;
      Xp.y=syssizey-1;
      Xp.z=syssizez/2.0;
      vnwall.x=0.0;
      vnwall.y=1.0;
      vnwall.z=0.0;
      particle[npart+isurf]=add_wall(npart+isurf,Xp,vnwall);
      particle[npart+isurf]=wall_properties(particle[npart+isurf]);
      box_size.y-=0.5;
    }
	
  if(state_obj[5]==1)
    {
      // add wall z=0	isurf=1
      isurf=1;
      Xp.x=syssizex/2.0;
      Xp.y=syssizey/2.0;
      Xp.z=0.0;
      vnwall.x=0.0;
      vnwall.y=0.0;
      vnwall.z=-1.0;
      particle[npart+isurf]=add_wall(npart+isurf,Xp,vnwall);
      particle[npart+isurf]=wall_properties(particle[npart+isurf]);
      box_size.z-=0.5;
    }
	
  if(state_obj[6]==1)
    {
      // add wall z=syssizez-1 - isurf=6
      isurf=6;
      Xp.x=syssizex/2.0;
      Xp.y=syssizey/2.0;
      Xp.z=syssizez-1;
      vnwall.x=0.0;
      vnwall.y=0.0;
      vnwall.z=1.0;
      particle[npart+isurf]=add_wall(npart+isurf,Xp,vnwall);
      particle[npart+isurf]=wall_properties(particle[npart+isurf]);
      box_size.z-=0.5;
    }
	
  if(state_obj[7]==1)
    {
      // add cylinder along the z axe - isurf=7
      //add_cylinder();
      isurf=7;
      particle[npart+isurf]=wall_properties(particle[npart+7]);
    }

  if(state_obj[8]==1)
    {
      // add cone along the z axe - isurf=8
      isurf=8;
      add_cone(npart+isurf);
      particle[npart+isurf]=wall_properties(particle[npart+isurf]);
    }

  if(state_obj[nsurf+1]==1)
    {
      add_obj_spheres(Nobjfile,npart+nsurf+1,false);
      for(i=npart+nsurf+1;i<=npart+nobj;i++)
	{
	  particle[i]=wall_properties(particle[i]);
	}

    }
  fprintf(stdout,"  system size : %d %d %d \n",syssizex,syssizey,syssizez);
  fprintf(stdout,"  box_size : %5.2f %5.2f %5.2f \n",box_size.x,box_size.y,box_size.z);
  fprintf(stdout,"  Unity : %e m \n",unity);
} 
	
	


	

		
/**
 *\fn void particles_properties(void)
 *\brief function to set the properties of the particles
 *\param void	
 *\return void
 */		
		
void particles_properties(void)
{
  unsigned int k;
  double massred;
  double hpaste;
	
  for(k=1;k<=npart;k++)
    {
      particle[k].mass=4.0/3.0*PI*(particle[k].radius)*(particle[k].radius)*(particle[k].radius)*unity*unity*unity;
      particle[k].mass=particle[k].mass*density;
      particle[k].inertia=0.4*particle[k].mass*particle[k].radius*particle[k].radius*unity*unity;
      if(forceHCSS)
	{
	  // In the case of the ForceHCCS, the mass of soft shell is taken into account 
	  hpaste=hmax;
	  particle[k].mass=4.0/3.0*PI*(particle[k].radius+hpaste/2.0)*(particle[k].radius+hpaste/2.0)*(particle[k].radius+hpaste/2.0)*unity*unity*unity*paste_density;
	  massred=4.0/3.0*PI*(particle[k].radius)*(particle[k].radius)*(particle[k].radius)*unity*unity*unity*(density-paste_density);
	  particle[k].inertia=0.4*massred*particle[k].radius*particle[k].radius*unity*unity;
	  particle[k].inertia+=0.4*particle[k].mass*(particle[k].radius+hpaste/2.0)*(particle[k].radius+hpaste/2.0)*unity*unity;
	  particle[k].mass=particle[k].mass+massred;
	}
		
      particle[k].Yn=kn;
      particle[k].Nu=nu;
      particle[k].Ndamp=cn;
      particle[k].Mu=mu_gg;
      particle[k].Mur=mu_roll_gg;
    }

  // Writing some information in stdout and logfile
	
  fprintf(stdout,"\n\n  Assignment of material properties for the particles \n");
  fprintf(stdout,"====================================================\n");
  fprintf(stdout,"  Only one material\n");
  fprintf(flogfile,"\n\n  Assignment of materialproperties for the particles \n");
  fprintf(flogfile,"====================================================\n");
  fprintf(flogfile,"  Radius %e",particle[1].radius);
  fprintf(flogfile,"\t Mass %e",particle[1].mass);
  fprintf(flogfile,"\t Inerty %e \n",particle[1].inertia);
  fprintf(flogfile,"  Normal reduced modulus %e",particle[1].Yn);
  fprintf(flogfile,"\t Poisson's ratio %1.3e",particle[1].Nu);
  fprintf(flogfile,"\t Normal damping coef. %e \n",particle[1].Ndamp);
  fprintf(flogfile,"\t  Friction coef. %e \n",particle[1].Mu);
  fprintf(flogfile,"  Rolling friction coef. %e",particle[1].Mur);
  fprintf(flogfile,"\n");
} /* Fin du programme particles_properties.c */


/**
 *\fn struct sphere wall_properties(struct sphere part_wall)
 *\brief function to set the properties of the wall
 *\param struct sphere	
 *\return struct sphere
 */		

struct sphere wall_properties(struct sphere part_wall)
{
  // 	Affectation des proprietes physiques, micromecaniques et physico-chimiques							
  part_wall.mass=4.0/3.0*PI*pow(part_wall.radius,3.0)*pow(unity,3.0);
  part_wall.mass=part_wall.mass*density;
  part_wall.inertia=0.4*part_wall.mass*pow(part_wall.radius,2.0)*pow(unity,2.0);
  part_wall.Yn=kn;
  part_wall.Nu=nu;
  part_wall.Ndamp=cn;
  part_wall.Mu=mu_gw;
  part_wall.Mur=mu_roll_gw;
	  						
  // Printing wall characteristics in the log file 
					
		
  fprintf(flogfile," Radius : %e \t",part_wall.radius);
  fprintf(flogfile," Mass : %e \t",part_wall.mass);
  fprintf(flogfile," Inerty : %e \n",part_wall.inertia);
  fprintf(flogfile," Module reduit n : %e \t",part_wall.Yn);
  fprintf(flogfile," Coef. Poisson : %e \t",part_wall.Nu);
  fprintf(flogfile," Amort normale : %e \n",part_wall.Ndamp);
  fprintf(flogfile," Coef. frottement  : %e \n",part_wall.Mu);
  fprintf(flogfile," Coef. roulement   : %e \t",part_wall.Mur);
  if((part_wall.Vi.x!=0.0)||(part_wall.Vi.y!=0.0)||(part_wall.Vi.z!=0.0))
    {fprintf(flogfile,"  Mobil wall\n");}
  else {fprintf(flogfile,"  fix wall\n");}
		
  return part_wall;
		
} /* Fin du programme wallsproperties.c */











// Old check function with id=3 - detect contacts between idpart and other grains

void output_data(int istep,double timet)
{
  int idparti; //!w id of the particle i
  double vit_norm,wit_norm;
  int npwmax;  //!< Identity of the particle with the largest angular velocity and with it coordination number
  double nzcoor;  //!< Averaged coordination number
  vector Rimax,Rimin,Rimean; //!< Position of the particle localized at the extremun of the granular media along each axe
  double vmoy;
  double ecr,ecl,ecp;
  double wmax,vwmax;
  double wmoy;
  vector fobjects;
  double compacite;
  double sxy;
  double sr_supwall;
   
   
   
  // Call analysing fonction 
  vit_norm=0.0;
  wit_norm=0.0;
  wmax=0.0;
  vwmax=0.0;
  npwmax=0;
  vmoy=0.0;
  wmoy=0.0;
  ecr=0.0;
  ecp=0.0;
  ecl=0.0;

  Rimin=(vector){syssizex,syssizey,syssizez};
  Rimax=vect0;
  Rimean=vect0;
  // Compute mean values for the step 1 (only check the interaction forces) 
  for(idparti=1;idparti<=npart;idparti++)
    {
      Rimean.x+=particle[idparti].Ri.x;
      Rimean.y+=particle[idparti].Ri.y;
      Rimean.z+=particle[idparti].Ri.z;
      if(particle[idparti].Ri.x+particle[idparti].radius>Rimax.x){Rimax.x=particle[idparti].Ri.x+particle[idparti].radius;}
      if(particle[idparti].Ri.x-particle[idparti].radius<Rimin.x){Rimin.x=particle[idparti].Ri.x-particle[idparti].radius;}
      if(particle[idparti].Ri.y+particle[idparti].radius>Rimax.y){Rimax.y=particle[idparti].Ri.y+particle[idparti].radius;}
      if(particle[idparti].Ri.y-particle[idparti].radius<Rimin.y){Rimin.y=particle[idparti].Ri.y-particle[idparti].radius;}
      if(particle[idparti].Ri.z+particle[idparti].radius>Rimax.z){Rimax.z=particle[idparti].Ri.z+particle[idparti].radius;}
      if(particle[idparti].Ri.z-particle[idparti].radius<Rimin.z){Rimin.z=particle[idparti].Ri.z-particle[idparti].radius;}
		
      vit_norm=sqrt(particle[idparti].Vi.x*particle[idparti].Vi.x+particle[idparti].Vi.y*particle[idparti].Vi.y+particle[idparti].Vi.z*particle[idparti].Vi.z);
      // Translatiobnal kinetic energy
      ecl+=0.5*particle[idparti].mass*pow(vit_norm,2.0)*pow(unity,2.0);
      // Mean velocity
      vmoy+=vit_norm;
      wit_norm=sqrt(particle[idparti].Wi.x*particle[idparti].Wi.x+particle[idparti].Wi.y*particle[idparti].Wi.y+particle[idparti].Wi.z*particle[idparti].Wi.z);
      //
      wmoy+=wit_norm;
      // Angular kinetic energy
      ecr+=0.5*particle[idparti].inertia*pow(wit_norm,2.0);
      // Potential gravitational energy
      ecp+=-particle[idparti].mass*(particle[idparti].Ri.x*gravity.x+particle[idparti].Ri.y*gravity.y+particle[idparti].Ri.z*gravity.z)*unity;
     
      if(wit_norm>wmax)
	{
	  wmax=wit_norm;
	  vwmax=vit_norm;
	  npwmax=idparti;
	}
    }
    
  Rimean.x=Rimean.x/npart;
  Rimean.y=Rimean.y/npart;
  Rimean.z=Rimean.z/npart;
  vmoy=vmoy/npart;
  wmoy=wmoy/npart;
  nzcoor=num_coor();
  sxy=box_size.x*box_size.y;
  if(ChoixTypeCalcul==3)
    {
      compacite=compacity_z(0.0,Rimax.z,sxy);
    }
  else
    {
      compacite=compacity_z(0.25*Rimax.z,0.75*Rimax.z,sxy);
    }

  // Total force applied on the objects
  if(Fsph_obj)
    {
      fobjects=vect0;
      for(idparti=npart+nsurf+1;idparti<=npart+nobj;idparti++)
	{
	  fobjects.x += particle[idparti].Fi.x;
	  fobjects.y += particle[idparti].Fi.y;
	  fobjects.z += particle[idparti].Fi.z;
	}
    }
  // Average stresses on the boundary conditions
			
  sr_supwall=particle[npart+6].Fi.z/(sxy*unity*unity);
   
  //**************************************
  // Screen printing and writing of the log file writing (screen -> some essential informarion // logfile -> more detailed information
  if(istep==1||istep%ndodisplay==0||istep==niter)
    {
      // print screen 
      fprintf(stdout,"\n");
      fprintf(stdout,"  Current Iteration : %d      \n",istep);
      fprintf(stdout,"    zmoy %lf zmin %lf zmax %lf \n",Rimean.z,Rimin.z,Rimax.z);
      fprintf(stdout,"    Overlap mean %e min %e max %e \n",overlap[0],overlap[1],overlap[2]);
      fprintf(stdout,"    Indfric %e - Id particle with Wmax %d\n",indfric,npwmax);
			
			
      //********
      // Print some information for the logfile - compute first the coordination number of the largest angular velocity particle
      fprintf(flogfile,"  *************************** \n");
      fprintf(flogfile,"  Current Iteration %d time(s) %e \n",istep,timet);
      fprintf(flogfile,"    <x> %lf xmin %lf xmax %lf \n",Rimean.x,Rimin.x,Rimax.x);
      fprintf(flogfile,"    <y> %lf ymin %lf ymax %lf \n",Rimean.y,Rimin.y,Rimax.y);
      fprintf(flogfile,"    <z> %lf zmin %lf zmax %lf \n",Rimean.z,Rimin.z,Rimax.z);
      fprintf(flogfile,"    vmoy %lf  wmoy %lf \n",vmoy,wmoy);
			
      fprintf(flogfile,"    Particle with the largest angular velocity id %d Wmax %f vwmax %f \n",npwmax,wmax,vwmax);
      //fprintf(flogfile,"    Overlap mean %e - min %e - max %e n1 %d n2 %d \n",overlap[0],overlap[1],overlap[2],(int)overlap[3],(int)overlap[4]);
      //fprintf(flogfile,"    Contact: tot %u pot %u gran %u slide %u - friction index %f\n");
    }
  // First line of results.dat file
  if(istep==1)
    {
      fresultfile=fopen("results.dat","w");
      // In all cases
      fprintf(fresultfile,"%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s","#Times_(s)","Iteration","Zmean_(u)","Zmin_(un)","Zmax_(un)","Compacity","Ecl_(J)","Ecr_(J)","Ecp_(J)","Ncoor_num ");
      /*
      // Depend of the boundary conditions
      if(Fcyl_bond)
      {
      fprintf(fresultfile,"\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s","Diamcyl(u)","sr_latcyl(Pa)","Fcyl.x(N)","Fcyl.y(N)","Fcyl.z(N)");
      }
				
      // Depend of the type of simulation
      if(ChoixTypeCalcul==3||(ChoixTypeCalcul==5))
      {
      fprintf(fresultfile,"\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s","Zplateau_(un)","Vplateau_(un/s)","Fsup.x_(N)","Fsup.y_(N)","Fsup.z_(N)","Finf.x_(N)","Finf.y_(N)","Finf.z_(N)","sr_supwall_(Pa)");
      }
      if(Fsph_obj)
      {
      fprintf(fresultfile,"\t%-12s\t%-12s\t%-12s","Fobj.x_(N)","Fobj.y_(N)","Fobj.z_(N)");
      if(ChoixTypeCalcul==5){fprintf(fresultfile,"\t%-12s","Fcis.y_(N)");}
      }

      if(Fclust)
      {
      fprintf(fresultfile,"\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s","<Npp>","<Npp-pond>","<Df>","<Df-pond>","<Comp>","<Comp-pond>","nc");
      }*/

      fprintf(fresultfile,"\n");
    }
		
  // Printing results
  // For all cases
  fprintf(fresultfile,"%e\t%-12d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%12.3f",timet,istep,Rimean.z,Rimin.z,Rimax.z,compacite,ecl,ecr,ecp,nzcoor);
			
  // Depend of the boundary conditions
  /*
    if(Fcyl_bond)
    {
    fprintf(fresultfile,"\t%e\t%e\t%e\t%e\t%e",diamcyl,sr_latcyl,particle[npart+7].Fi.x,particle[npart+7].Fi.y,particle[npart+7].Fi.z);
    }
			
			
    // Depend of the studied case
    if((ChoixTypeCalcul==3)||(ChoixTypeCalcul==5))
    {
    fprintf(fresultfile,"\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e",z_upwall,particle[npart+6].Vi.z,particle[npart+6].Fi.x,particle[npart+6].Fi.y,particle[npart+6].Fi.z,particle[npart+1].Fi.x,particle[npart+1].Fi.y,particle[npart+1].Fi.z,sr_supwall);
    }
			

    if(Fsph_obj)
    {
    fprintf(fresultfile,"\t%e\t%e\t%e",fobjects.x,fobjects.y,fobjects.z);
    if(ChoixTypeCalcul==5){fprintf(fresultfile,"\t%e",particle[npart+1].Fi.y+fobjects.y);}
    }
    if(Fclust)
    {
    freturn=cluster_properties(freturn);
    fprintf(fresultfile,"\t%e\t%e\t%e\t%e\t%e\t%e\t%12d",freturn[1],freturn[2],freturn[3],freturn[4],freturn[5],freturn[6],(int)freturn[0]);
    }
			
    if((ChoixTypeCalcul==2)&&(forceHCSS))
    {
    fprintf(fresultfile,"%e\t%e\t%e\t%e",particle[npart+2].Fi.y,particle[npart+3].Fi.y,particle[npart+1].Fi.y,particle[npart+2].Vi.y);
    }
  */
  fprintf(fresultfile,"\n");
  fflush(stdout);
  fflush(fresultfile);
	  
}
   

			
					
double num_coor()
{
  int idparti,idpartj;
  double xparti,yparti,zparti,xpartj,ypartj,zpartj,rparti,rpartj;
  double xdist,ydist,zdist,dcontact,dist2;
  int ncoor,icoor;
  double nu_coor;
  ncoor=0;
  //printf("nzoor\n");
  for(idparti=1;idparti<=npart;idparti++)
    {
      xparti=particle[idparti].Ri.x;
      yparti=particle[idparti].Ri.y;
      zparti=particle[idparti].Ri.z;
      rparti=particle[idparti].radius;
      //printf("%d ",idparti);
      icoor=0;
      for(idpartj=1;idpartj<=npart+nobj;idpartj++)
	{
		   
	  xpartj=particle[idpartj].Ri.x;
	  ypartj=particle[idpartj].Ri.y;
	  zpartj=particle[idpartj].Ri.z;
	  rpartj=particle[idpartj].radius;
	  // Need change for
      
	  if(idpartj<=npart)
	    {
	      //Periodic condition
				 
	      // periodic boundary conditions for x2
	      if((xpartj>xparti)&&(xpartj-xparti>xparti+syssizex-xpartj)){xpartj-=syssizex;}
	      else if((xpartj<xparti)&&(xparti-xpartj>xpartj+syssizex-xparti)){xpartj+=syssizex;}
				  
	      // periodic boundary conditions for x2
	      if((ypartj>yparti)&&(ypartj-yparti>yparti+syssizey-ypartj)){ypartj-=syssizey;}
	      else if((ypartj<yparti)&&(yparti-ypartj>ypartj+syssizey-yparti)){ypartj+=syssizey;}
				  
	      // periodic boundary conditions for x2
	      if((zpartj>zparti)&&(zpartj-zparti>zparti+syssizez-zpartj)){zpartj-=syssizez;}
	      else if((zpartj<zparti)&&(zparti-zpartj>zpartj+syssizez-zparti)){zpartj+=syssizez;}

				  
	    }
      
	  xdist=xparti-xpartj;
	  ydist=yparti-ypartj;
	  zdist=zparti-zpartj;
	  dist2=xdist*xdist+ydist*ydist+zdist*zdist;
	  dcontact=rparti+rpartj+hmax;
	  dcontact=dcontact*dcontact;
	  if((dist2 <= dcontact)&&(idparti!=idpartj))
	    {
	      if((particle[idparti].radius>0.0)&&(particle[idpartj].radius>0.0))
		{
		  ncoor++;
		  icoor++;
		  //printf("%d ",idpartj);
		}
	    }
		
	}
	  
      //printf("* %d\n",icoor);
    }
	
  nu_coor=((double)ncoor)/(npart);
  return nu_coor;
	
}
  
  

// If modular program, don't forget to include times.h and stdlib.h
/**
 *\fn void end_program(void)
 *\brief This function is called to stop job. It prints any error message, computes the elapsed time and displays some environnement variables 
 *\param void	
 *\return void
 */
void end_program(void)
{
  clock_t CPU_time;   
  time_t local_time,CPUt_s;
  struct tm CPUt_dhms;
  char *username,*currentpath;


  deallocate_variables();
  // Generate error message in standard error stream (stderr)

  if(ierror!=0)
    {
#if defined(WIN32) || defined(WIN64)
      HANDLE Ecran = GetStdHandle(STD_OUTPUT_HANDLE);
      SetConsoleTextAttribute(Ecran, FOREGROUND_GREEN);
# endif
      fprintf(stderr,"Job aborted : exit status %d \n",ierror);

      fflush(stderr);
    }


  // compute the local et elapsed CPU times (function clock must divided by CLOCKS_PER_SEC)
  CPU_time=clock();
  time(&local_time); 
  CPUt_s=CPU_time/CLOCKS_PER_SEC;
  CPUt_dhms=*localtime(&CPUt_s); 

  // Get some environnement variables
  username=getenv("USERNAME");
  currentpath=getenv("PWD");

  fprintf(stdout,"===============================================================================\n");
  fprintf(flogfile,"===============================================================================\n");	
  fprintf(stdout," Job done: %s \n",ctime(&local_time));
  fprintf(flogfile," Job done: %s \n",ctime(&local_time));
  fprintf(stdout,"Elapsed CPU time : %dD/%dh/%dm/%ds (%li sec)\n",CPUt_dhms.tm_mday-1,CPUt_dhms.tm_hour-1,CPUt_dhms.tm_min,CPUt_dhms.tm_sec,CPUt_s);
  fprintf(flogfile,"Elapsed CPU time : %dD/%dh/%dm/%ds (%li sec)\n",CPUt_dhms.tm_mday-1,CPUt_dhms.tm_hour-1,CPUt_dhms.tm_min,CPUt_dhms.tm_sec,CPUt_s);
  fprintf(flogfile,"User: %s \n",username);
  fprintf(flogfile,"Work directory: %s \n",currentpath);
  fprintf(stdout,"===============================================================================\n");
  fprintf(flogfile,"===============================================================================\n");
  fclose(flogfile);
}

	

/**
 *\fn void print_about(void)
 *\brief function displaying typical program information in terminal
 *\param void	
 *\return void
 */
	
void print_about(void)
{
  fprintf(stdout,"\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"       *********************************************************************************\n");
  fprintf(stdout,"       * demGCE: Discrete element modeling of Granular materials for Civil Engineering *\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"       * Numerical simulation of discrete spherical particles\n");
  fprintf(stdout,"       * Developped by : *\n");
  fprintf(stdout,"       *             Research Centre of Ecole Mines Douai                *\n");
  fprintf(stdout,"       *             Department of Civil and Environnemental Engineering        *\n");
  fprintf(stdout,"       *             Version: $Revision: 1.65 $ - $Date: 2013/07/16 18:07:23 $   *\n");
  fprintf(stdout,"       **************************************************************************\n");
  fprintf(stdout,"\n");
}
