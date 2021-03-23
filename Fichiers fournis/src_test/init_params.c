/**
*\file init_params.c
*\brief body of the function read_param
*/
#include <string.h>


#include "init_params.h"

void read_geom_param(geom_struct *geom_para)
  {
  geom_struct geom_para_read;
  char **tab_char;  //!< array of string
  int i,k,nstring;
  // Initialization and allocation
  geom_para_read=*geom_para;
  tab_char=malloc(50*sizeof(char*));  
  for(i=0;i<50;i++) {tab_char[i]=malloc(256*sizeof(char));}
  // Read string "tab_mat" and grep info
  
  nstring=get_line_char("syssize",tab_char);
  
  // Read syssize of the system alon x,y,z, define the limit poistion and the dimension of the box simulation and the sud domain backgrid
  for(k=0;k<nstring;k++) 
    {
    if (!strcmp ((const char *) tab_char[k], "syssize"))
      {
		geom_para_read.sizex = atoi (tab_char[k + 1]);
		geom_para_read.sizey = atoi (tab_char[k + 2]);
		geom_para_read.sizez = atoi (tab_char[k + 3]);
		break;
	  }
	}
	 // Read size of the system along l (max number of particles can be put in each voxel of backgrid)
  nstring=get_line_char("syssize_l",tab_char);
  if (!strcmp ((const char *) tab_char[0], "syssize_l")) {geom_para_read.sizel = atoi (tab_char[1]);}
	
  // Read unity length Conversion
  nstring=get_line_char("unity",tab_char);
  if (!strcmp ((const char *) tab_char[0], "unity")) {geom_para_read.unity = atof (tab_char[1]);}
  *geom_para=geom_para_read;
		
  }
  
 void adi_params(material_data *prop_mat_part,geom_struct geom)
   {
	material_data tab_mat;
	tab_mat=*prop_mat_part;
	// E Young modulus Pa (N.m-2 --> N.unity-2)
	tab_mat.E=tab_mat.E*geom.unity*geom.unity;
	// Bulk density
	tab_mat.density=tab_mat.density*geom.unity*geom.unity*geom.unity;
	*prop_mat_part=tab_mat;
	
	}


void read_table_mat(material_data *tab_mat_part)
  {
  material_data mat_params;
  char **tab_char;  //!< array of string (50 elements)
  int i,k,nstring;
  
  // Initialization and allocation
  mat_params=*tab_mat_part;
  tab_char=malloc(50*sizeof(char*));  
  for(i=0;i<50;i++) {tab_char[i]=malloc(256*sizeof(char));}
  // Read string "tab_mat" and grep info
  nstring=get_line_char("tab_mat",tab_char);
  
  for(k=0;k<nstring;k++) 
    {printf(" %d  %s \n",k,tab_char[k]);
    if (!strcmp ((const char *) tab_char[k], "density"))
		mat_params.density = atof (tab_char[k + 1]);
	if (!strcmp ((const char *) tab_char[k], "E"))
		mat_params.E = atof (tab_char[k + 1]);
	if (!strcmp ((const char *) tab_char[k], "nu"))
		mat_params.nu = atof (tab_char[k + 1]);
	if (!strcmp ((const char *) tab_char[k], "cn"))
		mat_params.cn = atof (tab_char[k + 1]);
    } 
  
  *tab_mat_part=mat_params;
  }
  
int get_line_char(char *detect,char **test)
  {
  char *Nparafile;
  FILE *fparafile;
  static const char separateur[] = " ,()=\t\n[]";
  int nbjalon = 0;
  char *s;
  char ligne[256];
  char *jalon[256];
  int i = 0;
  int nstring,lstring;
  
  Nparafile="params.in";
  
  fparafile = fopen (Nparafile, "r");
  if(fparafile != NULL)
	{
	while(!feof (fparafile))
		{
		fgets (ligne, 256, fparafile);
		if (feof (fparafile)) break;
		// First call of strtok function
		s = (char *) strtok (ligne, (const char *) separateur);
		nbjalon = 0;
		jalon[0] = s;
		while (s != NULL)
		  {
		  jalon[nbjalon++] = s;
		  s = (char *) strtok (NULL, (const char *) separateur);
		  }
		jalon[nbjalon] = 0;			
			
		// analyse de la ligne
		i = 0;
		if (!strcmp ((const char *) jalon[i],detect))
		 {
		 while (i < nbjalon)
			{
			if ((!strncmp ((const char *) jalon[i], "#", 1))||(!strncmp ((const char *) jalon[i], "!", 1))||(!strncmp ((const char *) jalon[i], ";", 1)))
			   {nstring=i;break;}
			i++;
		    }
		   break;
		   }
		 
	    
	    }
  
    fclose (fparafile);
    }  
  else
	{
	fprintf(stderr,"Error! No parameter file (%s) \n",Nparafile);
	exit(errno);
	}
  for(i=0;i<nstring;i++){strcpy(test[i],jalon[i]);}
  return nstring + 1;
  }

			
/*
 * void read_param(void)
	{    
	int i,j,k1;
	double shear_rate;		//!< Shear rate of fluid if Stokes forces
	static const char separateur[] = " ;,()=\t\n[]"; 
	char *jalon[64];
	int nbjalon = 0;
	Flag init_gravity,forcetangential;
	char *s;
	char ligne[256];
	 //	Step 1 : INITIALIZATION of values
	// Result iteration increment by default
	ndodisplay=1000; 
	ndowritemicro=10000; 
	ndowritevtk=10000; 
	ndowrite=500;

	// Flag by default
	forceVdw=false;
	forceElectro=false;
	forceFluid=false;
	forceContact=true;
	forcetangential=false;
	torqueRollRes=false;
	init_gravity=true;
	Freadmicro=false;   //<! Read initial microstructure
	Fclust=false;
	Fsph_obj=false;
	Fcyl_bond=false;
	Fcyl_mov=false;
	Fwall=false;
	Fwall_mov=false;
	Fupplane_conf=false;
	Ftriaxial=false;
	// Preinitialization
	nsurf=20;  // nsurf -> Default number of surfaces for the  boundary conditions (spherical  objects not included)
	tan_force_model=1;
	z_upwall_ini=0.0;
	mu_gg=0.0;
	mu_gw=0.0;
	// Number of object 8 by default 6 planes + one cylinder
	shear_rate=0.0;
	FILE *fichierPARAM;
	fichierPARAM=NULL;
  
	// name of initial file parameters for simulation
   	fichierPARAM = fopen (Nparafile, "r");
   	if(fichierPARAM != NULL)
		{	
		fprintf(stdout,"\n  Read parameter file: %s \n",Nparafile );		
		fprintf(stdout,"=============================================\n");
		// Reading and analysing of parameter file line per line 
		while (!feof (fichierPARAM))
			{			
			fgets (ligne, 256, fichierPARAM);
			if (feof (fichierPARAM))
				break;
			// First call of strtok function
			s = (char *) strtok (ligne, (const char *) separateur);
			nbjalon = 0;
			jalon[0] = s;
			while (s != NULL)
				{
				jalon[nbjalon++] = s;
				s = (char *) strtok (NULL, (const char *) separateur);
				}
			jalon[nbjalon] = 0;			
			
			// analyse de la ligne
			i = 0;

			while (i < nbjalon)
				{
				// test to check the fisrt letter or parameter is a capital letter or not. Disable for the moment
				if(i==0)
					{
					
					if(((strncmp ((const char *) jalon[i], "A", 1))>=0)&&((strncmp ((const char *) jalon[i], "Z", 1))<=0))
					fprintf(stdout,"Warning: '%s' parameter should be write in small letter \n",jalon[i]);
					}
					
				if ((!strncmp ((const char *) jalon[i], "#", 1))||(!strncmp ((const char *) jalon[i], "!", 1)))
					break;
												
				              
				// Section 1 : Definition of the packing and boundary conditions   
				
				//new : syssize=[10][22][50]
				if (!strcmp ((const char *) jalon[i], "syssize"))
					{
					syssizex = atoi (jalon[i + 1]);
					syssizey = atoi (jalon[i + 2]);
					syssizez = atoi (jalon[i + 3]);
					if((syssizex<=0)||(syssizey<=0)||(syssizez<=0))
						{
						fprintf(stdout,"Warning : problem in syssize %d %d %d\n",syssizex,syssizey,syssizez);
						exit(EXIT_FAILURE);
						}
					// set ncont parameter needed for mic array allocation
					ncont=50;
					}
				
				//function : PackingGeneration=method(Npart=3000,rpart)									
				if (!strcmp ((const char *) jalon[i], "PackingGeneration"))
					{
					if (!strcmp ((const char *) jalon[i + 1], "RSA")) { generationmodel=1;}
					else if (!strcmp ((const char *) jalon[i + 1], "CCA")){ generationmodel=2;}
					else if (!strcmp ((const char *) jalon[i + 1], "PCA")){ generationmodel=3;}
					
					if (!strcmp ((const char *) jalon[i + 2], "npart"))
						{
						npart = atoi (jalon[i + 3]);
						}
					if (!strcmp ((const char *) jalon[i + 4], "rpart"))
						{
						rpart = atof (jalon[i + 5]);
						}
					
						
					}
				
				//function readmicro : readmicro (micro_file,contact_file)
				if (!strcmp ((const char *) jalon[i], "readmicro"))
					{
					Freadmicro=true;
					strncpy (FichMicroInitiale, (const char *) jalon[i+1],strlen ((const char *) jalon[i+1]));
					strncpy (FichContInitiale, (const char *) jalon[i+2],strlen ((const char *) jalon[i+2]));
					}
					
				//Graine generateur
				//#nseed=-32
				if (!strcmp ((const char *) jalon[i], "nseed"))
					{
					nseed = atof (jalon[i + 1]);
					seed=(&nseed);
					
					if(nseed>=0){fprintf(stdout,"Warning, nseed should be negative");}
					} 
				
				// Unity of the system
				if (!strcmp ((const char *) jalon[i], "unity"))
					{
					unity = atof (jalon[i + 1]);
					} 	
				
				//Choix des conditions aux limites	
				// rigidwall
				if (!strcmp ((const char *) jalon[i], "rigidwall"))
					{
					for (k1=1;k1<=6;k1++)
						{
						state_obj[k1]= atof (jalon[i + k1]);
						}
					break;
					}
				
				// cylinder
				if (!strcmp ((const char *) jalon[i], "cylinder"))
					{
					state_obj[7]=1;
					if (!strcmp ((const char *) jalon[i+1], "diameter"))
						{
						diamcyl = atof (jalon[i + 2]);
						}
					}
				
				// cone
				if (!strcmp ((const char *) jalon[i], "cone"))
					{
					state_obj[8]=1;
					if (!strcmp ((const char *) jalon[i+1], "height"))
						{
						hcone = atof (jalon[i + 2]);
						}
					if (!strcmp ((const char *) jalon[i+3], "rbot"))
						{
						rbcone = atof (jalon[i + 4]);
						}
					if (!strcmp ((const char *) jalon[i+5], "rtop"))
						{
						rtcone = atof (jalon[i + 6]);
						}
					}
				
				// Add spherical objects
				if (!strcmp ((const char *) jalon[i], "add_objects_spheres"))
					{
					state_obj[nsurf+1]=1;
					strncpy (Nobjfile, (const char *) jalon[i+1],strlen ((const char *) jalon[i+1]));
					}
					
			
				
				// Section 2 : DEFINITION OF THE INTERACTIONS     
				// Gravity
				// By default, gravity along z axis gz = 9.81 m.s-2
				// #Gravity(gx,gy,gz)

				if (!strcmp ((const char *) jalon[i], "Gravity"))
					{
					gravity.x=atof (jalon[i + 1]);
					gravity.y=atof (jalon[i + 2]);
					gravity.z=atof (jalon[i + 3]);
					init_gravity=false;
					}
					
				// Contact force laws
				// - Normal force
				// Disable the contact force
				if (!strcmp ((const char *) jalon[i], "no_contact_force"))
					{
					forceContact=false;
					}
				
				// CARACTERISTIQUES DU MATERIAU
				if (!strcmp ((const char *) jalon[i], "density"))
					{
					density = atof (jalon[i + 1]);
					} 
		  
				if (!strcmp ((const char *) jalon[i], "ModuleYoung"))
					{
					kn = atof (jalon[i + 1]);
					}
		  
				if (!strcmp ((const char *) jalon[i], "nu"))
					{
					nu = atof (jalon[i + 1]);
					}
										
				if (!strcmp ((const char *) jalon[i], "cn"))
					{
					cn = atof (jalon[i + 1]);
					}
				
				// - Tangential force
				//#Force_tangential(tan_force_model,mu_gg,mu_gw)
				if (!strcmp ((const char *) jalon[i], "Force_tangential"))
					{
					forcetangential=true;
					tan_force_model= atoi (jalon[i + 1]);
					mu_gg= atof (jalon[i + 2]);
					mu_gw= atof (jalon[i + 3]);
					}
				
				// - Rolling resistant torque
				//#Torque_rolling_resistance(mu_roll_gg,mu_roll_gw)
				if (!strcmp ((const char *) jalon[i], "Torque_rolling_resistance"))
					{
					torqueRollRes=true;
					mu_roll_gg = atof (jalon[i + 1]);
					mu_roll_gw = atof (jalon[i + 2]);
					}
				

				// MODELE HARD CORE - SOFT SHELL
				// 
				if (!strcmp ((const char *) jalon[i], "HardCore_SoftShell"))
						{
						forceHCSS=true;
						paste_yield_stress= atof (jalon[i + 1]);
						paste_consistancy= atof (jalon[i + 2]);
						paste_exponentn=  atof (jalon[i + 3]);
						paste_density=  atof (jalon[i + 4]);
						paste_hmed= atof (jalon[i + 5]);
						hmax=atof (jalon[i + 6]);
						paste_kn=atof (jalon[i + 7]);
						}


				//FORCES DE VAN DER WAALS			  
				// #Force_vanderWalls(model,Hamaker_constant,distmin,distmax) // Hamaker -> J - Distance -> m
				//Force_vanderWalls(1,6.5E-20,0.4E-9,100.0E-9)
				if (!strcmp ((const char *) jalon[i], "Force_vanderWalls"))
					{
					forceVdw=true;
					modelefvdw= atof (jalon[i + 1]);
					Hamaker= atof (jalon[i + 2]);
					hmin= atof (jalon[i + 3]);
					hmax=atof (jalon[i + 4]);
					}

				//FORCES ELECTROSTATICS								
				//#Force_electrostatic(zetapotential,kappa) // Zeta potential in V, kappa in m-1
				//Force_electrostatic(0.03,1.0428e+09)
				if (!strcmp ((const char *) jalon[i], "Force_electrostatic"))
					{
					forceElectro=true;
					pzeta = atof (jalon[i + 1]);
					kappa= atof (jalon[i + 2]);
					}

		  
				//Fluid-particle interactions
				//#Forces_fluid(density_fluid,fluid_viscosity) // density_fluid -> kg/m^3 & dynamic viscosity Pa.s -> N.s/unity2
				//Forces_fluid(1.E+03,1.E-3)
				if (!strcmp ((const char *) jalon[i], "Forces_fluid"))
					{
					forceFluid=true;
					densiteFluide=atof (jalon[i + 1]);
					viscofluide=atof (jalon[i + 2]);
					}
				
		  
				
				
				// Section 3 :  TYPE OF PROBLEM,  NUMERICAL PARAMETERS and printing results
				//=========================================
				//ChoixTypeCalcul==1
				// Free systeme
				//sedimentation
				if (!strcmp ((const char *) jalon[i], "sedimentation"))
					{
					ChoixTypeCalcul=1;
					fprintf(stdout,"  Simulation case : sedimentation\n");
					fprintf(flogfile,"  Simulation case : sedimentation\n");
					} 
					
				//ChoixTypeCalcul==2
				// Moving the upper and bottom plan - shear system
				//******************
				// fluid_shear_y(shear_rate)
				if (!strcmp ((const char *) jalon[i], "fluid_shear_y"))
					{
					shear_rate= atof (jalon[i + 1]);
					ChoixTypeCalcul=2;
					fprintf(stdout,"  fluid_shear_y(shear_rate) \n");
					fprintf(flogfile,"  fluid_shear_y(shear_rate) \n");
					fprintf(flogfile,"  shear_rate %e \n",shear_rate);
					}
				
				
                // rheometry(vbase,npalier,niterpalier)
				if (!strcmp ((const char *) jalon[i], "rheometry"))
					{
					vbase=atof (jalon[i + 1]);
					npalier=atoi (jalon[i + 2]);
					niterpalier=atoi (jalon[i + 3]);
					fprintf(stdout,"  rheometry(vbase,npalier,niterpalier)\n");
					fprintf(flogfile,"  rheometry(vbase,npalier,niterpalier)\n");
					fprintf(flogfile,"  vbase %e npalier %d niterpalier %d\n",vbase,npalier,niterpalier);
					if(!forceHCSS){fprintf(stdout,"  Warning : configuration only available with the HCSS model forcefor the moment\n");}
					ChoixTypeCalcul=2;
					}
				
					
				//ChoixTypeCalcul==3 -- Moving the upper plane with other case (triaxial, stress-controled compression, cylindric condition)
				//*******************
				// compression_uni(vplateau,z_upwall_ini)
				if (!strcmp ((const char *) jalon[i], "compression_uni"))
					{
					cs_vplateau=atof (jalon[i + 1]);
					z_upwall_ini=atof (jalon[i + 2]);
					ChoixTypeCalcul=3;
					fprintf(stdout,"  compression_uni(vplateau,z_upwall_ini)\n");
					fprintf(flogfile,"  compression_uni(vplateau,z_upwall_ini)\n");
					fprintf(flogfile,"  cs_vplateau %e z_upwall_ini %e \n",cs_vplateau,z_upwall_ini);
					}
				
				
				// compression_uni_stress_cont(vplateau,cs_sr_upwall,z_upwall_ini)
				if (!strcmp ((const char *) jalon[i], "compression_uni_stress_cont"))
					{
					cs_vplateau=atof (jalon[i + 1]);
					cs_sr_upwall=atof (jalon[i + 2]);
					z_upwall_ini=atof (jalon[i + 3]);
					Fupplane_conf=true;
					ChoixTypeCalcul=3;
					fprintf(stdout,"  Simulation case : compression_uni(cs_vplateau,cs_sr_upwall,z_upwall_ini)\n");
					fprintf(flogfile,"  Simulation case : compression_uni(cs_vplateau,cs_sr_upwall,z_upwall_ini)\n");
					fprintf(flogfile,"  cs_vplateau %e cs_sr_upwall %e z_upwall_ini %e\n",cs_vplateau,cs_sr_upwall,z_upwall_ini);
					}
				
				
				
				// compression_cyl_iso(vplateau,cs_sr_latcyl,z_upwall_ini)
				if (!strcmp ((const char *) jalon[i], "compression_cyl_iso"))
					{
					cs_vplateau=atof (jalon[i + 1]);
					cs_sr_latcyl=atof (jalon[i + 2]);
					z_upwall_ini=atof (jalon[i + 3]);
					ChoixTypeCalcul=3;
					Fcyl_mov=true;
					fprintf(stdout,"  compression_cyl_iso(cs_vplateau,cs_sr_latcyl,z_upwall_ini)\n");
					fprintf(flogfile,"  compression_cyl_iso(cs_vplateau,cs_sr_latcyl,z_upwall_ini)\n");
					fprintf(flogfile,"  cs_vplateau %e cs_sr_latcyl %e z_upwall_ini %e\n",cs_vplateau,cs_sr_latcyl,z_upwall_ini);
					}	
				
				
					
				// triaxial_test(cs_vplateau,cs_sr_latcyl,z_upwall_ini)
				if (!strcmp ((const char *) jalon[i], "triaxial_test"))
					{
					cs_vplateau=atof (jalon[i + 1]);
					cs_sr_latcyl=atof (jalon[i + 2]);
					z_upwall_ini=atof (jalon[i + 3]);
					ChoixTypeCalcul=3;
					Fcyl_mov=true;
					Ftriaxial=true;
					fprintf(stdout,"  triaxial_test(cs_vplateau,cs_sr_latcyl,z_upwall_ini)\n");
					fprintf(flogfile,"  triaxial_test(cs_vplateau,cs_sr_latcyl,z_upwall_ini)\n");
					fprintf(flogfile,"  cs_vplateau %e cs_sr_latcyl %e z_upwall_ini %e \n",cs_vplateau,cs_sr_latcyl,z_upwall_ini);
					}
				
				// ChoixTypeCalcul==4  -- Moving (vibration) of the bottom plane
				//*******************
				// If vibration of the bottom wall
				//#vibration_bottom_wall(ampl,freq) // ampl -> unity / freq -> Hz
				if (!strcmp ((const char *) jalon[i], "vibration_bottom_wall"))
					{
					ampl=atof (jalon[i + 1]);
					freq=atof (jalon[i + 2]);
					ChoixTypeCalcul=4;
					fprintf(stdout,"  Simulation case : vibration_bottom_wall(ampl,freq)\n");
					fprintf(flogfile,"  Simulation case : vibration_bottom_wall(ampl,freq)\n");
					fprintf(flogfile,"  ampl %f freq %f\n",ampl,freq);
					}
				
				// ChoixTypeCalcul==5 -- 
				//*********************
				// tribology(vbase,z_upwall_ini,cs_vplateau,cs_srupwall)
				if (!strcmp ((const char *) jalon[i], "tribology"))
					{
					vbase=atof (jalon[i + 1]);
					z_upwall_ini=atof (jalon[i + 2]);
					cs_vplateau=atof (jalon[i + 3]);
					cs_sr_upwall=atof (jalon[i + 4]);
					Fupplane_conf=true;
					ChoixTypeCalcul=5;
					fprintf(stdout,"  tribology(vbase,z_upwall_ini,cs_vplateau,cs_sr_upwall) \n");
					fprintf(flogfile,"  tribology(vbase,z_upwall_ini,cs_vplateau,cs_sr_upwall) \n");
					fprintf(flogfile,"  vbase %e z_upwall_ini %e cs_vplateau %e cs_sr_upwall %e \n",vbase,z_upwall_ini,cs_vplateau,cs_sr_upwall);
					}
				fflush(flogfile);
				
				
				if (!strcmp ((const char *) jalon[i], "TauxCisaillement"))
					{
					Taux= atof (jalon[i + 1]);					
					Vflumax=0.5*syssizex*Taux; // Vitesse maximale a la paroi: Vmax=0.5*distance entre parois*taux cisaillement geany 
					} 
					
					
					// Numerical parameters
					// ---------------------------
									
				if (!strcmp ((const char *) jalon[i], "niter"))
					{
					//Number of iterations
					niter = atoi (jalon[i + 1]);
					} 

				if (!strcmp ((const char *) jalon[i], "deltat"))
					{
					deltat = atof (jalon[i + 1]);
					} 
					
				// Printing result frequencies
				// ---------------------------
				
				if (!strcmp ((const char *) jalon[i], "ndodisplay"))
					{
					ndodisplay = atof (jalon[i + 1]);
					}				
		
				// Nombre de fichier de sortie vtk pour Paraview
				if (!strcmp ((const char *) jalon[i], "ndowritevtk"))
						{
						ndowritevtk = atof (jalon[i + 1]);
						}

				if (!strcmp ((const char *) jalon[i], "ndowritemicro"))
					{
					ndowritemicro= atof (jalon[i + 1]);
					} 
				
								
				if (!strcmp ((const char *) jalon[i], "ndowrite"))
					{
					ndowrite= atof (jalon[i + 1]);
					}

				// contact history function
				if (!strcmp ((const char *) jalon[i], "contact_history"))
					{
					j=0;
					do	{
						j++;
						}
					while(jalon[j+1] != NULL);
					contact_history= calloc(j, sizeof(int));
					for (k1=1;k1<j;k1++)
						{
						contact_history[k1]=atoi (jalon[k1]);
						}
					}
				//FIN DE LA BOUCLE
				++i;
				}
			}
		fclose (fichierPARAM);
		}
	else
		{
		fprintf(stderr,"Warning! No parameter file %s in the directory \n",Nparafile);
		exit(ierror=10);
		}
		
	// Initialization of the parameters (variables and flags)
	//** Number of surfaces: 20 by default (isurf==1-6, isurf==7 : cylinder along z, isurf==8 ....)
	
	nobj=nsurf;

	for(k1=1;k1<=6;k1++)
		{
		if(state_obj[k1]==1)
			{
			Fwall=true;
			}
		}

	if(state_obj[7]==1)
		{
		// Set flag for the cylindrical boundary conditions
		Fcyl_bond=true;
		}

	if(state_obj[nsurf+1]==1)
		{
		// Read the number of objects
		nobj+=add_obj_spheres(Nobjfile,0,true);
		Fsph_obj=true;
		}
		
	if(!forceFluid)
		{
		densiteFluide=0.0;
		viscofluide = 0.0;
		}
	if(forceFluid&&(ChoixTypeCalcul==2))
	    {
		Vflumax=0.0;Vflumax=0.5*syssizex*shear_rate; // Vitesse maximale a la paroi: Vmax=0.5*distance entre parois*taux cisaillement
		}
	else if((!forceFluid)&&(ChoixTypeCalcul==2)){fprintf(stdout,"  Warning : no fluid force - shear rate  %f",shear_rate);}


	if(init_gravity)
		{
		gravity.x=0.0;
		gravity.y=0.0;
		gravity.z=-9.81;  // m.s-2
		}
					
	if(!forceContact)
		{
		nu=0.0;
		cn=0.0;
		}
	if((!forcetangential)||(!forceContact))
		{
		mu_gg=0.0;
		mu_gw=0.0;
		}
	if((!torqueRollRes)||(!forceContact))
		{
		mu_roll_gg=0.0;
		mu_roll_gw=0.0;
		}
			
	if(ChoixTypeCalcul!=4)
		{
		ampl=0.0;
		freq=0.0;
		}
	
	if((ChoixTypeCalcul==2)&&(forceVdw))
		{
		// Enable the computation of cluster properties
		Fclust=true;
		}

	// Initialization of coll array (used by choc function)
	for(k1=0;k1<MAXCONT;k1++)
		{
		coll[k1]=0;
		}
		
// initialisation of structure 

	vect0.x=0.0;
	vect0.y=0.0;
	vect0.z=0.0;
	

	// length Conversion of the parameter data
	kn=kn*unity*unity; // Young Modulus Pa -> N/unity²
	hmin=hmin/unity;  // m -> unity
	hmax=hmax/unity;  // m-> unity
	paste_hmed=paste_hmed/unity; // m -> unity
	Hamaker=Hamaker/unity;
	kappa=kappa*unity;
	densiteFluide=densiteFluide*pow(unity,3.0); // kg/m3 -> kg/unity3 				
	viscofluide=viscofluide*pow(unity,2.0); // Pa.s -> N.s/unity2
	paste_kn=paste_kn*unity; //N/m -> N/unity
	paste_yield_stress=paste_yield_stress*unity*unity; // Pa -> N/unity2
	paste_consistancy=paste_consistancy*unity*unity; // Pa.s -> N.s/unity2
	freq=2*PI*freq; // Hz -> rad.s-1
	hcone = hcone / unity; // m -> unity
	rbcone = rbcone / unity; // m -> unity
	rtcone= rtcone / unity; // m -> unity
		
	//Writing some information in the logfile

	fprintf(flogfile,"\n  Parameter file: %s \n",Nparafile );
	fprintf(flogfile,"=============================================\n");
	fprintf (flogfile,"  nseed %d - unity %e (m)\n",*seed,unity);
	fprintf (flogfile,"  System size %d %d %d \n",syssizex,syssizey,syssizez);
	fprintf (flogfile,"  npart %d nsurf %d nobj %d \n", npart,nsurf,nobj);
	fprintf (flogfile,"  Density of the material %lf (Kg/m3) \n", density);
	fprintf (flogfile,"\n  Contact force (1:enabled/0:disable): %d - Tangential force model %d \n",forceContact,tan_force_model);
	fprintf (flogfile,"  Young modulus %e (Pa) - Poisson coefficient %lf \n",kn/(unity*unity),nu);
	fprintf (flogfile,"  Friction coefficient: %f (grain/grain) %f (grain/wall)\n", mu_gg,mu_gw);
	fprintf(flogfile,"  Rolling friction (0:disable) %d - Rolling friction coefficient (grain/grain) %f (grain/wall) %f \n",torqueRollRes,mu_roll_gg,mu_roll_gw);
	fprintf(flogfile,"  Normal damping coefficient %lf (sec) \n",cn);

	fprintf (flogfile,"\n  Van der Waals force (1:enable/0:disable): %d \n",forceVdw);
	if(forceVdw){fprintf(flogfile,"Model %d - Hamaker constant %e (J)- hmin %e (unity) hmax %e (unity)\n",modelefvdw,Hamaker*unity,hmin,hmax);}
	fprintf (flogfile,"\n  Electrostatic force (1:enable/0:disable): %d \n",forceElectro);
	if(forceElectro){fprintf(flogfile,"Zeta potential %e (Volt) - kappa %e (m-1)- epsi0 %e (F.m-1 / m-3 kg-1 s4 A2) epsir %e (-)\n",pzeta,kappa/unity,epsi_vide,epsi_rel);}
	fprintf (flogfile,"\n  Hydrodynamic force (1:enable/0:disable): %d \n",forceFluid);
	fprintf (flogfile,"\n  Gravity %lf %lf %lf (m.s-1) \n",gravity.x,gravity.y,gravity.z);
	fprintf(flogfile,"\n  Time step (s): %lf - Number of iterations: %d - Total time : %lf (s)\n",deltat,niter,deltat*((double)niter));
	}
	*/
