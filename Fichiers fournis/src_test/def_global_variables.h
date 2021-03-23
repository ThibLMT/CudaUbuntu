/**
*\file def_global_variables.h
*\brief definition of global variables
*/
 
#ifndef __DEF_GLOBAL_VARIABLES_H
#define __DEF_GLOBAL_VARIABLES_H

#include "def_types.h"
#include <omp.h>
	
// Definitions of global variables
// General variables to run the program  
// Declaration of the general files
FILE *flogfile; //!< Log file
FILE *fhistcontfile; //!< File containing some contact information between particles
FILE *fresultfile;
unsigned int ierror;  	//!< Integer giving the error value in case of problem or bug (equal to 0 by default or errno if standard C error) 

// Name of the general files
char FichMicroInitiale[50];  			 	// nom fichier contenant la microstructure initiale agregats
char FichContInitiale[50];   			    // nom fichier contenant les contacts initiaux agregats
char *Nparafile;							// name of the parameter file (set to params.in or with the following command: demGCE -file)
char Nobjfile[50];


// geometry of system


// 
discrete_elt *particle;  //!< array containing the whole description of particles and forces
material_data prop_mat_part;
geom_struct geom;

//
double unity; 
unsigned int npart;		//!< Number of particles
unsigned int nsurf;		//!<Number of surfaces
unsigned int nobj;		//!< Number of objects (like wall, cylinder, etc,...)
int syssizex;		//!< size of the background grid along x axis (unity)	
int syssizey;       //!< size of the background grid along y axis (unity)	
int syssizez;       //!< size of the background grid along z axis (unity)
int syssizel;		//!< Maximum number of particles for each voxel of the background grid
int ncont;		//!< Maximum number of particles per voxel. Ncont is used to set the 4th colomn of mic array.
int state_obj[100];    //!< state_obj[nobj] state of the object 0: no object 1 : object
double rpart; //!< Radius of particles given by parameter file 
int nseed,*seed;	//!<global random number seed


// Arrays defining the particle assembly and the digitalization by voxels
unsigned int ****mic;   //!< Pointer to 4D mic array  describing the whole digitized system. Mic [][][][n] containes the number of the nth particle connected to every voxel
unsigned int ****mic_boundary; //!< store state of mic array right after boundary_conditions routine call
unsigned int ***mic_insert, ***mic_insert_boundary;
//struct sphere *particle;  //!< array containing the whole description of particles and forces
unsigned int *tab_mic;
unsigned int *tab_mic_insert;

unsigned int *backgrid;
int *backgrid_insert;


// Utiliser pour récupérer les identifiants, structure pas adéquat à changer
long int coll[MAXCONT]; //!< 1D array containing ID's of particles in contact with the particle checked with the check function 

 



// definition of vectors
vector vect0;     //!< Vanished vector
vector box_size;    //!< Size of box simulation
vector gravity;    //!< Gravitational acceleration vector

// movement of boundary conditions and computation condition
// *********************************************************
int generationmodel;   //!<  Model for initial system generation
Flag Freadmicro;   //!< Read initial microstructure
int ChoixTypeCalcul;      //!<  Type of simulation

double cs_vplateau;   //!< Desired value for the velocity of the superior plane
double z_upwall_ini; //!< Initial position of the upper plane
double vbase;	//!< Velocity of lower plane (z=0)
double diamcyl;  //!< Diameter of the cylindrical boundary conditions
double cs_sr_latcyl; //!< Required stress value (consigne) for the cylinder
double cs_sr_upwall; //!< Required stress value (consigne) for the upper plane
double hcone,rbcone,rtcone;    //!< height and radii of the conical boundary conditions
unsigned int npalier,niterpalier;  //!< Define step 
double freq; //!< angular frequency (vibration of bottom plane)
double ampl; //!< amplitude (vibration of bottom plane)

Flag Fsph_obj;  //!< Enable sphere objects
Flag Fcyl_bond;  //!< Enable cylinder boundary condition
Flag Fcyl_mov;  //!< The cylinder boundary condition moves during the simulation
Flag Fwall;   //!< Enable wall boundary condition
Flag Fwall_mov;   //!< The wall boundary condition moves during the simulation
Flag Fupplane_conf; //!< The upper wall moves during the simulation as  a function to the required stress (servo stress-controlled)
Flag Fclust;   //!< Enable cluster
Flag Ftriaxial; //!< Enable triaxial conditions

// Material parameters for the interaction laws
// ********************************************
// - Flags
Flag forceContact;  //!< contact force laws
Flag forceVdw;	//!< van der Walls interactions
Flag forceElectro; //!< electrostatic interactions
Flag forceFluid;  //!< particle-fluid interactions	
Flag forceHCSS; //!< HCSS model 
Flag torqueRollRes; //!< Rolling resistance torque

// Material parameters
double kn; //!< Young modulus
double nu; //!< Poisson coefficient
double cn; //!< Normal dashpot constant (sse walton model for the dissipative part of the normal contact force)
double density; //!< Bulk density

// - Parameters of interaction law
double hmax,hmin; //!< Layer thickness of the "shell" around grain, debend of the selected interaction law (used in, van der Walls, electrostatic and HCSS model
// -- Tangential law

int tan_force_model;	//!< Model for the tangential force (1 or 2)
double mu_gg; //!< Friction coefficient grain-grain
double mu_gw; //!< Friction coefficient grain-wall

// -- Rolling resistant torque
double mu_roll_gg; //!< Rolling friction coefficient grain-grain
double mu_roll_gw; //!< Rolling friction coefficient grain-wall


// -- Hard Core Soft Shell model - Paste parameters
double paste_density; //!< Paste density
double paste_yield_stress, paste_consistancy, paste_exponentn; //!< caractéristiques rhéologiques de Herschel-Bulkley attribuées aux SS pour le calcul des forces tangentielles
double paste_kn; //!< raideur normale pour les interactions SS/SS
double paste_hmed; //!< profondeur de recouvrement à partir de laquelle la force normale SS/SS (de raideur paste_kn) est activée

// -- van der Walls and electrostatic forces
int modelefvdw; //!< Model of van der Walls forces (1 or 2)
double Hamaker; //!< Hamaker constant (van der Walls parameters)
double kappa,pzeta; //!< electrostatic force parameters

// -- Fluid forces
double Vflumax; //!< Fluid velocity at the extremun defenition of the shear plane flowing
double densiteFluide; //!< Fluid density
double viscofluide; //!< Fluid viscosity

// Numerical parameters
// ********************
double cond_stab;	//!< Not used for teh moment, Stability condition to obtain a small portion of critical time step (see Cundall et al., 1979, Geotechnique, vol. 29, 47-65.
unsigned int niter;		//!< Number of iterations
unsigned int istep;      //!< Step (iteration) i
//double deltat;      //!< Time step in second

// Recording and displaying paremeters
//************************************
unsigned int ndodisplay; //!< Iteration increment for screen display
unsigned int ndowritemicro; //!< Iteration increment for microstructure and contact files writing
unsigned int ndowritevtk; //!< Iteration increment for vtk files writing
unsigned int ndowrite; //!< Iteration increment for result file

// Packing analyse values
// **********************
double overlap[5]; //!<Overlap mean value, min value, max value, n1 for max value,n2 for max value
unsigned int nconttot; //!<Number of the contacts (inside granular media + with boundary conditions) at each time step
unsigned int ncontgran; //!<Number of the contacts (inside granular media) at each time step
unsigned int ncontslide; //!<Number of the sliding contacts at each time step
double indfric; //!< Friction index
unsigned int *contact_history;  //!< Pointer to 1D array containing the number of particles to follow (see contact_history function)

	
#endif



