/*
 * def_types.h
 *
 *  Created on: 14 déc. 2012
 *      Author: patrick
 */

#ifndef DEF_TYPES_H_
#define DEF_TYPES_H_

#include "def_const.h"
#include <stdbool.h>


// Type definitions
//*************************

typedef bool Flag;

// Declaration of the structures
//*************************

typedef struct vect      //!< coordinate of vector in 3D
	{
	double x;								//!< following x axe
	double y;								//!< following y axe
	double z;								//!< following y axe
	} vector;

// geometry of system
typedef struct geom_sys {
	unsigned int nb_part; //!< Number of particles
	unsigned int nb_bc ; //!< Number of boundary conditions
	unsigned int max_cont_per_part; //!< Maximum number of contact per particle (equal to MAXCONT )
	double unity; //!< length Conversion rate 
	int sizex; //!< Size along x axis for the background grid
	int sizey; //!< Size along y axis for the background grid
	int sizez; //!< Size along z axis for the background grid
	int sizel; //!< Size of 4th dimension of  backgrid
	double deltat;
	} geom_struct;

// Structure used to define the discrete elemnt
typedef struct sphere     //!< data structure used for the description of spherical particles in 3D
	{
	vector Ri;						//!< position of the center of the particle
	vector Vi;						//!< linear velocity of the particle */
	/**/vector Ai;					//!< Accelerations lineaires
	vector Wi;						//!< angular velocity of the particle
	/**/vector Aroti; 				//!< Accelarations angulaires
	vector Fi;                      //!< total force exerted on the particle
	vector Mi;						//!< total momentum exerted on the particle
	double radius;					//!< radius of the particle 
	double mass;					//!< mass of the particle
	double inertia;					//!< inertia of the particle
	double Yn;						//!< Module de Young
	double Nu;						//!< coefficient de poisson
	double Ndamp;					//!<  coefficient d'amortissement normal particule-particule ou particule-paroi
	double Mu;						//!< coefficient de frottement particule-particule ou particule-paroi
	double Mur;						//!< rolling resistance friction coefficient
	int next;						//!< ID of next particle for the definition of clusters (-1 if end of cluster) 
	int clust;                      //!< Definit l'appartenance d'une particule à un cluster de particules
	int type[MAXCONT];			//!< contact particles  - Id des autre particules en contact 
	int contact[MAXCONT];			//!< contact particles  - Id des autre particules en contact 
	double ut[MAXCONT];				//!< cumulated tangential displacement at contact point (scalar)
	vector ftanold[MAXCONT];		//!< Tangential force  at contact point at the previous time step
	vector nijold[MAXCONT];			//!< Normal direction  at contact point at the previous time step
	} discrete_elt;
	
typedef struct mat     //!< Material parameter
	{
	// Material parameters
	double E; //!< Young modulus
	double nu; //!< Poisson coefficient
	double cn; //!< Normal dashpot constant (sse walton model for the dissipative part of the normal contact force)
	double density; //!< Bulk density

	// -- Tangential law
	int tan_force_model; //!< Model for the tangential force (1 or 2)
	double mu_gg; //!< Friction coefficient grain-grain
	double mu_gw; //!< Friction coefficient grain-wall
	
	// -- Rolling resistant torque
	double mu_roll_gg; //!< Rolling friction coefficient grain-grain
	double mu_roll_gw; //!< Rolling friction coefficient grain-wall
	} material_data;



#endif /* DEF_TYPES_H_ */
