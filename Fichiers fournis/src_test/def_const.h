/**
*\file def_const.h
*\brief definition of constants of the program
*/

#ifndef DEF_CONST_H_
#define DEF_CONST_H_

// Including of the standard libraries
//--------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>



// Definition of the macros
//*************************
 // DEM system constants

#define MAXCONT 			    20      /* maximum number of contacts per particle */

 // Universal constants
#define PI 				   		3.14159265358979		//!< valeur de PI o e
#define DTOR                  0.017453        			//!< convert degrees to radians
#define RTOD            	   57.29578        			//!< convert radians to degrees
# define planck               6.62606896E-34      		//!< constante de Planck
# define cboltzmann           1.3806504E-23         	//!< constante de Boltzmann
# define epsi_vide            8.854187E-12          	//!< permitivity of vacum C^2/J/m = (C^2/J.m)
# define epsi_rel             80.0                  	//!< relative permitivity of water [-] [M.Yang,1997]
# define Kelvin               273.15                 	//!< absolute temperature [K]
# define Avogadro_nbre        6.02214179E+23        	//!< Avogadro's number
# define charge_elec          1.602176487E-19      	    //!< charge of electron [C]
# define NmaxEspece           10                   		//!< Nombre maximal d'espece ions

// Specific macro
#define SQR(a)          ((a)*(a))      								//!< a squared

// Boolean values
#define NON 			0
#define OUI 			1
#define TRUE            1
#define FALSE           0
#define true            1
#define false           0
#define yes             1
#define no              0


#endif /* DEF_CONST_H_ */
