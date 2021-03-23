/**
*\file normal_force_hertz.c
*\brief body of the function normal_force_hertz
*/
#include "normal_force_hertz.h"



double normal_force_hertz(double deltan,double deltandot, double rij, double Eij, double Aij)
	{
	double fcnij;

	// The interaction force is given by two parts. The First one concerns the elastic contribution given by Hertz model
	// The second part is dissipative force allows taking into account the energy dissipation at contact during a collision.
	// The model used here is described in T. Pöschel and T. Schwager, Computational Granular dynamics Springer book
	 //(equation 2.14 p20)


	fcnij=(-4.0/3.0)*sqrt(rij)*Eij*(pow(deltan,3.0/2.0)+Aij*deltandot*sqrt(deltan));


	// If the force is positive due to an attractive contribution of dissipative force during the contact separation,
	// the force must vanish to avoid an attractive effect (detail in Computational Granular Dynamics book, p. 22).

	if (fcnij>0.0) {fcnij=0;}

	return(fcnij);
	}

