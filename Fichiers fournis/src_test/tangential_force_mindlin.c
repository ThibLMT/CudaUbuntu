/**
*\file tangential_force_mindlin.c
*\brief body of the function tangential_force_mindlin
*/
#include "def_global_variables.h"
#include "tangential_force_mindlin.h"


double tangential_force_mindlin(double deltanl,double deltatanl,double fnijl,double mul, double Eijl, double Gijl)
	{
	double ftij;
	double deltatanmax;
	double indfricij;

	// The tangential force is given by the equation of Mindlin and Deriewsic in the case where the normal force is constant
	// with the increasing of the tangential force. deltatamax is the value of maximal tangential displacement
	// corresponding to total sliding at the contact (Coulomb criteria is activated)

	deltatanmax=mul*Eijl*deltanl/(4*Gijl);


	deltatanl=fmin(deltatanl,deltatanmax);

	// Checking the status of contact (sliding or not)
	if(deltatanl==deltatanmax){ncontslide++;}
	// Computing the tangential force

	ftij=-1.0*mul*fnijl*(1.0-pow((1.0-deltatanl/deltatanmax),3.0/2.0));
	indfricij = (1.0-pow((1.0-deltatanl/deltatanmax),3.0/2.0));
	indfric= indfric + indfricij;

	return(ftij);
	}

