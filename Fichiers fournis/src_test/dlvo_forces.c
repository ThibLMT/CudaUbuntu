/**
*\file van_der_waals_force.c
*\brief body of the function van_der_waals_force and functions 
*/

#include "dlvo_forces.h"

extern double Hamaker;
extern double hmin,pzeta,kappa;


double van_der_waals_force(double hij,double ri,double rj,double rij,int model_vdw)
	{
	double numfvdw,den1fvdw,den2fvdw;
	double attracvdw;

	attracvdw=0.0;
	switch(model_vdw)
		{
		case 1 :
			//======================================================
			//Model 1 :
			//J. N. Israelachvili, Intermolecular and Surface Forces, 2nd ed. (Academic, London,1991)
			//H. C. Hamaker, Physica (Amsterdam), 4, 1058 (1937)
			//======================================================
			numfvdw=64.0*pow(ri,3.0)*pow(rj,3.0)*(hij+ri+rj);
			den1fvdw=pow(hij,2.0)+2.0*ri*hij+2.0*rj*hij;
			den2fvdw=den1fvdw+4.0*ri*rj;
			attracvdw=(1.0/6.0)*(Hamaker*numfvdw/(pow(den1fvdw,2.0)*pow(den2fvdw,2.0)));
			break;

		case 2 :
			//==================================================
			//Model 2: M.L. Eggersdorfer and Herrmann et al. 342(2010) 261-268
			//=====================================================
			attracvdw=(Hamaker*rij)/(6.0*pow(hij,2.0));
			break;
			}

	// sign convention
	// The unit vector points from grain n1 toward grain n2. We compute the force excerced by grain n2 on Grain n1.
	// The van der Wall forces are attractives, the scalar product n.F must be positive.

	if(attracvdw<0.0){attracvdw=0.0;}
	return(attracvdw);
	}
	
// Electrostatic forces
double NormalElectrostaticForce(double hij,double ri, double rj,double rij)
	{
	double felectro,prefact,term1,term2;
	// ref ?

	if(hij<hmin){hij=hmin;}
	if((ri*kappa>1.0)&&(rj*kappa>1.0))
		{
		if((hij<ri)&&(hij<rj))
			{

			prefact=4.0*PI*epsi_rel*epsi_vide*rij*pzeta*pzeta;
			term1=kappa*exp(-kappa*hij);
			term2=1.0+exp(-kappa*hij);
			felectro=-prefact*(term1/term2);
			}
		}
	else{felectro=0.0;}
	//The electrostatic force is repulsive force. This normal component must be negative.
	if(felectro>0.0){felectro=0.0;}
	return(felectro);
	}


