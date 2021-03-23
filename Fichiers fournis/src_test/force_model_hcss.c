/**
*\file rolling_resistance_torque.c
*\brief body of the function tangential_force_mindlin
*\version $Revision: 1.4 $  $Date: 2013/06/13 14:17:32 $
*/

#include <stdio.h>
#include <stdlib.h>

#include "force_model_hcss.h"
#include "def_global_variables.h"

double normal_force_hcss(double deltan)
	{
	//
	double fnpaste;

	if(deltan+paste_hmed>0)
		{
		fnpaste=-paste_kn*(deltan+paste_hmed);
		}
	else
		{
		fnpaste=0.0;
		}

	return fnpaste;
	}



vector tangential_force_hcss(double deltan,vector nij,double r1,double r1eq,vector v1,vector w1,double r2,double r2eq,vector v2,vector w2)
	{
	double ep; //!<
	double rs,surface;
	double r1eq2,r2eq2;
	double deltanpdot;
	double normet;
	double dist;
	double gap;
	double norme_ftpaste;
	vector vijrelt,vijrel,ftpaste;

	if(deltan+hmax>=0)
		{
		dist=r1+r2-deltan;
		ep=hmax/2.0;
		r1eq=(dist*dist+(r1+ep)*(r1+ep)-(r2+ep)*(r2+ep))/(2.0*dist);
		r2eq=dist-r1eq;

		// Computation of the contact radius between
		rs=(r1+ep)*(r1+ep)-r1eq*r1eq;
		surface=rs*PI;
		rs=sqrt(rs);
		//
		if(r1eq>=r1){r1eq2=r1eq;}
		else{r1eq2=r1;}
		if(r2eq>=r2){r2eq2=r2eq;}
		else{r2eq2=r2;}


		vijrel.x=v1.x*(r1/r1eq2)*(r1/r1eq2)-v2.x*(r2/r2eq2)*(r2/r2eq2)+((w1.y*(r1/r1eq2)*(r1/r1eq2)*(r1eq2)+w2.y*(r2/r2eq2)*(r2/r2eq2)*(r2eq2))*nij.z-(w1.z*(r1/r1eq2)*(r1/r1eq2)*(r1eq2)+w2.z*(r2/r2eq2)*(r2/r2eq2)*(r2eq2))*nij.y);
		vijrel.y=v1.y*(r1/r1eq2)*(r1/r1eq2)-v2.y*(r2/r2eq2)*(r2/r2eq2)+((w1.z*(r1/r1eq2)*(r1/r1eq2)*(r1eq2)+w2.z*(r2/r2eq2)*(r2/r2eq2)*(r2eq2))*nij.x-(w1.x*(r1/r1eq2)*(r1/r1eq2)*(r1eq2)+w2.x*(r2/r2eq2)*(r2/r2eq2)*(r2eq2))*nij.z);
		vijrel.z=v1.z*(r1/r1eq2)*(r1/r1eq2)-v2.z*(r2/r2eq2)*(r2/r2eq2)+((w1.x*(r1/r1eq2)*(r1/r1eq2)*(r1eq2)+w2.x*(r2/r2eq2)*(r2/r2eq2)*(r2eq2))*nij.y-(w1.y*(r1/r1eq2)*(r1/r1eq2)*(r1eq2)+w2.y*(r2/r2eq2)*(r2/r2eq2)*(r2eq2))*nij.x);

		deltanpdot=vijrel.x*nij.x+vijrel.y*nij.y+vijrel.z*nij.z;

		vijrelt.x=vijrel.x-deltanpdot*nij.x;
		vijrelt.y=vijrel.y-deltanpdot*nij.y;
		vijrelt.z=vijrel.z-deltanpdot*nij.z;
		normet=sqrt(vijrelt.x*vijrelt.x+vijrelt.y*vijrelt.y+vijrelt.z*vijrelt.z);

		if(normet>0.0)
			{
			gap=-deltan;
			if(gap<paste_hmed){gap=paste_hmed;}
			
			norme_ftpaste=-surface*(paste_yield_stress+paste_consistancy*pow(2.0*normet/(gap),paste_exponentn));

			ftpaste.x=norme_ftpaste*vijrelt.x/normet;
			ftpaste.y=norme_ftpaste*vijrelt.y/normet;
			ftpaste.z=norme_ftpaste*vijrelt.z/normet;
			}
		else{ftpaste.x=0.0;ftpaste.y=0.0;ftpaste.z=0.0;}
		}
	return ftpaste;
	}

