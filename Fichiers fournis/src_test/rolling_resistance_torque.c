/**
*\file rolling_resistance_torque.c
*\brief body of the function tangential_force_mindlin
*/

#include <stdio.h>
#include <stdlib.h>

#include "rolling_resistance_torque.h"

vector rolling_resistance_torque(double fnijl,double murol, double Rij, vector resis_tor,vector wi,vector wj)
	{
	double normewij;
	vector wij;

	// The expression of the resistance torque is given by Balevicius and al (Powder technology 2011, num 22, pp226-235)
	// It is not a true rolling resistance torque because the result vector is along wij, so rolling and twisting resistance torque.

	wij.x=wi.x-wj.x;
	wij.y=wi.y-wj.y;
	wij.z=wi.z-wj.z;

	normewij=sqrt(wij.x*wij.x+wij.y*wij.y+wij.z*wij.z);

	if(normewij==0)
		{
		resis_tor.x=0.0;
		resis_tor.y=0.0;
		resis_tor.z=0.0;
		}
	else
		{
		wij.x=wij.x/normewij;
		wij.y=wij.y/normewij;
		wij.z=wij.z/normewij;

		resis_tor.x=murol*fnijl*Rij*wij.x;
		resis_tor.y=murol*fnijl*Rij*wij.y;
		resis_tor.z=murol*fnijl*Rij*wij.z;
		}

	return resis_tor;
	}

