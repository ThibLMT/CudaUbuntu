/**
*\file fluid_forces.c
*\brief body of the fluid_force and fluid_torque functions
*/
#include "fluid_forces.h"
extern double viscofluide,densiteFluide;
extern vector gravity;
extern double unity;


vector fluid_force(double radius,vector vpart,vector vflu)
	{
	vector fluid_force;
	
	// Buoyancy force
	fluid_force.x=-(4.0/3.0)*PI*pow(radius,3.0)*densiteFluide*gravity.x;
	fluid_force.y=-(4.0/3.0)*PI*pow(radius,3.0)*densiteFluide*gravity.y;
	fluid_force.z=-(4.0/3.0)*PI*pow(radius,3.0)*densiteFluide*gravity.z;
	
	// Stokes force
	fluid_force.x+=-6.0*PI*radius*viscofluide*(vpart.x-vflu.x);
	fluid_force.y+=-6.0*PI*radius*viscofluide*(vpart.y-vflu.y);
	fluid_force.z+=-6.0*PI*radius*viscofluide*(vpart.z-vflu.z);
	//printf("ftot %e ftot %e ftot %e \n",fluid_force.x,fluid_force.y,fluid_force.z);
	return fluid_force;
	}
	
vector fluid_torque(double radius,vector wpart,vector vorflu)
	{
	vector fluid_torque;
	fluid_torque.x=-8.0*PI*viscofluide*unity*pow(radius,3.0)*(wpart.x-vorflu.x);
	fluid_torque.y=-8.0*PI*viscofluide*unity*pow(radius,3.0)*(wpart.y-vorflu.y);
	fluid_torque.z=-8.0*PI*viscofluide*unity*pow(radius,3.0)*(wpart.z-vorflu.z);
	//printf("fluid_torque %e ftot %e ftot %e \n",fluid_torque.x,fluid_torque.y,fluid_torque.z);
	return 	fluid_torque;
	}

/*
double HydrodynamicDragForces(int idparticule)
	{
	double Vflux,Vfluy,Vfluz;
	double vx,vy,vz;
	double Fdragx,Fdragy,Fdragz;
	double Parchimedx,Parchimedy,Parchimedz;
					
	Parchimedx=0.0;Parchimedy=0.0;Parchimedz=0.0;
	Fdragx=0.0;Fdragy=0.0;Fdragz=0.0;
					
	if(FLUIDE)
		{			
		// Calcul de la force de Poussee
		Parchimedx=-(4.0/3.0)*PI*pow(particle[idparticule].radius,3.0)*densiteFluide*gravity.x;
		Parchimedy=-(4.0/3.0)*PI*pow(particle[idparticule].radius,3.0)*densiteFluide*gravity.y;
		Parchimedz=-(4.0/3.0)*PI*pow(particle[idparticule].radius,3.0)*densiteFluide*gravity.z;
						 
		particle[idparticule].Fi.x+=Parchimedx;
		particle[idparticule].Fi.y+=Parchimedy;
		particle[idparticule].Fi.z+=Parchimedz;
					
		// Vecteur vitesse solide : vitesse de la particule
		vx=particle[idparticule].Vi.x;
		vy=particle[idparticule].Vi.y;
		vz=particle[idparticule].Vi.z;
					
					
		if(ChoixTypeCalcul==2){
						 
			
						 
			
			// Vecteur vitesse fluide : ecoulement suivant l'axe y 
			Vflux=0.0;
			Vfluy=Vflumax*((2.0/syssizex)*particle[idparticule].Ri.x-1.0);
			Vfluz=0.0;
									
			// Calcul de la force de train�e : depend du signe de la vitesse relative 
				entre la particule et le fluide
				Fdrag>0 si vi-vflu>0; 
				Fdrag<0 si vi-vflu<0;
				Fdrag=0 si vi-vflu=0;
			//
									
			Fdragx=-6.0*PI*particle[idparticule].radius*viscofluide*(vx-Vflux);
			Fdragy=-6.0*PI*particle[idparticule].radius*viscofluide*(vy-Vfluy);
			Fdragz=-6.0*PI*particle[idparticule].radius*viscofluide*(vz-Vfluz);
											
						
				}
			else{	
				// Forces de Stokes : ici vitesse du fluide nulle 
				Fdragx=-6.0*PI*particle[idparticule].radius*viscofluide*vx;
				Fdragy=-6.0*PI*particle[idparticule].radius*viscofluide*vy;
				Fdragz=-6.0*PI*particle[idparticule].radius*viscofluide*vz;		
				}
						
			//================================== 
			 //  Contribution de la force de trainee 
			//=====================================
									
			particle[idparticule].Fi.x+=Fdragx;
			particle[idparticule].Fi.y+=Fdragy;
			particle[idparticule].Fi.z+=Fdragz;
			
			} // Fin de la condition si fluide
	return 1.0;				
	}


void HydrodynamicTorque(int idparticule)
	{
	double rotflux, rotfluy,rotfluz;
	double wx,wy,wz,Mdragx,Mdragy,Mdragz;

	Mdragx=0.0;
	Mdragy=0.0;
	Mdragz=0.0;
					
	if(FLUIDE){
					
		// Vecteur vitesse rotation particule			
		wx=particle[idparticule].Wi.x;
		wy=particle[idparticule].Wi.y;
		wz=particle[idparticule].Wi.z;
					
		if(ChoixTypeCalcul==2)
			{			
			//Vecteur vitesse rotation fluide : vorticite omaga=1/2 gradient vectoriel vitesse fluide
					   
			rotflux=0.0;
			rotfluy=0.0;
			rotfluz=Vflumax/syssizex;
					   
			// Calcul du moment de trainee			
			Mdragx=-8.0*PI*viscofluide*unity*pow(particle[idparticule].radius,3.0)*(wx-rotflux);
			Mdragy=-8.0*PI*viscofluide*unity*pow(particle[idparticule].radius,3.0)*(wy-rotfluy);
			Mdragz=-8.0*PI*viscofluide*unity*pow(particle[idparticule].radius,3.0)*(wz-rotfluz);				
			}
		else
			{				
			Mdragx=-8.0*PI*viscofluide*unity*pow(particle[idparticule].radius,3.0)*wx;
			Mdragy=-8.0*PI*viscofluide*unity*pow(particle[idparticule].radius,3.0)*wy;
			Mdragz=-8.0*PI*viscofluide*unity*pow(particle[idparticule].radius,3.0)*wz;			
			}
		} 
								  
	particle[idparticule].Mi.x=Mdragx;
	particle[idparticule].Mi.y=Mdragy;
	particle[idparticule].Mi.z=Mdragz;
	}
	*/

