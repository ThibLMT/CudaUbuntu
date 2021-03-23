#include "def_global_variables.h"
#include "cal_compacity.h"

double compacity_z(double zmin, double zmax,double sxy)
	{
	unsigned int k1;
	double dist,density,vols;			
	density=0.0;
	for(k1=1;k1<=npart;k1++)
		{
		vols=4.0/3.0*PI*particle[k1].radius*particle[k1].radius*particle[k1].radius;
		if((particle[k1].Ri.z+particle[k1].radius<=zmax)&&(particle[k1].Ri.z-particle[k1].radius>=zmin))
			{
			density+=vols;
			}
		else if((particle[k1].Ri.z<zmax)&&(particle[k1].Ri.z+particle[k1].radius>zmax))
			{
			dist=(particle[k1].Ri.z+particle[k1].radius)-zmax;
			density+=vols-PI*dist*dist*(particle[k1].radius-dist/3.0);
			}
		else if((particle[k1].Ri.z>zmax)&&(particle[k1].Ri.z-particle[k1].radius<zmax))
			{
			dist=zmax-(particle[k1].Ri.z-particle[k1].radius);
			density+=PI*dist*dist*(particle[k1].radius-dist/3.0);
			}
		else if((particle[k1].Ri.z<zmin)&&(particle[k1].Ri.z+particle[k1].radius>zmin))
			{
			dist=(particle[k1].Ri.z+particle[k1].radius)-zmin;
			density+=PI*dist*dist*(particle[k1].radius-dist/3.0);
			}
		else if((particle[k1].Ri.z>zmin)&&(particle[k1].Ri.z-particle[k1].radius<zmin))
			{
			dist=zmin-(particle[k1].Ri.z-particle[k1].radius);
			density+=vols-PI*dist*dist*(particle[k1].radius-dist/3.0);
			}
		}

	density=density/(double)(sxy*(zmax-zmin));
	return(density);
	// Note : Volume of a portion of sphere : pi*H*H*(R-H/3) with R radius and H height of the portion
	}

//*******************
//** compacity_profil_z

void compacity_profil_z(const char *Nprofilfile,double sxy)
	{
	unsigned int k1,k3;
	double lmin,h,comps,compb,rs,rb,hempmax,dz,r;	
	FILE *fprofil;
	
	// Writing of the microstructure file
	fprofil=fopen(Nprofilfile,"w");
	fprintf(fprofil,"%-12s\t%-12s\t%-12s\n","#Lmin","Comps","Compb");
	lmin=0.0;
	rs=particle[1].radius;
	hempmax=0.0;
	h=0.0;
	r=0.0;
	
	for(k1=1;k1<=npart;k1++)
		{
		// calculating hmax
		if((particle[k1].Ri.z+particle[k1].radius)>hempmax)
			{
			hempmax=particle[k1].Ri.z+particle[k1].radius;
			}
		// rs (radius small) and rb (big radius) 
		
		if(particle[k1].radius<rs)
			{
			rs=particle[k1].radius;
			}
		else
			{
			rb=particle[k1].radius;
			}
		}
		
	dz=0.1*2*rs;
	
	do
		{
		comps=0.0;
		compb=0.0;
		for(k3=1;k3<=npart;k3++)
		
			{
		
			if((particle[k3].Ri.z+particle[k3].radius>lmin)&&(particle[k3].Ri.z-particle[k3].radius<lmin))
		{	
			h=particle[k3].Ri.z-lmin;
				if(particle[k3].radius==rs)
				{
					r=sqrt(pow(rs,2)-pow(h,2));
					comps+=PI*pow(r,2);
				}
				else if(particle[k3].radius==rb)
				{
					r=sqrt(pow(rb,2)-pow(h,2));
					compb+=PI*pow(r,2);
				}
			
		}
			
		}
		comps=comps/sxy;
		compb=compb/sxy;
		fprintf(fprofil,"%lf\t%lf\t%lf\n",lmin,comps,compb);	
		lmin+=dz;
	
	} while (lmin<=hempmax);
	fclose(fprofil);
	fprintf(flogfile,"    Writing %s rs %f rb %f \n",Nprofilfile,rs,rb);
	}




