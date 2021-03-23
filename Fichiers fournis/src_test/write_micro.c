/**
*\file write_micro.c
*\brief body of the function write_micro
*/
#include "write_micro.h"
#include "string.h"
//#include "def_global_variables.h"

void microfile_write(const char *Nfile,discrete_elt *particle,geom_struct geom)
  {
  FILE *fmicrofile;
  FILE *fmicrocontfile;
  char Ncontactfile[256];
  int j;
  unsigned int nb_conttot,i;
  
  // Write Coordinate of particle (micro file)
  fmicrofile=fopen(Nfile,"w");
	if (fmicrofile != NULL)
		{
		// writing of particle number and system sizes on the first line
		fprintf(fmicrofile,"\t%d\t%d\t%d\t%d\t%lf\n",geom.nb_part,geom.sizex,geom.sizey,geom.sizez,geom.unity);
		// writing of particle number and system sizes on the first line
		for(i=1;i<=geom.nb_part;i++)
			{
			
			fprintf(fmicrofile,"%u\t%16.9lf\t%16.9lf\t%16.9lf\t",i,particle[i].Ri.x,particle[i].Ri.y,particle[i].Ri.z);
			fprintf(fmicrofile,"%16.9lf\t%16.9lf\t%16.9lf\t",particle[i].Vi.x,particle[i].Vi.y,particle[i].Vi.z);
			fprintf(fmicrofile,"%16.9lf\t%16.9lf\t%16.9lf\t%16.9lf\n",particle[i].Wi.x,particle[i].Wi.y,particle[i].Wi.z,particle[i].radius);
						
			}
		fclose(fmicrofile);
		}
	else
		{
		perror(Nfile);
		exit(errno);
		}
		
	// Write details of contacts between particles (micro%d_cont file)
	sprintf(Ncontactfile,"%s_cont",Nfile);
	
	// Fisrt compute the number of contact
	nb_conttot=0;
	for(i=1;i<=geom.nb_part;i++)
		{
		for(j=0;j<MAXCONT;j++)
			{
			if(particle[i].contact[j]>0.0)
				{
				nb_conttot++;
			    }
		}
	}
	// Write file
	fmicrocontfile=fopen(Ncontactfile,"w");
		if (fmicrocontfile != NULL)
			{
			// Writing firstly the number of contacts	
			fprintf(fmicrocontfile,"\t%d\n",nb_conttot);
			// 
			
			for(i=1;i<=geom.nb_part;i++)
				{
				for(j=0;j<MAXCONT;j++)
					{
					if(particle[i].contact[j]>0.0)
						{
						
						// basic Information of the contact i (n1,n3,n2,type)
						fprintf(fmicrocontfile,"%d\t%d\t%d\t%d\t",i,j,particle[i].contact[j],particle[i].type[j]);
						// Tangential cumulative displacement (scalar value)
						fprintf(fmicrocontfile,"%e\t",particle[i].ut[j]);
						// Tangential force and normal direction
						fprintf(fmicrocontfile,"%e\t%e\t%e\t",particle[i].ftanold[j].x,particle[i].ftanold[j].y,particle[i].ftanold[j].z);
						fprintf(fmicrocontfile,"%e\t%e\t%e\n",particle[i].nijold[j].x,particle[i].nijold[j].y,particle[i].nijold[j].z);
						
						}
					}
				
				}
			fclose(fmicrocontfile);
			}
		else
			{
			perror(Ncontactfile);					
			exit(errno);
			}	
		
		}
  
/*

void write_micro(const char *Nmicrofile,_Bool write_contact)
	{
	// write_contact : writing condition 1 -> Writing of the contact and micro files -> 0 only micro file
	unsigned int i,j,ncontact=0;
	FILE *fmicro,*fcontacts;
	char Ncontactfile[256];
	
	// Writing of the microstructure file
	fmicro=fopen(Nmicrofile,"w");
	if (fmicro != NULL)
		{
		// writing of particle number and system sizes on the first line
		fprintf(fmicro,"\t%d\t%d\t%d\t%d\t%lf\n",npart,syssizex,syssizey,syssizez,unity);
		// writing of particle number and system sizes on the first line
		for(i=0;i<=npart;i++)
			{
			if(particle[i].radius>0.0)
				{
				fprintf(fmicro,"%u\t%16.9lf\t%16.9lf\t%16.9lf\t",i,particle[i].Ri.x,particle[i].Ri.y,particle[i].Ri.z);
				fprintf(fmicro,"%16.9lf\t%16.9lf\t%16.9lf\t",particle[i].Vi.x,particle[i].Vi.y,particle[i].Vi.z);
				fprintf(fmicro,"%16.9lf\t%16.9lf\t%16.9lf\t%16.9lf\n",particle[i].Wi.x,particle[i].Wi.y,particle[i].Wi.z,particle[i].radius);
				}				
			}
		fclose(fmicro);
		}
	else
		{
		perror(Nmicrofile);
		exit(ierror=errno);
		}
	// Writing of the contact file
	
	if(write_contact)
		{
		sprintf(Ncontactfile,"%s_cont",Nmicrofile);
		
		fcontacts=fopen(Ncontactfile,"w");
		
		if (fcontacts != NULL)
			{
			// Writing firstly the number of contacts	
			fprintf(fcontacts,"\t%d\n",nconttot);
			// 
			
			for(i=1;i<=npart;i++)
				{
				for(j=0;j<MAXCONT;j++)
					{
					if(particle[i].contact[j]>0.0)
						{
						ncontact++;
						// basic Information of the contact i (n1,n3,n2,type)
						fprintf(fcontacts,"%d\t%d\t%d\t%d\t",i,j,particle[i].contact[j],particle[i].type[j]);
						// Tangential cumulative displacement (scalar value)
						fprintf(fcontacts,"%e\t",particle[i].ut[j]);
						// Tangential force and normal force at the previous time step
						fprintf(fcontacts,"%e\t%e\t%e\t",particle[i].ftanold[j].x,particle[i].ftanold[j].y,particle[i].ftanold[j].z);
						fprintf(fcontacts,"%e\t%e\t%e\n",particle[i].nijold[j].x,particle[i].nijold[j].y,particle[i].nijold[j].z);
						
						}
					}
				
				}
			fclose(fcontacts);
			}
		else
			{
			perror(Ncontactfile);					
			exit(ierror=errno);
			}	
		
		}
	}

*/
