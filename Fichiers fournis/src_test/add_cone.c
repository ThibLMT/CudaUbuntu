/**
*\file add_cone.c
*\brief body of the function add_cone
*/


#include "def_global_variables.h"
#include "add_cone.h"

void add_cone(int nwall)
	{
	int i,j,k,l,l1;
	double rk,dist0,dist1,dist2,im,ip,jm,jp;
	//struct sphere part_wall;
	
	// Scan all the pixels
	for(i=0;i<syssizex;i++)
		{
		for(j=0;j<syssizey;j++)
			{
				for(k=0;(k<syssizez)&&(k<=hcone);k++)
					{
					// For a given height k, compute the radius of cone at that height
					rk=(hcone-(double)k)/hcone*(rbcone-rtcone)+rtcone-hmax/2;
					// Compute horizontal distance between center of pixel and axix of cone
					dist0=sqrt((i-syssizex/2.0)*(i-syssizex/2.0)+(j-syssizey/2.0)*(j-syssizey/2.0));
					
					//remove cylindre condition inside the cone
					if(dist0<rk+1.0+hmax/2.0)
						{
						//printf(stdout,"mic0=%d mic1=%d mic2=%d\n",
						l=0;
						while(mic[i][j][k][l]!=0)
							{
							if(mic[i][j][k][l]==npart+7)
								{
								mic[i][j][k][l]=0;
								l1=l+1;
								while(mic[i][j][k][l1]!=0)
									{
									mic[i][j][k][l1-1]=mic[i][j][k][l1];
									l1++;
									}
								}
							l++;
							}
						}

					// if pixel is able to be cut by cone
					if((dist0>rk-1.0)&&(dist0<rk+1.0))
						{
							// Compute radius of cone for the height k-0.5
							rk=(hcone-(double)k+0.5)/hcone*(rbcone-rtcone)+rtcone;
							// Compute mini and maxi distances from axis of cone to edges of pixel on the face located in k-0.5
							im=i;ip=i;
							jm=j;jp=j;
							if(i>syssizex/2.0){im=i-0.5;ip=i+0.5;}
							else if(i<syssizex/2.0){im=i+0.5;ip=i-0.5;}
							if(j>syssizey/2.0){jm=j-0.5;jp=j+0.5;}
							else if(j<syssizey/2.0){jm=j+0.5;jp=j-0.5;}
							dist1=sqrt((im-syssizex/2.0)*(im-syssizex/2.0)+(jm-syssizey/2.0)*(jm-syssizey/2.0));
							dist2=sqrt((ip-syssizex/2.0)*(ip-syssizex/2.0)+(jp-syssizey/2.0)*(jp-syssizey/2.0));
							if((dist1<=rk)&&(dist2>=rk))
								{
									l=0;
									while(mic[i][j][k][l]!=0) {l++;}
									mic[i][j][k][l]=nwall;

									if(k>0)
										{
										l=0;
										while(mic[i][j][k][l]!=0) {l++;}
										mic[i][j][k-1][l]=nwall;
										}
								}
						}
					}
			}
		}
		
	// Geometrie 
	particle[nwall].radius=1000000.0;

	particle[nwall].Ri=vect0;
	particle[nwall].Vi=vect0;
	particle[nwall].Wi=vect0;

	particle[nwall].Fi=vect0;
	particle[nwall].Mi=vect0;
	


	fprintf(stdout,"cone id : %d - height : %lf - lower radius : %lf - upper radius : %lf \n",nwall,hcone,rbcone,rtcone);	
	
	fprintf(flogfile," =======================================\n");
	fprintf(flogfile,"cone id : %d - height : %lf - lower radius : %lf - upper radius : %lf \n",nwall,hcone,rbcone,rtcone);	

	} 


