/**
*\file write_cyl_vtk.c
*\brief body of the function write_cyl_vtk
*/

#include "def_global_variables.h"
#include "write_cyl_vtk.h"

void write_cyl_vtk(const char *Nvtkfile)
	{
	int i,resolu,tab[5];
	float m;
	struct vecteur
		{
		double x;
		double y;
		double z;
		};
	struct vecteur point[200];
	resolu=50; 			//resolution du cylindre
	FILE *fparaview;
	fparaview = fopen(Nvtkfile,"w");

	if (fparaview != NULL) 
		{
		// Writing tyPIcal header for a vtk file
		fprintf(fparaview,"# vtk DataFile Version 3.4.0 \n");
		fprintf(fparaview,"Cylindre file\n");
		fprintf(fparaview,"ASCII\n");
		fprintf(fparaview,"DATASET POLYDATA\n");
		fprintf(fparaview,"POINTS %d float\n",resolu*4);

		// Calcul des coordonn�es des points dans le structure
		for (i=0;i<2*resolu-1;i=i+2)
			{
			point[i].x=cos(2*PI-PI/resolu*i)*diamcyl/2+0.5*box_size.x;
			point[i].y=sin(2*PI-PI/resolu*i)*diamcyl/2+0.5*box_size.y;
			point[i].z=box_size.z;
			}
		for (i=1;i<2*resolu;i=i+2)
			{
			point[i].x=cos(2*PI-PI/resolu*i+PI/resolu)*diamcyl/2+0.5*box_size.x;
			point[i].y=sin(2*PI-PI/resolu*i+PI/resolu)*diamcyl/2+0.5*box_size.y;
			point[i].z=0;
			}

		//Imprimer les coordonn�es de chaque point
		for (i=0;i<2*resolu;i=i+1)
			{
			fprintf(fparaview,"%f %f %f\n ",point[i].x,point[i].y,point[i].z);
			}
		//Imprimer deux facettes
		for (i=0;i<2*resolu-1;i=i+2)
			{
			fprintf(fparaview,"%f %f %f\n ",point[i].x,point[i].y,point[i].z);
			}
		for (i=1;i<2*resolu;i=i+2)
			{
			fprintf(fparaview,"%f %f %f\n ",point[i].x,point[i].y,point[i].z);
			}

		fprintf(fparaview,"\n");
		fprintf(fparaview,"POLYGONS");
		fprintf(fparaview," %d %d\n", resolu+2,7*resolu+2);

		tab[0]=4;
		tab[1]=0;
		tab[2]=1;
		tab[3]=3;
		tab[4]=2;
		for (i=0;i<resolu-1;i++)
			{
			fprintf(fparaview,"%d %d %d %d %d\n ",tab[0],tab[1],tab[2],tab[3],tab[4]);
			tab[0]=4;
			tab[1]=tab[1]+2;
			tab[2]=tab[2]+2;
			tab[3]=tab[3]+2;
			tab[4]=tab[4]+2;
			}
		
		fprintf(fparaview,"%d %d %d %s %s ",tab[0],tab[1],tab[2],"1","0\n");
		fprintf(fparaview,"%d ",resolu);
		
		for (i=1;i<resolu+1;i=i+1)
			{
			fprintf(fparaview,"%d ",2*resolu-1+i);
			}
		
		fprintf(fparaview,"\n");
		fprintf(fparaview,"%d ",resolu);
		
		for (i=1;i<resolu+1;i=i+1)
			{
			fprintf(fparaview,"%d ",3*resolu-1+i);
			}
		
		fprintf(fparaview,"\n");
		
		//Imprimer normals

		fprintf(fparaview,"POINT_DATA %d \n",resolu*4); 
		fprintf(fparaview,"NORMALS Normals float\n");
		for (i=0;i<2*resolu;i=i+1)
			{
			fprintf(fparaview,"%f %f %d\n ",(point[i].x-0.5*box_size.x)/(diamcyl/2),(point[i].y-0.5*box_size.y)/(diamcyl/2),0);
			}
		for (i=0;i<resolu;i=i+1)
			{
			fprintf(fparaview,"0 0 1\n");
			}
		for (i=0;i<resolu;i=i+1)
			{
			fprintf(fparaview,"0 0 -1\n");
			}
		fprintf(fparaview,"\n");

		//Imprimer texture coordinate
		fprintf(fparaview,"TEXTURE_COORDINATES TCoords 2 float\n");
		for (i=resolu;i>-1;i=i-2)
			{
			m=(float)i/(float)resolu;
			fprintf(fparaview,"%f %s\n ",m,"0");
			fprintf(fparaview,"%f %s\n ",m,"1");
			}
		for (i=2;i<resolu;i=i+2)
			{
			m=(float)i/(float)resolu;
			fprintf(fparaview,"%f %s\n ",m,"0");
			fprintf(fparaview,"%f %s\n ",m,"1");
			}
		for (i=0;i<resolu;i=i+1)
			{
			fprintf(fparaview,"%f %f\n",cos(2*PI-2*PI/resolu*i),sin(2*PI-2*PI/resolu*i));
			}
		for (i=resolu-1;i>-1;i=i-1)
			{
			fprintf(fparaview,"%f %f\n",cos(2*PI-2*PI/resolu*i),sin(2*PI-2*PI/resolu*i));
			}

		fprintf(fparaview,"\n");
		}
	}


