//
// Created by ThibLMT on 30/03/2021.
//

#include "write_micro.cuh"
#include "string.h"

void microfile_write(const char *Nfile,discrete_elt *particle,geom_struct *geom)
{
    FILE *fmicrofile;
    FILE *fmicrocontfile;
    char Ncontactfile[256];
    int j;
    unsigned int nb_conttot,i;

    // Write Coordinate of particle (micro file)
    fmicrofile=fopen(Nfile,"w");
    if (fmicrofile != nullptr)
    {
        // writing of particle number and system sizes on the first line
        fprintf(fmicrofile,"\t%d\t%d\t%d\t%d\t%lf\n",geom->nb_part,geom->sizex,geom->sizey,geom->sizez,geom->unity);
        // writing of particle number and system sizes on the first line
        for(i=1;i<=geom->nb_part;i++)
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
    for(i=1;i<=geom->nb_part;i++)
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
    if (fmicrocontfile != nullptr)
    {
        // Writing firstly the number of contacts
        fprintf(fmicrocontfile,"\t%d\n",nb_conttot);
        //

        for(i=1;i<=geom->nb_part;i++)
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