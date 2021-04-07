//
// Created by ThibLMT on 23/03/2021.
//

#include "read_micro.cuh"

unsigned int microfile_read_npart(const char *Nfile)
{
    FILE *fmicrofile;
    unsigned int npart;
    fmicrofile=fopen(Nfile,"r");
    if (fmicrofile != NULL)
    {
        (void)fscanf(fmicrofile,"%d",&npart);
    }
    else
    {
        perror(Nfile);
        exit(errno);
    }
    return npart;
}


void microfile_read_particle(const char *Nfile,discrete_elt *particle,geom_struct *geom)
{
    FILE *fmicrofile;
    unsigned int ipart;
    int i;
    int npart,syssizex,syssizey,syssizez;
    double unity;



    fmicrofile=fopen(Nfile,"r");
    if (fmicrofile != NULL)
    {
        // Read the first line which contain information about number of particle, syssize of system and unity.
        (void)fscanf(fmicrofile,"%d",&npart);
        (void)fscanf(fmicrofile,"%d",&syssizex);
        (void)fscanf(fmicrofile,"%d",&syssizey);
        (void)fscanf(fmicrofile,"%d",&syssizez);
        (void)fscanf(fmicrofile,"%lf",&unity);
        (void)fscanf(fmicrofile,"\n");

        geom->nb_part = npart;

        // Scanning of the coordinates of particles and placement in the mic array
        for(i=1;i<=npart;i++)
        {
            (void)fscanf(fmicrofile,"%d",&ipart);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Ri.x);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Ri.y);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Ri.z);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Vi.x);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Vi.y);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Vi.z);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Wi.x);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Wi.y);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].Wi.z);
            (void)fscanf(fmicrofile,"%lf",&particle[ipart].radius);
            (void)fscanf(fmicrofile,"\n");
        }
        fclose(fmicrofile);
    }
    else
    {
        perror(Nfile);
        exit(errno);
    }

    fprintf(stdout,"\n  Read microstructure file : %s \n",Nfile);
    fprintf(stdout,"=============================================\n");
    fprintf(stdout,"  npart %d syssizes %d %d %d unity %lf\n",npart,syssizex,syssizey,syssizez,unity);
}




void microcontfile_read_contact(const char *Nfile,discrete_elt *particle,geom_struct *geom)
{
    FILE *fmicrocontfile;
    char Ncontactfile[256];
    unsigned int i;
    int nconttot,ipart,scan_kij,kij;

    // Fist add extension
    sprintf(Ncontactfile,"%s_cont",Nfile);

    // Open file
    fmicrocontfile=fopen(Ncontactfile,"r");
    if (fmicrocontfile != NULL)
    {
        fscanf(fmicrocontfile,"%d",&nconttot);
        if(nconttot==0){fprintf(stdout,"Warning! No contact - ncontot = %d \n",nconttot);}
        for(i=0;i<=nconttot;i++)
        {
            fflush(stdin);
            fscanf(fmicrocontfile,"%d",&ipart);

            fscanf(fmicrocontfile,"%d",&scan_kij);
            kij=scan_kij;
            fscanf(fmicrocontfile,"%d",&particle[ipart].contact[kij]);
            fscanf(fmicrocontfile,"%d",&particle[ipart].type[kij]);
            fscanf(fmicrocontfile,"%lf",&particle[ipart].ut[kij]);
            fscanf(fmicrocontfile,"%lf %lf %lf",&particle[ipart].ftanold[kij].x,&particle[ipart].ftanold[kij].y,&particle[ipart].ftanold[kij].z);
            fscanf(fmicrocontfile,"%lf %lf %lf",&particle[ipart].nijold[kij].x,&particle[ipart].nijold[kij].y,&particle[ipart].nijold[kij].z);
        }
        fclose(fmicrocontfile);
    }
    else
    {
        perror(Ncontactfile);
        exit(errno);
    }

}