//
// Created by ThibLMT on 24/03/2021.
//

#include "sub_domain.cuh"
#include <stdio.h>
#include <math.h>

__global__ void initialize_backgrid(unsigned int *backgrid,int *backgrid_insert,geom_struct *geom)
{
    int size_backgrid,size_backgrid_insert;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    size_backgrid_insert=geom->sizex*geom->sizey*geom->sizez;
    size_backgrid=size_backgrid_insert*geom->sizel;

    for (int i = index; i < size_backgrid_insert; i+= stride)
    {
        backgrid_insert[i] = 0;
    }

    for (int i = index; i < size_backgrid; i+= stride)
    {
        backgrid[i] = 0;
    }
}

__device__ unsigned int get_id_backgrid(int xv,int yv,int zv,int lv,unsigned int *backgrid,geom_struct *geom)
{
    int index_backgrid;
    index_backgrid=lv+zv*geom->sizel+yv*geom->sizel*geom->sizez+xv*geom->sizel*geom->sizez*geom->sizey;
    return backgrid[index_backgrid];
}

__host__ void set_id_backgrid(int xv,int yv,int zv,unsigned int idp,unsigned int *backgrid, int *backgrid_insert,geom_struct *geom)
{
    int index_backgrid,index_backgrid_insert,lv;
    // Keep the last position lv
    index_backgrid_insert=zv+yv*geom->sizez+xv*geom->sizez*geom->sizey;
    lv=backgrid_insert[index_backgrid_insert];
    index_backgrid=lv+zv*geom->sizel+yv*geom->sizel*geom->sizez+xv*geom->sizel*geom->sizez*geom->sizey;
    backgrid[index_backgrid]=idp;
    backgrid_insert[index_backgrid_insert]++;
}

__global__ void insert_sph_backgrid(discrete_elt *particle,unsigned int *backgrid,int *backgrid_insert,geom_struct *geom)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    // xpart,ypart,zpart --> sphere center position (double)
    // radpart --> radius of the sphre
    // idpart --> id number of the sphere
    int kp, index_insert, index_backgrid;
    int cont, xp, yp, zp, lp, i, j, k, ip1, jp1, kp1, ip2, jp2, kp2;
    double d, dist, xdist, ydist, zdist, x1, y1, z1, x2, y2, z2;
    int syssizex, syssizey, syssizez, syssizel;
    syssizex = geom->sizex;
    syssizey = geom->sizey;
    syssizez = geom->sizez;
    syssizel = geom->sizel;

    lp = 0;
    cont = 0;
    for (int idpart = index; idpart <= geom->nb_part; idpart+= stride) {
        double xpart = particle[idpart].Ri.x;
        double ypart = particle[idpart].Ri.y;
        double zpart = particle[idpart].Ri.z;
        double radpart = particle[idpart].radius;
        /* identification of the pixels that have to be scanned */
        /* xp,yp,zp are located in the middle of the pixel */
        /* if the double position is greater than xp.5 the pixel */

        /* bearing that point is located at xp+1 */

        ip1 = (int) (xpart - radpart);
        ip2 = (int) (xpart + radpart);

        if (xpart - radpart - (double) ip1 > 0.5) { ip1 += 1; }
        else if (xpart - radpart - (double) ip1 < -0.5) { ip1 -= 1; }

        if (xpart + radpart - (double) ip2 > 0.5) { ip2 += 1; }

        jp1 = (int) (ypart - radpart);
        jp2 = (int) (ypart + radpart);

        if (ypart - radpart - (double) jp1 > 0.5) { jp1 += 1; }
        else if (ypart - radpart - (double) jp1 < -0.5) { jp1 -= 1; }

        if (ypart + radpart - (double) jp2 > 0.5) { jp2 += 1; }

        kp1 = (int) (zpart - radpart);
        kp2 = (int) (zpart + radpart);

        if (zpart - radpart - (double) kp1 > 0.5) { kp1 += 1; }
        else if (zpart - radpart - (double) kp1 < -0.5) { kp1 -= 1; }

        if (zpart + radpart - (double) kp2 > 0.5) { kp2 += 1; }

        /* Check all pixels within the digitized cube volume */
        for (i = ip1; i <= ip2; i++) {
            x2 = xpart;
            xp = i;

            // use periodic boundary conditions for sphere placement
            if (xp < 0) {
                xp += syssizex;
                x2 += syssizex;
            }
            else if (xp > syssizex - 1) {
                xp -= syssizex;
                x2 -= syssizex;
            }

            x1 = (int) xp;

            d = x2 - (double) xp;
            if ((d <= 0.5) && (d >= -0.5)) { x1 = x2; }
            else if (x2 > xp) { x1 = xp + 0.5; }
            else if (x2 < xp) { x1 = xp - 0.5; }
            xdist = (x2 - x1);

            for (j = jp1; j <= jp2; j++) {
                y2 = ypart;
                yp = j;

                // use periodic boundary conditions for sphere placement
                if (yp < 0) {
                    yp += syssizey;
                    y2 += syssizey;
                }
                else if (yp > syssizey - 1) {
                    yp -= syssizey;
                    y2 -= syssizey;
                }

                y1 = (int) yp;

                d = y2 - (double) yp;
                if ((d <= 0.5) && (d >= -0.5)) { y1 = y2; }
                else if (y2 > yp) { y1 = yp + 0.5; }
                else if (y2 < yp) { y1 = yp - 0.5; }
                ydist = (y2 - y1);

                for (k = kp1; k <= kp2; k++) {
                    z2 = zpart;
                    zp = k;

                    /* use periodic boundary conditions for sphere placement */
                    if (zp < 0) {
                        zp += syssizez;
                        z2 += syssizez;
                    }
                    else if (zp > syssizez - 1) {
                        zp -= syssizez;
                        z2 -= syssizez;
                    }

                    /* compute minimum distance from center of particle to pixel */
                    z1 = (int) zp;

                    d = z2 - (double) zp;

                    if ((d <= 0.5) && (d >= -0.5)) { z1 = z2; }
                    else if (z2 > zp) { z1 = zp + 0.5; }
                    else if (z2 < zp) { z1 = zp - 0.5; }

                    zdist = (z2 - z1);
                    dist = sqrt(xdist * xdist + ydist * ydist + zdist * zdist);


                    // check distance between voxel and the center of particle
                    if (dist < radpart) {

                        index_insert = zp + yp * syssizez + xp * syssizez * syssizey;
                        // Perte d'identifiants de particules au bout de 100000 itérations  race data sur le tableau backgrid_insert

                        lp=++backgrid_insert[index_insert];
                        lp--;

                        index_backgrid = lp + zp * syssizel + yp * syssizel * syssizez + xp * syssizel * syssizez * syssizey;
                        backgrid[index_backgrid] = idpart;




                        // TODO try to make fprintf work in the GPU
                        /*if ((lp+1)==syssizel) {
                            fprintf(stdout,"Error during the inscription of the id sphere in mic array \n");
                            fprintf(stdout,"  lp > syssizel lp : %d syssizel: %d \n",lp,syssizel);
                            for(lp=0;lp<syssizel;lp++) {
                                fprintf(stdout,"%u ",get_id_backgrid(xp,yp,zp,lp,backgrid,geom));
                            }
                            fprintf(stdout,"\n"),
                                    exit(EXIT_FAILURE);
                        }*/
                    } // end if dist<radpart

                } // end for k
            } // end for j
        } // end for i
    }

    // return coordination number (number of contacts per particle)
    /*return(cont);*/
}

__device__ int detect_contact_sph_backgrid(discrete_elt *particle, double rad, geom_struct *geom,unsigned int idpart,unsigned int* list_part,unsigned int *backgrid)
{
    unsigned int idparttest,npart;
    int ncontact,ncontactpot,icontact,icontactpot;
    int size_list_part_pot;
    int cont,e,xp,yp,zp,i,j,k,k1,k4,ip1,jp1,kp1,ip2,jp2,kp2;
    double d,dist,xdist,ydist,zdist,x1,y1,z1,r1,x2,y2,z2;
    double nor_cyl,nx_cyl,ny_cyl,nz_cyl;
    double xin,yin,zin;
    int syssizex,syssizey,syssizez;
    Flag new_pot_contact;


    double diamcyl,rbcone,rtcone,hcone,hmax;
    hmax=0.0;

    unsigned int list_part_pot[500];

    size_list_part_pot=500;
    memset(list_part_pot,0,size_list_part_pot*sizeof(unsigned int));

    xin=particle[idpart].Ri.x;
    yin=particle[idpart].Ri.y;
    zin=particle[idpart].Ri.z;
    syssizex=geom->sizex;
    syssizey=geom->sizey;
    syssizez=geom->sizez;
    npart=geom->nb_part;

    /* identification of the pixels that have to be scanned */
    /* xp,yp,zp are located in the middle of the pixel */
    /* if the double position is greater than xp.5 the pixel */

    /* bearing that point is located at xp+1 */

    ip1=(int)(xin-rad);
    ip2=(int)(xin+rad);

    if(xin-rad-(double)ip1>0.5){ip1+=1;}
    else if(xin-rad-(double)ip1<-0.5){ip1-=1;}

    if(xin+rad-(double)ip2>0.5){ip2+=1;}

    jp1=(int)(yin-rad);
    jp2=(int)(yin+rad);

    if(yin-rad-(double)jp1>0.5){jp1+=1;}
    else if(yin-rad-(double)jp1<-0.5){jp1-=1;}

    if(yin+rad-(double)jp2>0.5){jp2+=1;}

    kp1=(int)(zin-rad);
    kp2=(int)(zin+rad);

    if(zin-rad-(double)kp1>0.5){kp1+=1;}
    else if(zin-rad-(double)kp1<-0.5){kp1-=1;}

    if(zin+rad-(double)kp2>0.5){kp2+=1;}



    /* Check all voxel within the digitized cube volume and keep the id of particles in list_*/
    /*
    icontactpot=0;
    for(i=ip1;i<=ip2;i++)
      {
      xp=i;
      if(xp<0) {xp+=syssizex;}
      else if(xp>syssizex-1) {xp-=syssizex;}

      for(j=jp1;j<=jp2;j++)
      {
      yp=j;
      if(yp<0) {yp+=syssizey;}
      else if(yp>syssizey-1) {yp-=syssizey;}
       for(k=kp1;k<=kp2;k++)
          {
          zp=k;

            //use periodic boundary conditions for sphere placement
            if(zp<0) {zp+=syssizez;}
            else if(zp>syssizez-1) {zp-=syssizez;}
          k1=0;
                while(get_idp_backgrid(xp,yp,zp,k1,backgrid)!=0)
              {

              list_part_pot[icontactpot]=get_idp_backgrid(xp,yp,zp,k1,backgrid);
              icontactpot++;
               k1++;

              }
          }

      }
      }
      */
// Boucle sur les voxles avec détection des distances
    icontactpot=0;
    for(i=ip1;i<=ip2;i++)
    {
        x2=xin;
        xp=i;

        // use periodic boundary conditions for sphere placement
        if(xp<0) {xp+=syssizex;x2+=syssizex;}
        else if(xp>syssizex-1) {xp-=syssizex;x2-=syssizex;}
        x1=(int)xp;
        for(j=jp1;j<=jp2;j++)
        {
            y2=yin;
            yp=j;

            // use periodic boundary conditions for sphere placement
            if(yp<0) {yp+=syssizey;y2+=syssizey;}
            else if(yp>syssizey-1) {yp-=syssizey;y2-=syssizey;}
            y1=(int)yp;

            for(k=kp1;k<=kp2;k++)
            {
                z2=zin;
                zp=k;

                // use periodic boundary conditions for sphere placement
                if(zp<0) {zp+=syssizez;z2+=syssizez;}
                else if(zp>syssizez-1) {zp-=syssizez;z2-=syssizez;}
                z1=(int)zp;

                // compute minimum distance from center of particle to pixel

                d=x2-(double)xp;


                if ((d<=0.5)&&(d>=-0.5)){x1=x2;}
                else if(x2>xp){x1=xp+0.5;}
                else if (x2<xp) {x1=xp-0.5;}

                xdist=(x2-x1);
                d=y2-(double)yp;

                if ((d<=0.5)&&(d>=-0.5)){y1=y2;}
                else if(y2>yp){y1=yp+0.5;}
                else if (y2<yp) {y1=yp-0.5;}

                ydist=(y2-y1);
                d=z2-(double)zp;

                if ((d<=0.5)&&(d>=-0.5)){z1=z2;}
                else if(z2>zp){z1=zp+0.5;}
                else if (z2<zp) {z1=zp-0.5;}

                zdist=(z2-z1);
                dist=sqrt(xdist*xdist+ydist*ydist+zdist*zdist);

                // test distance
                if(dist<rad)
                {
                    // flag=3 detect contacts between idpart and other grains cut the voxel and then set the list_part array

                    k1=0;
                    while(get_id_backgrid(xp,yp,zp,k1,backgrid,geom)!=0)
                    {
                        //if(icontactpot>=size_list_part_pot){printf("problem in icontactpot %d %d\n",icontactpot,size_list_part_pot);fflush(stdout);exit(EXIT_FAILURE);}

                        list_part_pot[icontactpot]=get_id_backgrid(xp,yp,zp,k1,backgrid,geom);
                        icontactpot++;
                        k1++;

                    }
                }
            }
        }
    }

    //if(idpart==1){for(k=0;k<icontactpot;k++) {printf("%d ",list_part_pot[k]);}printf("\n %d \n",icontactpot);}


    // Sort potential contact arrays --> insertion sorting and remove the double ide
    ncontactpot=list_part_pot_sorting(list_part_pot,icontactpot);
    //if(idpart==1){for(k=0;k<icontactpot;k++) {printf("%d ",list_part_pot[k]);}printf("\n %d \n",icontactpot);}

    // Test if contact between the particle idpart and other potential particles (id of the listpart_ppt arrayBoundary conditions or particles (with boundary conditions)
    ncontact=0;
    for(icontact=0;icontact<ncontactpot;icontact++)
    {
        idparttest=list_part_pot[icontact];
        if(idparttest!=idpart){

            // Keep the position of boundaray condition of the tested particle (boundary condition or periodic condition)
            //  Walls
            if((idparttest>npart)&&(idparttest<=npart+6))
            {
                x1=particle[idparttest].Ri.x;
                y1=particle[idparttest].Ri.y;
                z1=particle[idparttest].Ri.z;
                r1=particle[idparttest].radius;
                xdist=(xin-x1);
                ydist=(yin-y1);
                zdist=(zin-z1);

            }
                // Cylinder
            else if(idparttest==npart+7)
            {
                nx_cyl=xin-syssizex/2.0;
                ny_cyl=yin-syssizey/2.0;
                nor_cyl=sqrt(nx_cyl*nx_cyl+ny_cyl*ny_cyl);
                nx_cyl=nx_cyl/nor_cyl;
                ny_cyl=ny_cyl/nor_cyl;
                x1=syssizex/2.0+nx_cyl*(particle[idparttest].radius+diamcyl/2.0);
                y1=syssizey/2.0+ny_cyl*(particle[idparttest].radius+diamcyl/2.0);
                z1=z2;
                r1=particle[idparttest].radius;
                xdist=(xin-x1);
                ydist=(yin-y1);
                zdist=(zin-z1);
            }
                // Cone
            else if(idparttest==npart+8)
            {
                r1=particle[idparttest].radius;
                z1=z2-(rbcone-rtcone)/hcone*sqrt((x2-syssizex/2.0)*(x2-syssizex/2.0)+(y2-syssizey/2.0)*(y2-syssizey/2.0));
                nx_cyl=x2-syssizex/2.0;
                ny_cyl=y2-syssizey/2.0;
                nz_cyl=z2-z1;
                nor_cyl=sqrt(nx_cyl*nx_cyl+ny_cyl*ny_cyl+nz_cyl*nz_cyl);
                nx_cyl=nx_cyl/nor_cyl;
                ny_cyl=ny_cyl/nor_cyl;
                nz_cyl=nz_cyl/nor_cyl;
                dist=(rbcone/(rbcone-rtcone)*hcone-z1)*sin(atan((rbcone-rtcone)/hcone));
                dist=dist+r1;
                x1=syssizex/2.0+dist*nx_cyl;
                y1=syssizey/2.0+dist*ny_cyl;
                z1=z1+dist*nz_cyl;
                xdist=(xin-x1);
                ydist=(yin-y1);
                zdist=(zin-z1);
            }
                // spherical particles
            else
            {
                x1=particle[idparttest].Ri.x;
                y1=particle[idparttest].Ri.y;
                z1=particle[idparttest].Ri.z;
                r1=particle[idparttest].radius;
                /*
                // periodic boundary conditions for x2
                if((x1>x2)&&(x1-x2>x2+syssizex-x1)){x1-=syssizex;}
                else if((x1<x2)&&(x2-x1>x1+syssizex-x2)){x1+=syssizex;}

                // periodic boundary conditions for y2
                //if((y1>y2)&&(y1-y2>y2+syssizey-y1)){y1-=syssizey;}
                //else if((y1<y2)&&(y2-y1>y1+syssizey-y2)){y1+=syssizey;}

                // periodic boundary conditions for z2
                //if((z1>z2)&&(z1-z2>z2+syssizez-z1)){z1-=syssizez;}
                //else if((z1<z2)&&(z2-z1>z1+syssizez-z2)){z1+=syssizez;}
                */

                xdist=(xin-x1);
                ydist=(yin-y1);
                zdist=(zin-z1);
                // periodic boundary conditions
                if(fabs(xdist)>(((double)syssizex)/2.0)) {xdist=fabs(xdist)-syssizex;}
                if(fabs(ydist)>(((double)syssizey)/2.0)) {ydist=fabs(ydist)-syssizey;}
                if(fabs(zdist)>(((double)syssizez)/2.0)) {zdist=fabs(zdist)-syssizez;}
            }

            // distance
            dist=sqrt(xdist*xdist+ydist*ydist+zdist*zdist);

            // Warning, we take into account hmax/2 around each grain
            if(dist<=r1+hmax/2.0+rad)
            {
                if(ncontact>=15)printf("idpart %d dist %f ncontact %d \n",idpart,dist,ncontact);
                //if(ncontact>=geom->max_cont_per_part){fprintf(stdout,"** Warning!\n   increase MAXCONT idpart %d \n",idpart);fflush(stdout);exit(EXIT_FAILURE);} TODO make fprintf work
                list_part[ncontact]=idparttest;
                ncontact++;
            }
        }
    }
    //if(idpart==1){for(k=0;k<icontactpot;k++) {printf("%d ",list_part[k]);}printf("\n %d \n",ncontactpot);}
    //printf("\n");
    return(ncontact);
}

__device__ int list_part_pot_sorting(unsigned int* array, int nelement)
{
    int i, j,iele;
    unsigned int elementInsere;

    // Insertion sorting of the unsigned interger array
    for (i = 1; i < nelement; i++) {
        /* Stockage de la valeur en i */
        elementInsere = array[i];
        /* Decale les elements situes avant array[i] vers la droite
           jusqu'Ã  trouver la position d'insertion */
        for (j = i; j > 0 && array[j - 1] > elementInsere; j--) {
            array[j] = array[j - 1];
        }
        /* Insertion de la valeur stockee   la place vacante */
        array[j] = elementInsere;
    }
    //** Remove the double and 0 values
    iele=0;
    // Pass the single element at the left of the array
    for(i=0;i<nelement;i++){
        if(array[i]!=array[i+1]){
            array[iele]=array[i];
            iele++;}
    }
    // Erase the double value array remainded in the array at the right side
    for(i=iele;i<nelement;i++) array[i]=0;

    // return the number of sorted element
    return iele;
}