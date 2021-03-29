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

void set_id_backgrid(int xv,int yv,int zv,unsigned int idp,unsigned int *backgrid, int *backgrid_insert,geom_struct *geom)
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
    for (int loopI = index; loopI < geom->nb_part; loopI+= stride) {
        double xpart = particle[loopI].Ri.x;
        double ypart = particle[loopI].Ri.y;
        double zpart = particle[loopI].Ri.z;
        double radpart = particle[loopI].radius;
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
                        backgrid[index_backgrid] = loopI;




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
