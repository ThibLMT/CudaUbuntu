#include <iostream>

#include "def_types.h"
#include "def_global_variables.h"
#include "allocate_variables.cuh"
#include "read_micro.cuh"
#include "init_params.cuh"
#include "sub_domain.cuh"

#define NB_PART 200
#define SYSSIZEX 5
#define SYSSIZEY 5
#define SYSSIZEZ 14
#define UNITY 0.010000
#define BOUNDARYCONT 20


int main() {
    char Nlogfile[50]="logfile";
    char Nmicrofile[50]="init_micro";
    // Initialization of some global variables
    //** Flag variables
    Nparafile=(char*)"params.in";
    int i,j,k,k1;
    int nb_taches;
    int iter,niter,imicro;
    char filename[50];
    discrete_elt *particle;
    geom_struct *geom;
    unsigned int *backgrid = nullptr;
    int *backgrid_insert = nullptr;

    // Initialization of Ierror
    ierror=EXIT_SUCCESS;

    cudaMallocManaged(&geom, sizeof(geom_struct));

    // * Get the number of particles
    geom->nb_part = NB_PART;
    // * Get the number of boundary contacts
    geom->nb_bc = BOUNDARYCONT;
    read_geom_param(geom);

    // Allocate the particle array
    int nb_elements = geom->nb_part + geom->nb_bc + 1;
    cudaMallocManaged(&particle,nb_elements * sizeof(discrete_elt));

    int blockSize = 256;
    int numBlocks = (nb_elements + blockSize - 1)/blockSize;

    // Sets all particle members to 0
    initialize_particle<<<numBlocks,blockSize>>>(particle,geom);
    cudaDeviceSynchronize();

    microfile_read_particle(Nmicrofile,particle,geom);
    microcontfile_read_contact(Nmicrofile,particle,geom);

    //********************
    // Initialize parameters
    // Bulk parameters
    read_table_mat(&prop_mat_part);

    // Friction parameters
    prop_mat_part.mu_gg=0.3; //!< Friction coefficient grain-grain
    prop_mat_part.mu_gw=0.3; //!< Friction coefficient grain-wall
    // -- Rolling resistant parameters
    prop_mat_part.mu_roll_gg=0.01; //!< Rolling friction coefficient grain-grain
    prop_mat_part.mu_roll_gw=0.01;
    // Set bulk forces
    gravity.x=0.0;
    gravity.y=0.0;
    gravity.z=-9.81;  // m.s-2

    //adimention of length
    adi_params(&prop_mat_part,geom);
    // give properties for particle
    for(i=1;i<=geom->nb_part;i++)
    {
        give_properties_particle(&particle[i],geom->unity,prop_mat_part);
    }

    // set parameter for the plane
    i=geom->nb_part+1;
    give_properties_particle(&particle[i],geom->unity,prop_mat_part);
    particle[i].radius=1000000.0;
    particle[i].Ri.x=geom->sizex/2.0;
    particle[i].Ri.y=geom->sizey/2.0;
    particle[i].Ri.z=0.0 + -1.0 * particle[i].radius;


    // Start DEM computation
    // Allocation of subdomain backgrid
    nb_elements=geom->sizex*geom->sizey*geom->sizez*geom->sizel;
    cudaMallocManaged(&backgrid,nb_elements * sizeof(backgrid));
    nb_elements=geom->sizex*geom->sizey*geom->sizez;
    cudaMallocManaged(&backgrid_insert,nb_elements * sizeof(backgrid_insert));


    initialize_backgrid<<<numBlocks,blockSize>>>(backgrid,backgrid_insert,geom);
    cudaDeviceSynchronize();

    for(i=0;i<geom->sizex;i++)
    {
        for(j=0;j<geom->sizey;j++)
        {
            set_id_backgrid(i,j,0,geom->nb_part+1,backgrid,backgrid_insert,geom);
        }

    }

    geom->deltat=0.000001;
    niter=100000;
    imicro=0;
    iter=0;

    do{
        // Reset the forces and moments on the particles
        set_forces_0<<<numBlocks,blockSize>>>(particle,geom);
        cudaDeviceSynchronize();

        // Apply gravity
        apply_gravity<<<numBlocks,blockSize>>>(particle,geom,gravity);
        cudaDeviceSynchronize();

        // Insert the spheres in the background
        insert_sph_backgrid<<<numBlocks,blockSize>>>(particle,backgrid,backgrid_insert,geom);
        cudaDeviceSynchronize();
        iter++;
    }
    while(iter<=niter);

    // Frees allocated memory
    cudaFree(particle);
    cudaFree(geom);
    cudaFree(backgrid);
    cudaFree(backgrid_insert);
    return 0;
}

