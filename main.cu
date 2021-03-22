#include <iostream>

#include "def_types.h"
#include "def_global_variables.h"
#include "allocate_variables.cuh"

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
    geom_sys *geom;

    cudaMallocManaged(&geom, sizeof(geom_sys));

    // * Get the number of particles
    geom->nb_part = NB_PART;
    // * Get the number of boundary contacts
    geom->nb_bc = BOUNDARYCONT;

    // Allocate the particle array
    int nb_elements = geom->nb_part + geom->nb_bc + 1;
    cudaMallocManaged(&particle,nb_elements * sizeof(discrete_elt));

    int blockSize = 256;
    int numBlocks = (nb_elements + blockSize - 1)/blockSize;

    // Sets all particle members to 0
    initialize_particle<<<numBlocks,blockSize>>>(particle,geom);

    // Frees allocated memory
    cudaFree(particle);
    cudaFree(geom);
    return 0;
}
