#include <iostream>

#include "def_types.h"
#include "def_global_variables.h"
#define NB_PART 200
#define SYSSIZEX 5
#define SYSSIZEY 5
#define SYSSIZEZ 14
#define UNITY 0.010000
#define BOUNDARYCONT 20;


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
    discrete_elt **particle;
    geom_sys *geom;

    cudaMallocManaged(&geom, sizeof(geom));

    // * Get the number of particles
    geom->nb_part = NB_PART;
    // * Get the number of boundary contacts
    geom->nb_bc = BOUNDARYCONT;

    // TODO Implement set_vector function in CUDA

    std::cout << "Hello, World!" << std::endl;
    return 0;
}
