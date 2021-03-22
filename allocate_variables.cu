//
// Created by ThibLMT on 22/03/2021.
//

#include "allocate_variables.cuh"

__global__ void initialize_particle(discrete_elt *particle, geom_struct *geom)
{
    unsigned int nelement=geom->nb_part+geom->nb_bc; //!< Number of elements of particle arrays
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    unsigned int k;
    for (int i = index; i < nelement; i+= stride)
    {
        particle[i].Ri.x=0.0;	// vectoru
        particle[i].Ri.y=0.0;
        particle[i].Ri.z=0.0;
        particle[i].Vi.x=0.0;	// vector
        particle[i].Vi.y=0.0;
        particle[i].Vi.z=0.0;
        particle[i].Ai.x=0.0;	// vector
        particle[i].Ai.y=0.0;	// vecto
        particle[i].Ai.z=0.0;	// vecto
        particle[i].Wi.x=0.0;		// vector
        particle[i].Wi.y=0.0;		// vector
        particle[i].Wi.z=0.0;		// vector
        particle[i].Aroti.x=0.0; 	// vector
        particle[i].Aroti.y=0.0;
        particle[i].Aroti.z=0.0;
        particle[i].Fi.x=0.0;   // vector
        particle[i].Fi.y=0.0;
        particle[i].Fi.z=0.0;
        particle[i].Mi.x=0.0;		// vector
        particle[i].Mi.x=0.0;
        particle[i].Mi.x=0.0;
        particle[i].radius=0.0;		// double
        particle[i].mass=0.0;		// double
        particle[i].inertia=0.0;	// double
        particle[i].Yn=0.0;	//double
        particle[i].Nu=0.0;		// double
        particle[i].Ndamp=0.0;		// double
        particle[i].Mu=0.0;	// double
        particle[i].Mur=0.0;	// double
        particle[i].next=0;	// int
        particle[i].clust=0;    // int
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].type[k]=0;}	// I
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].contact[k]=0;}	// Int array
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].ut[k]=0.0;}			// double array
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].ftanold[k].x=0.0;}	//vector array
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].ftanold[k].y=0.0;}	//vector array
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].ftanold[k].z=0.0;}	//vector array
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].nijold[k].x=0.0;}	//vector array
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].nijold[k].y=0.0;}
        for(k=0;k<geom->max_cont_per_part;k++){particle[i].nijold[k].z=0.0;}
    }

}