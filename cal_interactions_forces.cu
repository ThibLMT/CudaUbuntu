//
// Created by ThibLMT on 01/04/2021.
//

#include "cal_interactions_forces.cuh"
#include "particle_interactions.cuh"
#include "sub_domain.cuh"
#include "Lock.cuh"

__device__ volatile int sem = 0;

__device__ void acquire_semaphore(volatile int *lock){
    while (atomicCAS((int *)lock, 0, 1) != 0);
}

__device__ void release_semaphore(volatile int *lock){
    *lock = 0;
    __threadfence();
}

__global__ void cal_interaction_forces(discrete_elt *particle,geom_struct *geom,unsigned int *backgrid)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    double xpartx,xparty,xpartz,radpart;
    double hmax=0;
    unsigned int idpartj;
    unsigned int ncontactperpart,icontact,icontactbis,ilist;
    unsigned int *list_part;
    vector forceji,torqueji;
    Flag test_contact,test_new_contact;
    for (int idparti = index; idparti <= geom->nb_part; idparti += stride) {


        radpart = particle[idparti].radius;


        // detect the particles in contact with the particle idparti and store the id in list_part[] array and the number of contacts per particle
        list_part = static_cast<unsigned int *>(malloc(geom->max_cont_per_part * sizeof(unsigned int)));
        memset(list_part, 0, geom->max_cont_per_part * sizeof(unsigned int));

        ncontactperpart = detect_contact_sph_backgrid(particle, radpart + hmax / 2.0, geom, idparti, list_part,
                                                      backgrid);
        // ncontactperpart=detect_contact_sphere_backgrid(xpartx,xparty,xpartz,radpart+hmax/2.0,idparti,list_part,backgrid);
//printf("id test ncontact %d\n",idparti,ncontactperpart);
        // Loop foor on the if of particles in contact with the particle	idparti
        // Computing the interaction forces for each contact
        for (ilist = 0; ilist < ncontactperpart; ilist++) {
            icontact = 0;
            idpartj = list_part[ilist];
            // test if new contact and grep the position i contact and the identifiant in particle[ni].contact[icontact] array
            test_new_contact = false;
            test_contact = false;
            do {
                if (particle[idparti].contact[icontact] == idpartj) {
                    // Contact already exist - mark the id as negative number to keep the contact in the second loop
                    particle[idparti].contact[icontact] = -idpartj;
                    test_contact = true;
                } else if (particle[idparti].contact[icontact] == 0) {
                    // New contact - mark the id as negative number to keep the contact in the second loop
                    test_contact = true;
                    test_new_contact = true;
                }
                //printf("\n icont %d particle[idparti].contact[icontact] %d \n",icontact,particle[idparti].contact[icontact]);
                icontact += 1;

            } while (!test_contact);
            // decriment icontact
            icontact -= 1;
            // mark the id of new contact as negative number to keep the contact in the second loop
            if (test_new_contact) { particle[idparti].contact[icontact] = -idpartj; }
            // Compute interaction force between parti and partj
            //!< particle_interactions(idparti,idpartj,icontact,&(forceji),&(torqueji),particle,geom);

            // adpate the sum Fext->i
            particle[idparti].Fi.x += forceji.x;
            particle[idparti].Fi.y += forceji.y;
            particle[idparti].Fi.z += forceji.z;

            particle->Mi.x += torqueji.x;
            particle->Mi.y += torqueji.y;
            particle->Mi.z += torqueji.z;

            // Add the balanced force for boundary conditions
            //if(idpartj>geom.nb_part){printf("i %d j %d f %e %e %f \n",idparti,idpartj,forceji.x,forceji.y,forceji.z);}
            if (idpartj > geom->nb_part) {
                particle[idpartj].lock.lock();
                particle[idpartj].Fi.x -= forceji.x;
                particle[idpartj].Fi.y -= forceji.y;
                particle[idpartj].Fi.z -= forceji.z;

                particle[idpartj].Mi.x += torqueji.x;
                particle[idpartj].Mi.y += torqueji.y;
                particle[idpartj].Mi.z += torqueji.z;
                particle[idpartj].lock.unlock();
            }
        }
        free(list_part);

        // remove contact from list, for particles that have been separated during the last time step
        icontact = 0;
        do {
            if (particle[idparti].contact[icontact] < 0) {
                particle[idparti].contact[icontact] = -particle[idparti].contact[icontact];
                icontact += 1;
            } else if (particle[idparti].contact[icontact] > 0) {
                icontactbis = icontact;

                do {
                    particle[idparti].contact[icontactbis] = particle[idparti].contact[icontactbis + 1];
                    particle[idparti].ut[icontactbis] = particle[idparti].ut[icontactbis + 1];
                    particle[idparti].ftanold[icontactbis] = particle[idparti].ftanold[icontactbis + 1];
                    particle[idparti].nijold[icontactbis] = particle[idparti].nijold[icontactbis + 1];
                    icontactbis += 1;
                } while (particle[idparti].contact[icontactbis] != 0);
            }

        } while (particle[idparti].contact[icontact] != 0);
    }
}