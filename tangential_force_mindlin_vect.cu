//
// Created by bluet on 06/04/2021.
//

/**
*\file tangential_force_mindlin.c
*\brief body of the function tangential_force_mindlin
*/


#include "def_global_variables.h"
#include "tangential_force_mindlin_vect.cuh"

__device__ vector tangential_force_mindlin_vect(double deltanl,double fnijl,double mul, double Rij, double Gijl,vector ftan,vector vtanrel,double deltat)
{
    double kstiff;
    double fcoulomb;
    double normefreturn;
    double normeftan;

    vector freturn;

    // The tangential force is given by the equation of Mindlin and Deriewsic in the case where the normal force is constant
    // with the increasing of the tangential force. deltatamax is the value of maximal tangential displacement
    // corresponding to total sliding at the contact (Coulomb criteria is activated)
    //First compute the Coulomb criteria
    fcoulomb=fabs(mul*fnijl);
    kstiff=sqrt(deltanl*Rij);
    kstiff=kstiff*8*Gijl;

    freturn.x=ftan.x-kstiff*vtanrel.x*deltat;
    freturn.y=ftan.y-kstiff*vtanrel.y*deltat;
    freturn.z=ftan.z-kstiff*vtanrel.z*deltat;
    normefreturn=sqrt(freturn.x*freturn.x+freturn.y*freturn.y+freturn.z*freturn.z);

    normeftan=sqrt(ftan.x*ftan.x+ftan.y*ftan.y+ftan.z*ftan.z);
// test the Coulomb criteria - Need to be sure that ftan is not vanished (like first contact)
    if ((normefreturn > fcoulomb)&&(normeftan>0.0))
    {
        /*
        freturn.x=fcoulomb*freturn.x/normefreturn;
        freturn.y=fcoulomb*freturn.y/normefreturn;
        freturn.z=fcoulomb*freturn.z/normefreturn;
        */

        freturn.x=fcoulomb*ftan.x/normeftan;
        freturn.y=fcoulomb*ftan.y/normeftan;
        freturn.z=fcoulomb*ftan.z/normeftan;
        ncontslide++;
    }
    normefreturn=sqrt(freturn.x*freturn.x+freturn.y*freturn.y+freturn.z*freturn.z);

    indfric= indfric + (normefreturn/fcoulomb);
    ftan=freturn;

    return(ftan);
}