//
// Created by ThibLMT on 01/04/2021.
//

//#include "particle_interactions.cuh"
// TODO implement functions here

/**
 *\file particle_interactions.c
 *\brief body of the function particle_interactions - (previous name choc())
 */




#include "particle_interactions.cuh"
#include "dlvo_forces.cuh"
#include "tangential_force_mindlin_vect.cuh"
#include "small_vect_rot.cuh"
#include "rolling_resistance_torque.cuh"
#include  "normal_force_hertz.cuh"

__device__ void set_vect_0(vector *vec)
{
    vec->x=0.0;
    vec->y=0.0;
    vec->z=0.0;
}
// Using these following interaction force described in main.c for the moment. The other are declared in list header files

__device__ void particle_interactions(unsigned int ni,unsigned int nj,int k_ij,vector* forceji_pointer,vector* torqueji_pointer,discrete_elt *particle,geom_struct geom)
/* Calls: no other routines */
/* Called by: integrate() */
{
    double dist_ij,dist_iwall; //!< Distance between particle i and j and between i and wall (ca
    double deltan,deltandot; //!< Overlap (deflexion) and d/dt of the overlap (local basis)
    double hij; //!< Overlap after ranging from hmin and hmax
    double normevrtij; //!< Norm of the relative velocity
    double zpr; //!<
    double ri,rj; //!< Radius of the particles
    double rij; //!< Effective radius
    // Effective material parameters
    double Eij; //!< Effective Young Modulus radius (Pa)
    double Aij; //!< Effective damping coefficient (sec)
    double Gij; //!< coefficient for the calculation of the tangential force
    double muij; //!< Friction coefficient between grain i and grain j
    double murollij; //!< Rolling friction coefficient between grain i and grain j
    // Parameters for the hcss force model
    double rieq,rjeq; //!< Intermediate geometrical parameters for the HCSS force model
    double paste_e; //!< Thickness of the paste for the HCSS force model
    // Norms of the forces. Can be used to print some history of the contacts
    double fvdw; //!< Van der Walls force (local basis)
    double felec; //!< Electrostatic force (local basis)
    double fncontact=0.0; //!< Normal contact force (local basis)
    double ftcontact=0.0; //!< Tangential contact force (local basis)
    double fnpaste=0.0; //!< Normal contact force for hcss model(local basis)
    double ftpaste_loc=0.0; //!< Normal contact force for hcss model(local basis)
    // Cinematic parameters
    vector nij; //!< Normal and tangential vectors of the local basis
    vector niwall; //!< Intermediate normal vector between particle and boundary conditions walls, cylinder, cone 1<iobj<nsurf
    vector xi,xj; //!< Coordinates particles in contact
    vector vi,vj,wi,wj; //!< Linear and angular velocities of particles
    vector vrij; //!< relative velocity vector (global basis)
    vector vrtij; //!< tangential components in the local basis of the relative velocity vector projected in the global basis
    // Force parameters
    vector ftan; //!< Tangential contact force (global basis)
    vector ftpaste; //!< Tangential contact force for the paste model (global basis)
    vector torq_rol; //!< Rolling resistance torque (global basis)
    vector forceji; //!< Total force exerced by partj on part i (global basis)
    vector torqueji; //!< Total torque exerced by partj on part i (global basis)
    unsigned int npart,nsurf;
    npart=geom.nb_part;
    nsurf=geom.nb_bc;
    double deltat;
    deltat=geom.deltat;
    int syssizex,syssizey,syssizez,syssizel;
    syssizex=geom.sizex;
    syssizey=geom.sizey;
    syssizez=geom.sizez;
    syssizel=geom.sizel;
    // Initialization of forceij and torqie ij vectors
    set_vect_0(&forceji);
    set_vect_0(&torqueji);

    // Treating the contact between grain ni and grain nj
    //**************************************************
    // Bringing information about the state of each particle in contact (ni and nj)
    // Grain ni : only particles constituting the granular media
    xi.x=particle[ni].Ri.x;
    xi.y=particle[ni].Ri.y;
    xi.z=particle[ni].Ri.z;
    ri=particle[ni].radius;

    vi.x=particle[ni].Vi.x;
    vi.y=particle[ni].Vi.y;
    vi.z=particle[ni].Vi.z;

    wi.x=particle[ni].Wi.x;
    wi.y=particle[ni].Wi.y;
    wi.z=particle[ni].Wi.z;

    // Grain nj. This grain can be a particle or a limit condition as object, cylinder or wall.
    xj.x=particle[nj].Ri.x;
    xj.y=particle[nj].Ri.y;
    xj.z=particle[nj].Ri.z;
    rj=particle[nj].radius;

    /*

    // If grain nj is a cylinder condition iobj=7
    if(nj==npart+7)
      {
        niwall.x=xi.x-syssizex/2.0;
        niwall.y=xi.y-syssizey/2.0;
        dist_iwall=sqrt(niwall.x*niwall.x+niwall.y*niwall.y);
        niwall.x=niwall.x/dist_iwall;
        niwall.y=niwall.y/dist_iwall;
        xj.x=syssizex/2.0+niwall.x*(particle[nj].radius+diamcyl/2.0);
        xj.y=syssizey/2.0+niwall.y*(particle[nj].radius+diamcyl/2.0);
        xj.z=xi.z;
        rj=particle[nj].radius;
      }

    // If grain nj is a cone condition iobj=8
    if(nj==npart+8)
      {
        rj=particle[nj].radius;

        //zpr=hcone of point P' located at the intersection of cone axis and line going through ni and perpendicular to cone surface
        zpr=xi.z-(rbcone-rtcone)/hcone*sqrt((xi.x-syssizex/2.0)*(xi.x-syssizex/2.0)+(xi.y-syssizey/2.0)*(xi.y-syssizey/2.0));

        //n=unit vector perpendicular to cone surface and going through ni (vector colinear to P'ni)
        niwall.x=xi.x-syssizex/2.0;
        niwall.y=xi.y-syssizey/2.0;
        niwall.z=xi.z-zpr;
        dist_iwall=sqrt(niwall.x*niwall.x+niwall.y*niwall.y+niwall.z*niwall.z);
        niwall.x=niwall.x/dist_iwall;
        niwall.y=niwall.y/dist_iwall;
        niwall.z=niwall.z/dist_iwall;

        //distance between P' and point on cone surface in the direction n
        dist_iwall=(rbcone/(rbcone-rtcone)*hcone-zpr)*sin(atan((rbcone-rtcone)/hcone));

        //distance between P' and center of nj
        dist_iwall=dist_iwall+rj;

        xj.x=syssizex/2.0+dist_iwall*niwall.x;
        xj.y=syssizey/2.0+dist_iwall*niwall.y;
        xj.z=zpr+dist_iwall*niwall.z;
      }
  */
    // Checking the periodic boundary condition excepted for all surfaces (walls, cylinder, cone 1<iobj<nsurf
    if((nj<=npart)||(nj>npart+nsurf))
    {
        /* periodic boundary conditions for xj.x */
        if((xj.x>xi.x)&&(xj.x-xi.x>xi.x+syssizex-xj.x)){xj.x-=syssizex;}
        else if((xj.x<xi.x)&&(xi.x-xj.x>xj.x+syssizex-xi.x)){xj.x+=syssizex;}
        /* periodic boundary conditions for xj.y */
        if((xj.y>xi.y)&&(xj.y-xi.y>xi.y+syssizey-xj.y)){xj.y-=syssizey;}
        else if((xj.y<xi.y)&&(xi.y-xj.y>xj.y+syssizey-xi.y)){xj.y+=syssizey;}
        /* periodic boundary conditions for xj.z */
        if((xj.z>xi.z)&&(xj.z-xi.z>xi.z+syssizez-xj.z)){xj.z-=syssizez;}
        else if((xj.z<xi.z)&&(xi.z-xj.z>xj.z+syssizez-xi.z)){xj.z+=syssizez;}
    }

    vj.x=particle[nj].Vi.x;
    vj.y=particle[nj].Vi.y;
    vj.z=particle[nj].Vi.z;

    wj.x=particle[nj].Wi.x;
    wj.y=particle[nj].Wi.y;
    wj.z=particle[nj].Wi.z;

    // computation of the normal unit vector (ni,nj), pointing from the center of the grain ni toward the center of the grain nj
    nij.x=xj.x-xi.x;
    nij.y=xj.y-xi.y;
    nij.z=xj.z-xi.z;

    dist_ij=sqrt(nij.x*nij.x+nij.y*nij.y+nij.z*nij.z);

    nij.x=nij.x/dist_ij;
    nij.y=nij.y/dist_ij;
    nij.z=nij.z/dist_ij;

    // computation of the relative velocity at the contact point I1 versus I2
    vrij.x=vi.x-vj.x+((wi.y*ri+wj.y*rj)*nij.z-(wi.z*ri+wj.z*rj)*nij.y);
    vrij.y=vi.y-vj.y+((wi.z*ri+wj.z*rj)*nij.x-(wi.x*ri+wj.x*rj)*nij.z);
    vrij.z=vi.z-vj.z+((wi.x*ri+wj.x*rj)*nij.y-(wi.y*ri+wj.y*rj)*nij.x);

    // computation of normal indentation (overlap)
    deltan=ri+rj-dist_ij;

    // computation of the normal component of relative velocity
    deltandot=vrij.x*nij.x+vrij.y*nij.y+vrij.z*nij.z;

    // Effective radius
    rij=(ri*rj)/(ri+rj);

    // Effective Young modulus
    Eij=((1.0-pow(particle[ni].Nu,2.0))/particle[ni].Yn)+((1.0-pow(particle[nj].Nu,2.0))/particle[nj].Yn);
    Eij=1.0/Eij;

    // Effective damping coefficient
    Aij=(particle[ni].Ndamp+particle[nj].Ndamp)/2.0;


    // Computation of the interaction between particles
    // sign convention
    // The unit vector points from grain ni toward grain nj. We compute the force excerced by grain nj on Grain ni.

    // Van der Waals forces
    /*
    if (forceVdw)
      {
        // The van der Wall forces are attractives, the scalar product n.F must be positive.
        // The van der Walls forces are computed are significant at a distance ranging from hmin to hmax.
        if(-deltan<=hmax)
      {
        hij=fmax(-deltan,hmin);
        fvdw=van_der_waals_force(hij,ri,rj,rij,modelefvdw);

      }
        else
      {fvdw=0.0;}

        //Contribution of the van der Walls forces on particle ni
        forceji.x+=fvdw*nij.x;
        forceji.y+=fvdw*nij.y;
        forceji.z+=fvdw*nij.z;
      }
    else
      {fvdw=0.0;}
    //Electrostatic force

    if(forceElectro)
      {
        //The electrostatic force is repulsive force. This normal component must be negative.

        felec=NormalElectrostaticForce(-deltan,ri,rj,rij);

        // add normal force for ni and nj
        forceji.x+=felec*nij.x;
        forceji.y+=felec*nij.y;
        forceji.z+=felec*nij.z;


      }
    else{felec=0.0;}

    // HCSS force model (Hard Core Soft Shell model)
    if (forceHCSS)
      {
        // Firstly, computation of the geometrical parameters to determine the equivalent radius needed for the torque
        if(deltan+hmax>=0)
      {
        // Warning fnpaste and ftpaste_loc are computed in the local basis and ftpaste in the global basis
        // Thickness of the paste equal to hmax/2 by default.
        // Warning, define also as local variable in tangential_force_hcss function
        paste_e=hmax/2.0;
        // Equivalent radius inner the paste at the contact between two particles
        rieq=(dist_ij*dist_ij+(ri+paste_e)*(ri+paste_e)-(rj+paste_e)*(rj+paste_e))/(2.0*dist_ij);
        rjeq=dist_ij-rieq;

        // Normal component of the HCSS force
        fnpaste=normal_force_hcss(deltan);

        // Tangential component of the HCSS force
        ftpaste=tangential_force_hcss(deltan,nij,particle[ni].radius,rieq,particle[ni].Vi,particle[ni].Wi,particle[nj].radius,rjeq,particle[nj].Vi,particle[nj].Wi);
        ftpaste_loc=sqrt(ftpaste.x*ftpaste.x+ftpaste.y*ftpaste.y+ftpaste.z*ftpaste.z);

        //exclude sliding on the bottom plane if rheometer simulation
        //exclude sliding for cylinder and cone boundary conditions
        if(((ChoixTypeCalcul==2)&&(nj==npart+1))||(nj==npart+7)||(nj==npart+8)){
          ftpaste=vect0;
        }

        forceji.x+=fnpaste*nij.x+ftpaste.x;
        forceji.y+=fnpaste*nij.y+ftpaste.y;
        forceji.z+=fnpaste*nij.z+ftpaste.z;

        torqueji.x+=rieq*(nij.y*ftpaste.z-nij.z*ftpaste.y)-4.0*paste_consistancy*PI*particle[ni].radius*particle[ni].Wi.x*particle[ni].radius*particle[ni].radius;
        torqueji.y+=rieq*(nij.z*ftpaste.x-nij.x*ftpaste.z)-4.0*paste_consistancy*PI*particle[ni].radius*particle[ni].Wi.y*particle[ni].radius*particle[ni].radius;
        torqueji.z+=rieq*(nij.x*ftpaste.y-nij.y*ftpaste.x)-4.0*paste_consistancy*PI*particle[ni].radius*particle[ni].Wi.z*particle[ni].radius*particle[ni].radius;

      }
        else{fnpaste=0.0;ftpaste_loc=0.0;ftpaste=vect0;}
      }
    else{fnpaste=0.0;ftpaste_loc=0.0;ftpaste=vect0;}*/

    // computation of the contact forces

    ftan=particle[ni].ftanold[k_ij];

    if (/*(forceContact)&&*/(deltan>0.0))
    {
        // computation of the normal contact force exerted by nj on ni
        // Computing the Hertz non-linear elastic force + dissipative component proposed by T. Pöschel and T. Schwager,eq 2.14
        fncontact=normal_force_hertz(deltan,deltandot,rij,Eij,Aij);

        // Adding normal force for ni and nj
        forceji.x+=fncontact*nij.x;
        forceji.y+=fncontact*nij.y;
        forceji.z+=fncontact*nij.z;

        // Computation of the tangential component of relative velocity */

        vrtij.x=vrij.x-deltandot*nij.x;
        vrtij.y=vrij.y-deltandot*nij.y;
        vrtij.z=vrij.z-deltandot*nij.z;
        normevrtij=sqrt(vrtij.x*vrtij.x+vrtij.y*vrtij.y+vrtij.z*vrtij.z);

        if((normevrtij>0.0)&&(fncontact!=0.0))
        {

            // Note that fncontact can be wanish during small unloading contact due to dashpot term
            // Increment the number of contact only if fncontact<0
            // computation of the scalar accumulated tangential displacement of particles ni and nj
            particle[ni].ut[k_ij]+=normevrtij*deltat;

            // Gij coefficient
            Gij=((2.0-particle[ni].Nu)*(1.0+particle[ni].Nu)/particle[ni].Yn)+((2.0-particle[nj].Nu)*(1.0+particle[nj].Nu)/particle[nj].Yn);
            Gij=1.0/(2.0*Gij);

            // Friction coefficient and check if it is contact grain-grain or grain-wall
            muij=particle[ni].Mu;
            if(nj>npart){muij=particle[nj].Mu;}

            // Computation of the tangential (Mindlin and Deresiewik expression, limited by the Coulomb criteria
/*
	  if(tan_force_model==1)
	    {
	      ftcontact=tangential_force_mindlin(deltan,particle[ni].ut[k_ij],fabs(fncontact),muij,Eij,Gij);

	      ftan.x=ftcontact*vrtij.x/normevrtij;
	      ftan.y=ftcontact*vrtij.y/normevrtij;
	      ftan.z=ftcontact*vrtij.z/normevrtij;
	    }
	  else if (tan_force_model==2)
	    {*/
            ftan=small_vect_rot(particle[ni].ftanold[k_ij],nij,particle[ni].nijold[k_ij],particle[ni].Wi,particle[nj].Wi,geom.deltat);

            ftan=tangential_force_mindlin_vect(deltan,fabs(fncontact),muij,rij,Gij,ftan,vrtij,geom.deltat);

            ftcontact=sqrt(ftan.x*ftan.x+ftan.y*ftan.y+ftan.z*ftan.z);
            // }

            // add tangential force for ni and nj
            forceji.x+=ftan.x;
            forceji.y+=ftan.y;
            forceji.z+=ftan.z;

            // computation of the momentum exerted by Ft on particle ni
            torqueji.x+=particle[ni].radius*(nij.y*ftan.z-nij.z*ftan.y);
            torqueji.y+=particle[ni].radius*(nij.z*ftan.x-nij.x*ftan.z);
            torqueji.z+=particle[ni].radius*(nij.x*ftan.y-nij.y*ftan.x);

            // computation of the momentum exerted by Ft on particle nj

        }

        else
        {
            ftcontact=0.0;
            ftan.x=0;
            ftan.y=0;
            ftan.z=0;
        }

        // If rolling friction
/*
     if(torqueRollRes!=0)
	{
	  murollij=particle[ni].Mur;
	  if(nj>npart){murollij=particle[nj].Mur;}
	  torq_rol=rolling_resistance_torque(fabs(fncontact),murollij,rij,torq_rol,particle[ni].Wi,particle[nj].Wi);
	  torqueji.x-=torq_rol.x;
	  torqueji.y-=torq_rol.y;
	  torqueji.z-=torq_rol.z;
	}
      else{set_vect0(torq_rol);}
*/
    }
    else
    {
        // No contact force
        fncontact=0.0;
        ftcontact=0.0;
        ftan.x=0;
        ftan.y=0;
        ftan.z=0;
    }


    // Return force and torque exerced by part j on part i
    *forceji_pointer=forceji;
    *torqueji_pointer=torqueji;

    // Save the tangential force at the local state for the next time step
    particle[ni].ftanold[k_ij]=ftan;
    particle[ni].nijold[k_ij]=nij;


    return ;
}