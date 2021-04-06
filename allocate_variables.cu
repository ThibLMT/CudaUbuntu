//
// Created by ThibLMT on 22/03/2021.
//

#include "allocate_variables.cuh"
#include "Lock.cuh"

__global__ void initialize_particle(discrete_elt *particle, geom_struct *geom) {
    unsigned int nelement = geom->nb_part + geom->nb_bc; //!< Number of elements of particle arrays
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    unsigned int k;
    for (int i = index; i < nelement; i += stride) {
        particle[i].Ri.x = 0.0;    // vectoru
        particle[i].Ri.y = 0.0;
        particle[i].Ri.z = 0.0;
        particle[i].Vi.x = 0.0;    // vector
        particle[i].Vi.y = 0.0;
        particle[i].Vi.z = 0.0;
        particle[i].Ai.x = 0.0;    // vector
        particle[i].Ai.y = 0.0;    // vecto
        particle[i].Ai.z = 0.0;    // vecto
        particle[i].Wi.x = 0.0;        // vector
        particle[i].Wi.y = 0.0;        // vector
        particle[i].Wi.z = 0.0;        // vector
        particle[i].Aroti.x = 0.0;    // vector
        particle[i].Aroti.y = 0.0;
        particle[i].Aroti.z = 0.0;
        particle[i].Fi.x = 0.0;   // vector
        particle[i].Fi.y = 0.0;
        particle[i].Fi.z = 0.0;
        particle[i].Mi.x = 0.0;        // vector
        particle[i].Mi.x = 0.0;
        particle[i].Mi.x = 0.0;
        particle[i].radius = 0.0;        // double
        particle[i].mass = 0.0;        // double
        particle[i].inertia = 0.0;    // double
        particle[i].Yn = 0.0;    //double
        particle[i].Nu = 0.0;        // double
        particle[i].Ndamp = 0.0;        // double
        particle[i].Mu = 0.0;    // double
        particle[i].Mur = 0.0;    // double
        particle[i].next = 0;    // int
        particle[i].clust = 0;    // int
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].type[k] = 0; }    // I
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].contact[k] = 0; }    // Int array
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].ut[k] = 0.0; }            // double array
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].ftanold[k].x = 0.0; }    //vector array
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].ftanold[k].y = 0.0; }    //vector array
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].ftanold[k].z = 0.0; }    //vector array
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].nijold[k].x = 0.0; }    //vector array
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].nijold[k].y = 0.0; }
        for (k = 0; k < geom->max_cont_per_part; k++) { particle[i].nijold[k].z = 0.0; }
    }

}

void give_properties_particle(discrete_elt *particle, double unity, material_data properties) {
    double radius;
    double mass;
    discrete_elt particle_i;
    // Mass of particle (Kg)

    particle_i = *particle;
    radius = particle_i.radius;
    mass = 4.0 / 3.0 * PI * radius * radius * radius * properties.density;
    particle_i.mass = mass;
    // Internia of particle (Kg.m2)
    particle_i.inertia = 2.0 / 5.0 * mass * radius * radius * unity * unity;
    // Elastic properties
    particle_i.Yn = properties.E;
    particle_i.Nu = properties.nu;
    particle_i.Ndamp = properties.cn;
    particle_i.Mu = properties.mu_gg;
    particle_i.Mur = properties.mu_roll_gg;
    *particle = particle_i;
}

__global__ void set_forces_0(discrete_elt *particle, geom_struct *geom) {
    unsigned int nelement = geom->nb_part + geom->nb_bc;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < nelement; i += stride) {
        particle[i].Fi.x = 0.0;   // vector
        particle[i].Fi.y = 0.0;
        particle[i].Fi.z = 0.0;
        particle[i].Mi.x = 0.0;    // vector
        particle[i].Mi.x = 0.0;
        particle[i].Mi.x = 0.0;
    }
}

__global__ void apply_gravity(discrete_elt *particle, geom_struct *geom, vect gravity) {
    unsigned int nelement = geom->nb_part;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < nelement; i += stride) {
        // Gravity force
        particle[i].Fi.x += gravity.x * particle[i].mass;
        particle[i].Fi.y += gravity.y * particle[i].mass;
        particle[i].Fi.z += gravity.z * particle[i].mass;
    }
}

__global__ void update_particle(discrete_elt *particle, geom_struct *geom) {
    unsigned int nelement = geom->nb_part;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < nelement; i += stride) {
        double deltat;
        deltat = geom->deltat;
        //**********
        //  Update particle position for the next ime step - loop on the particle
        // Semi-implicite euler scheme
        // ai(t)=1/mi*sum(Fext->i)
        // vi(t+dt)=vi(t)+ai(t)*dt
        // xi(t+dt)=xi(t)+vi(t+dt)*dt
        // awi(t)=1/Ii*sum(Text->i)
        // wi(t+dt)=wi(t)+awi(t)*dt

        // ai(t)=1/mi*sum(Fext->i)
        particle[i].Ai.x = particle[i].Fi.x / geom->unity / particle[i].mass;
        particle[i].Ai.y = particle[i].Fi.y / geom->unity / particle[i].mass;
        particle[i].Ai.z = particle[i].Fi.z / geom->unity / particle[i].mass;

        // vi(t+dt)=vi(t)+ai(t)*dt
        particle[i].Vi.x = (particle[i].Vi.x + particle[i].Ai.x * deltat);
        particle[i].Vi.y = (particle[i].Vi.y + particle[i].Ai.y * deltat);
        particle[i].Vi.z = (particle[i].Vi.z + particle[i].Ai.z * deltat);

        // Computing of angular acceleration -> awi(t)=1/Ii*sum(Text->i)
        particle[i].Aroti.x = particle[i].Mi.x * geom->unity / particle[i].inertia;
        particle[i].Aroti.y = particle[i].Mi.y * geom->unity / particle[i].inertia;
        particle[i].Aroti.z = particle[i].Mi.z * geom->unity / particle[i].inertia;

        // actualization of angular velocity -> wi(t+dt)=wi(t)+awi(t)*dt
        particle[i].Wi.x = particle[i].Wi.x + particle[i].Aroti.x * deltat;
        particle[i].Wi.y = particle[i].Wi.y + particle[i].Aroti.y * deltat;
        particle[i].Wi.z = particle[i].Wi.z + particle[i].Aroti.z * deltat;

        // actualization of position by checking the periodic condtition
        // xi(t+dt)=xi(t)+vi(t+dt)*dt
        particle[i].Ri.x += particle[i].Vi.x * deltat;
        if (particle[i].Ri.x > geom->sizex) { particle[i].Ri.x -= geom->sizex; }
        else if (particle[i].Ri.x < 0.0) { particle[i].Ri.x += geom->sizex; }

        particle[i].Ri.y += particle[i].Vi.y * deltat;
        if (particle[i].Ri.y > geom->sizey) { particle[i].Ri.y -= geom->sizey; }
        else if (particle[i].Ri.y < 0.0) { particle[i].Ri.y += geom->sizey; }

        particle[i].Ri.z += particle[i].Vi.z * deltat;
        if (particle[i].Ri.z > geom->sizez) { particle[i].Ri.z -= geom->sizez; }
        else if (particle[i].Ri.z < 0.0) { particle[i].Ri.z += geom->sizez; }
    }
}