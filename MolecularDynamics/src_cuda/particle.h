#pragma once 
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>

__global__ void initCurandStates(curandState* states, Domain domain){
    int ti = blockDim.x*blockIdx.x + threadIdx.x;
    int tj = blockDim.y*blockIdx.y + threadIdx.y;
    int tk = blockDim.z*blockIdx.z + threadIdx.z;
    
    if (ti >= domain.n_particles_x || tj >= domain.n_particles_y || tk >= domain.n_particles_z) {
        return;
    }
    const int idx = ti + domain.n_particles_x*tj + domain.n_particles_x*domain.n_particles_y*tk;
    if (idx>=domain.n_particles_total) return;
    // printf("idx: %d\n",idx);
    curand_init(12345ULL,idx,0,&states[idx]);
}

__global__ void generateParticles(const Domain domain, double* __restrict__ position, 
                                                 double* __restrict__ velocity,
                                                 double* __restrict__ acceleration){
    
    const int ti = blockDim.x*blockIdx.x+threadIdx.x;
    const int tj = blockDim.y*blockIdx.y+threadIdx.y;
    const int tk = blockDim.z*blockIdx.z+threadIdx.z;
    
    //if (ti<domain.n_particles_x && tj<domain.n_particles_x && tk<domain.n_particles_x ){
    if (ti >= domain.n_particles_x || tj >= domain.n_particles_y || tk >= domain.n_particles_z) {
        return;
    }
    const int idx = ti + domain.n_particles_x*tj + domain.n_particles_x*domain.n_particles_y*tk;
    // printf("%d,%d,%d\n",ti,tj,tk);
    if (domain.domainDimension==1){
        position[3*idx+0] = ti*domain.delta_x;
        position[3*idx+1] = 0.0;
        position[3*idx+2] = 0.0;
        
        acceleration[3*idx+0] = 0.0;
        acceleration[3*idx+1] = domain.gravity;
        acceleration[3*idx+2] = 0.0;
        }
    else if (domain.domainDimension==2){
        position[3*idx+0] = ti*domain.delta_x;
        position[3*idx+1] = tj*domain.delta_x;
        position[3*idx+2] = 0.0;
        
        acceleration[3*idx+0] = 0.0;
        acceleration[3*idx+1] = domain.gravity;
        acceleration[3*idx+2] = 0.0;
        }
    else{
        position[3*idx+0] = ti*domain.delta_x;
        position[3*idx+1] = tj*domain.delta_x;
        position[3*idx+2] = tk*domain.delta_x;
        
        acceleration[3*idx+0] = 0.0;
        acceleration[3*idx+1] = domain.gravity;
        acceleration[3*idx+2] = 0.0;
        } 
    //printf("Particle Id: %d, x: %f, y: %f, z: %f\n",idx,position[3*idx+0],position[3*idx+1],position[3*idx+2]);
    }

__global__ void gaussianVelocityDistribution(curandState* globalState, double* __restrict__ velocity, const Domain domain){
    const int ti = blockDim.x*blockIdx.x+threadIdx.x;
    const int tj = blockDim.y*blockIdx.y+threadIdx.y;
    const int tk = blockDim.z*blockIdx.z+threadIdx.z;
    
    //if (ti<domain.n_particles_x && tj<domain.n_particles_x && tk<domain.n_particles_x ){
    if (ti >= domain.n_particles_x || tj >= domain.n_particles_y || tk >= domain.n_particles_z) {
        return;
    }
    const int idx = ti + domain.n_particles_x*tj + domain.n_particles_x*domain.n_particles_y*tk;
    curandState localState = globalState[idx];
    const double stddev = sqrt(domain.T);
    double vx, vy, vz;
    double max_vel = 2.5*stddev;
    do {
        vx = curand_normal_double(&localState)*stddev;
    }while (vx>max_vel || vx<-max_vel);
    do {
        vy = curand_normal_double(&localState)*stddev;
    }while (vy>max_vel || vy<-max_vel);
    do {
        vz = curand_normal_double(&localState)*stddev;
    }while (vz>max_vel || vz<-max_vel);

    if (domain.domainDimension==1){
        velocity[3*idx+0] = vx;
        velocity[3*idx+1] = 0.0;
        velocity[3*idx+2] = 0.0;
    }
    else if (domain.domainDimension==2){
        velocity[3*idx+0] = vx;
        velocity[3*idx+1] = vy;
        velocity[3*idx+2] = 0.0;
    }
    else{
        velocity[3*idx+0] = vx;
        velocity[3*idx+1] = vy;
        velocity[3*idx+2] = vz;
    }
    }

__global__ void updatePositionVelocityFirstHalf(const Domain domain, double* __restrict__ position,
                                                                    double* __restrict__ velocity,
                                                                    double* __restrict__ acceleration
                                                                    ){
    const int ti = blockDim.x*blockIdx.x+threadIdx.x;
    const int tj = blockDim.y*blockIdx.y+threadIdx.y;
    const int tk = blockDim.z*blockIdx.z+threadIdx.z; 
    // const int dim = domain.domainDimension;
    double deltaTime = domain.deltaTime;
    double box_len = domain.box_len;
    if (ti >= domain.n_particles_x || tj >= domain.n_particles_y || tk >= domain.n_particles_z) {
        return;
    }
    const int idx = ti + domain.n_particles_x*tj + domain.n_particles_x*domain.n_particles_y*tk;
    
    if (domain.domainDimension==1){
        position[3*idx+0] += velocity[3*idx+0]*deltaTime+acceleration[3*idx+0]*deltaTime*deltaTime*0.5;
        velocity[3*idx+0] += acceleration[3*idx+0]*deltaTime*0.5; 
    }
    else if (domain.domainDimension==2){
        position[3*idx+0] += velocity[3*idx+0]*deltaTime+acceleration[3*idx+0]*deltaTime*deltaTime*0.5;
        position[3*idx+1] += velocity[3*idx+1]*deltaTime+acceleration[3*idx+1]*deltaTime*deltaTime*0.5;
        velocity[3*idx+0] += acceleration[3*idx+0]*deltaTime*0.5;
        velocity[3*idx+1] += acceleration[3*idx+1]*deltaTime*0.5;
    }
    else{
        position[3*idx+0] += velocity[3*idx+0]*deltaTime+acceleration[3*idx+0]*deltaTime*deltaTime*0.5;
        position[3*idx+1] += velocity[3*idx+1]*deltaTime+acceleration[3*idx+1]*deltaTime*deltaTime*0.5;
        position[3*idx+2] += velocity[3*idx+2]*deltaTime+acceleration[3*idx+2]*deltaTime*deltaTime*0.5;
        velocity[3*idx+0] += acceleration[3*idx+0]*deltaTime*0.5;
        velocity[3*idx+1] += acceleration[3*idx+1]*deltaTime*0.5;
        velocity[3*idx+2] += acceleration[3*idx+2]*deltaTime*0.5;    
    }
    if (position[3*idx+0]<0.0) position[3*idx+0] += box_len; 
    if (position[3*idx+0]>=box_len) position[3*idx+0] -= box_len; 
    if (position[3*idx+1]<0.0) position[3*idx+1] += box_len; 
    if (position[3*idx+1]>=box_len) position[3*idx+1] -= box_len;
    if (position[3*idx+2]<0.0) position[3*idx+2] += box_len; 
    if (position[3*idx+2]>=box_len) position[3*idx+2] -= box_len;
    }

