#include "kernels.cuh"
#include <cuda_runtime.h>

__global__ void updatePosVelFirstHalfKernel(const Domain domain, double* pos,
                                                                    double* vel,
                                                                    double* acc){
    int idx = blockDim.x*blockIdx.x + threadIdx.x;
    if (idx >= domain.n_particles_total) return;

    double dt = domain.deltaTime;
    pos[3*idx+0] += vel[3*idx+0]*dt + acc[3*idx+0]*dt*dt*0.5;
    pos[3*idx+1] += vel[3*idx+1]*dt + acc[3*idx+1]*dt*dt*0.5;
    pos[3*idx+2] += vel[3*idx+2]*dt + acc[3*idx+2]*dt*dt*0.5;
    vel[3*idx+0] += acc[3*idx+0]*dt*0.5;
    vel[3*idx+1] += acc[3*idx+1]*dt*0.5;
    vel[3*idx+2] += acc[3*idx+2]*dt*0.5;  

    if (domain.bc == "periodic"){
            if (pos[3*idx+0]<0.0) pos[3*idx+0] += box_len_x; 
            if (pos[3*idx+0]>=box_len) pos[3*idx+0] -= box_len_x; 
            if (pos[3*idx+1]<0.0) pos[3*idx+1] += box_len_y; 
            if (pos[3*idx+1]>=box_len) pos[3*idx+1] -= box_len_y;
            if (pos[3*idx+2]<0.0) pos[3*idx+2] += box_len_z; 
            if (pos[3*idx+2]>=box_len) pos[3*idx+2] -= box_len_z;
    }
}

void launchUpdatePosVelFirstHalfKernel(const Domain& domain, ParticleSystem& device_ps,int blocks, int threads){
    updatePosVelFirstHalfKernel<<<blocks,threads>>>(domain,device_ps.pos,device_ps.vel,device_ps.acc);
}