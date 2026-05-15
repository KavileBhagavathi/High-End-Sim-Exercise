#include <cstdio>
#include <iostream>
#include <fstream>
#include <thrust/device_vector.h>
#include "domain.h"
#include "cudautils.h"
#include "vtuWriter.h"
__global__ void generateParticles(Domain domain, double* __restrict__ position, 
                                                 double* __restrict__ velocity,
                                                 double* __restrict__ acceleration){
    int ti = blockDim.x*blockIdx.x+threadIdx.x;
    int tj = blockDim.y*blockIdx.y+threadIdx.y;
    int tk = blockDim.z*blockIdx.z+threadIdx.z;

    // printf("%d,%d,%d\n",ti,tj,tk);
    
    if (ti<domain.n_particles_x && tj<domain.n_particles_x && tk<domain.n_particles_x){
        int idx = ti + domain.n_particles_x*tj + 
                        domain.n_particles_x*domain.n_particles_x*tk;
        position[3*idx+0] = ti*domain.delta_x;
        position[3*idx+1] = tj*domain.delta_x; 
        position[3*idx+2] = tk*domain.delta_x; 

        velocity[3*idx+0] = (idx==0)?1:0;
        velocity[3*idx+1] = (idx==0)?1:0; 
        velocity[3*idx+2] = (idx==0)?1:0; 

        acceleration[3*idx+0] = 0.0;
        acceleration[3*idx+1] = 0.0; 
        acceleration[3*idx+2] = 0.0; 
    }
}


int main(){
    std::ofstream pvd("particles.pvd");
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd << "  <Collection>\n";

    Domain domain;
    
    double box_len = (double)domain.n_particles_x*domain.delta_x; //assuming equidistant distribution in all directions
    std::cout<<"Box len: "<<box_len<<"\n";

    //host arrays
    double* h_position = new double[3*domain.n_particles_total];
    double* h_velocity = new double[3*domain.n_particles_total];   
    double* h_acceleration = new double[3*domain.n_particles_total];   

    //device arrays
    double* d_position = new double[3*domain.n_particles_total];
    double* d_velocity = new double[3*domain.n_particles_total];   
    double* d_acceleration = new double[3*domain.n_particles_total];   
    int baseDim = 3;
    dim3 threads (baseDim,baseDim,baseDim);
    dim3 blocks (1,1,1);
    //allocate device memory and copy to device
    checkCudaError(cudaMalloc(&d_position, 3*domain.n_particles_total*sizeof(double)));
    checkCudaError(cudaMalloc(&d_velocity, 3*domain.n_particles_total*sizeof(double)));
    checkCudaError(cudaMalloc(&d_acceleration, 3*domain.n_particles_total*sizeof(double)));
    checkCudaError(cudaMemcpy(d_position,h_position,3*domain.n_particles_total*sizeof(double),cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_velocity,h_velocity,3*domain.n_particles_total*sizeof(double),cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_acceleration,h_acceleration,3*domain.n_particles_total*sizeof(double),cudaMemcpyHostToDevice));
    generateParticles<<<blocks,threads>>>(domain,d_position,d_velocity,d_acceleration);
    cudaDeviceSynchronize();
    checkCudaError(cudaMemcpy(h_position,d_position,3*domain.n_particles_total*sizeof(double),cudaMemcpyDeviceToHost));
    // double startTime = 0.0;
    // double endTime = 2.0;
    // double deltaTime = 0.01;
    // double radius = particles[0].p_radius;
    // domain.num_cells_x = (int)std::ceil(box_len/domain.rad_cutoff);
    // domain.num_cells_total = pow(domain.num_cells_x,3);
    std::string fname = "particle_"+std::to_string(0)+".vtu";
    writeVTU(fname,h_position,domain.n_particles_total,1);
    checkCudaError(cudaFree(d_position));
    checkCudaError(cudaFree(d_velocity));
    checkCudaError(cudaFree(d_acceleration));
    delete h_position;
    delete h_velocity;
    delete h_acceleration;
    printf("Program exited without errors\n");
}