#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include "domain.h"
#include "cudautils.h"
#include "vtuWriter.h"
#include "particle.h"
#include "force.h"

struct GpuLaunchConfig {
    dim3 blocks;
    dim3 threads;
};

// 3D Configuration Function
GpuLaunchConfig getLaunchConfig(const Domain& domain) {
    
    
    GpuLaunchConfig config;
    const int threadsPerAxis = 8;
    // 8 threads per axis = 512 threads per block. Highly efficient for warps.
     
    // Compute independent ceiling division blocks for every single spatial axis
    int bx = (domain.n_particles_x + threadsPerAxis - 1) / threadsPerAxis;
    int by = (domain.n_particles_y + threadsPerAxis - 1) / threadsPerAxis;
    int bz = (domain.n_particles_z + threadsPerAxis - 1) / threadsPerAxis;
    if (domain.domainDimension==1){
        config.threads = dim3(threadsPerAxis,1,1);
        config.blocks = dim3(bx,1,1);
    }
    else if(domain.domainDimension==2){
        config.threads = dim3(threadsPerAxis,threadsPerAxis,1);
        config.blocks = dim3(bx,by,1);
    }
    else{
        config.threads = dim3(threadsPerAxis, threadsPerAxis, threadsPerAxis);
        config.blocks  = dim3(bx, by, bz);
    }
    return config;
}

void domainInitialization(Domain& domain){
    if (domain.domainDimension==1){
        domain.n_particles_total = domain.n_particles_x;
        domain.n_particles_y = domain.n_particles_z = 1;
    }
    else if (domain.domainDimension==2) {
        domain.n_particles_total = domain.n_particles_x*domain.n_particles_x;
        domain.n_particles_y = domain.n_particles_x;
        domain.n_particles_z = 1;

    }
    else{
        domain.n_particles_total = domain.n_particles_x*domain.n_particles_x*domain.n_particles_x;
        domain.n_particles_y = domain.n_particles_z = domain.n_particles_x;
    }
    domain.delta_x = domain.rad_cutoff;
    double box_len = (double)(domain.n_particles_x)*domain.delta_x; //assuming equidistant distribution in all directions
    std::cout<<"Box len: "<<box_len<<"\n";
    domain.n_cells_x = (int)std::ceil(box_len/domain.rad_cutoff);

    if (domain.domainDimension==1) domain.n_cells_total = domain.n_cells_x;
    else if (domain.domainDimension==2) domain.n_cells_total = pow(domain.n_cells_x,2);
    else domain.n_cells_total = pow(domain.n_cells_x,3);
}

int main(){

    Domain domain;
    domainInitialization(domain);


    //host arrays
    double* h_position = new double[3*domain.n_particles_total](); //multiplied by three since particle is in Cartesian space
    double* h_velocity = new double[3*domain.n_particles_total](); //multiplied by three since particle is in Cartesian space  
    double* h_acceleration = new double[3*domain.n_particles_total](); //multiplied by three since particle is in Cartesian space  
    
    //device arrays
    double* d_position = nullptr;
    double* d_velocity = nullptr;   
    double* d_acceleration = nullptr;   
    GpuLaunchConfig config = getLaunchConfig(domain);
    
    // cudaDeviceProp prop;
    // cudaGetDeviceProperties(&prop, 0);
    // std::cout << "Max threads per block: " << prop.maxThreadsPerBlock << "\n";
    // std::cout << "Max threads per SM: " << prop.maxThreadsPerMultiProcessor << "\n";
    // std::cout << "Max grid size X: " << prop.maxGridSize[0] << "\n";

    //allocate device memory and copy to device
    checkCudaError(cudaMalloc(&d_position, 3*domain.n_particles_total*sizeof(double)));
    checkCudaError(cudaMalloc(&d_velocity, 3*domain.n_particles_total*sizeof(double)));
    checkCudaError(cudaMalloc(&d_acceleration, 3*domain.n_particles_total*sizeof(double)));
    
    curandState* d_states;
    checkCudaError(cudaMalloc(&d_states,sizeof(curandState)*domain.n_particles_total));
    
    initCurandStates<<<config.blocks,config.threads>>>(d_states,domain);
    cudaDeviceSynchronize();
    generateParticles<<<config.blocks,config.threads>>>(domain,d_position,d_velocity,d_acceleration);
    cudaDeviceSynchronize();
    
    gaussianVelocityDistribution<<<config.blocks,config.threads>>>(d_states,d_velocity,domain);
    checkCudaError(cudaMemcpy(h_position,d_position,3*domain.n_particles_total*sizeof(double),cudaMemcpyDeviceToHost));
    checkCudaError(cudaGetLastError());
    
    checkCudaError(cudaFree(d_states));
    double p_mass = 1.0;
    double p_radius = 0.5;

    printf("Domain dimension: %d\n",domain.domainDimension);
    printf("Num cells total: %d\n",domain.n_cells_total);
    printf("Num cells x: %d\n",domain.n_cells_x);
    printf("Num of particles in x: %d\n", domain.n_particles_x);
    printf("Num of total particles: %d\n",domain.n_particles_total);
    printf("Num of threads in x: %d, y: %d, z: %d\n",config.threads.x,config.threads.y,config.threads.z);
    printf("Num of blocks in x: %d, y: %d, z: %d\n",config.blocks.x,config.blocks.y,config.blocks.z);

    //initializing host cells and particles array for neighbourlist
    int* h_cells_arr = new int[domain.n_cells_total]();
    int* h_particles_arr = new int[domain.n_particles_total]();

    //initializing device cells and particles array for neighbourlist
    int* d_cells_arr = nullptr;
    int* d_particles_arr = nullptr;
    checkCudaError(cudaMalloc(&d_cells_arr,sizeof(int)*domain.n_cells_total));
    checkCudaError(cudaMalloc(&d_particles_arr,sizeof(int)*domain.n_particles_total));
    checkCudaError(cudaMemcpy(d_cells_arr,h_cells_arr,sizeof(int)*domain.n_cells_total,cudaMemcpyHostToDevice));
    checkCudaError(cudaMemcpy(d_particles_arr,h_particles_arr,sizeof(int)*domain.n_particles_total,cudaMemcpyHostToDevice));
    
    std::ofstream pvd("particles.pvd");
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd << "  <Collection>\n";
    std::string fname = "particle_"+std::to_string(0)+".vtu";
    // writeVTU(fname,h_position,domain.n_particles_total,p_radius);
    int count = 0;
    // solver loop starts
    for (double it=domain.startTime; it<=domain.endTime; it += domain.deltaTime){
        cudaMemset(d_cells_arr, -1, domain.n_cells_total*sizeof(int));
        cudaMemset(d_particles_arr, -1, domain.n_particles_total*sizeof(int));
        updatePositionVelocityFirstHalf<<<config.blocks,config.threads>>>(domain,d_position,d_velocity,d_acceleration);
        cudaDeviceSynchronize();
        checkCudaError(cudaGetLastError());
        
        createNeighbourList<<<config.blocks,config.threads>>>(domain,d_cells_arr,d_particles_arr,d_position);
        cudaDeviceSynchronize();
        checkCudaError(cudaGetLastError());
        
        LennardJones_VelocitySecondHalf<<<config.blocks,config.threads>>>(domain,d_cells_arr,d_particles_arr,d_position,d_velocity,d_acceleration,p_mass);
        cudaDeviceSynchronize();
        checkCudaError(cudaGetLastError());
        
        checkCudaError(cudaMemcpy(h_position,d_position,3*domain.n_particles_total*sizeof(double),cudaMemcpyDeviceToHost));
        
        fname = "particle_"+std::to_string(count)+".vtu";
        writeVTU(fname,h_position,domain.n_particles_total,p_radius);
        pvd<<"    <DataSet timestep=\"" << count << "\" file=\"" << fname << "\"/>\n";
        count++;
    }
    pvd<<"    <DataSet timestep=\"" << count << "\" file=\"" << fname << "\"/>\n";
    
    checkCudaError(cudaFree(d_position));
    checkCudaError(cudaFree(d_velocity));
    checkCudaError(cudaFree(d_acceleration));
    checkCudaError(cudaFree(d_cells_arr));
    checkCudaError(cudaFree(d_particles_arr));
    
    delete[] h_position;
    delete[] h_velocity;
    delete[] h_acceleration;
    delete[] h_cells_arr;
    delete[] h_particles_arr;

    printf("Program exited without errors\n");
    pvd<<"  </Collection>\n";
    pvd<<"</VTKFile>\n";
}
