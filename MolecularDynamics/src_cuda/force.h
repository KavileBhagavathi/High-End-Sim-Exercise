#pragma once
#include <cuda_runtime.h>

__global__ void createNeighbourList(const Domain domain, 
    int* __restrict__ cells_arr, int* __restrict__ neigh_particle_arr,                                
    const double* __restrict__ position
                                    ){
    //make sure to reset the cells and particle arrays
    //if (position[3*ti+0]<0 || position[3*tj+1]<0 || position[3*tk+2]<0) continue;

    const int ti = blockDim.x*blockIdx.x+threadIdx.x;
    const int tj = blockDim.y*blockIdx.y+threadIdx.y;
    const int tk = blockDim.z*blockIdx.z+threadIdx.z;
    if (ti>=domain.n_particles_x || tj>=domain.n_particles_y || tk>=domain.n_particles_z) return;

    const int idx = ti + tj*domain.n_particles_x + tk*domain.n_particles_x*domain.n_particles_y;

    int cell_x = (int)(position[3*idx+0]/domain.rad_cutoff);
    int cell_y = (int)(position[3*idx+1]/domain.rad_cutoff);
    int cell_z = (int)(position[3*idx+2]/domain.rad_cutoff);

    if (cell_x < domain.n_cells_x && cell_y < domain.n_cells_x && cell_z < domain.n_cells_x
        && cell_x >= 0 && cell_y >= 0 && cell_z >=0){
        int c = cell_x + domain.n_cells_x*cell_y + domain.n_cells_x*domain.n_cells_x*cell_z;
        int neigh_particle_head = atomicExch(&cells_arr[c],idx);
        neigh_particle_arr[idx] = neigh_particle_head;
    }

}

__global__ void LennardJones_VelocitySecondHalf(const Domain domain, 
                            const int* __restrict__ cells_arr,
                            const int* __restrict__ neigh_particle_arr,
                            double* __restrict__ position,
                            double* __restrict__ velocity,
                            double* __restrict__ acceleration, const double mass){
    
    const int ti = blockDim.x*blockIdx.x+threadIdx.x;
    const int tj = blockDim.y*blockIdx.y+threadIdx.y;
    const int tk = blockDim.z*blockIdx.z+threadIdx.z;
    if (ti>=domain.n_particles_x || tj>=domain.n_particles_y || tk>=domain.n_particles_z) return;
    int idx = ti + tj*domain.n_particles_x + tk*domain.n_particles_x*domain.n_particles_y;     
          
    acceleration[3*idx+0] = 0.0;
    acceleration[3*idx+1] = 0.0;
    acceleration[3*idx+2] = 0.0;
    double Fx = 0.0;
    double Fy = 0.0;
    double Fz = 0.0;
    
    double sigma_pow_six = domain.sigma*domain.sigma*domain.sigma*domain.sigma*domain.sigma*domain.sigma;

    double px = position[3*idx+0];
    double py = position[3*idx+1];
    double pz = position[3*idx+2];

    int i = (int)(px/domain.rad_cutoff);
    int j = (domain.domainDimension==2)?(int)(pz/domain.rad_cutoff):0;
    int k = (domain.domainDimension==3)?(int)(pz/domain.rad_cutoff):0;

    int min_dx = -1, max_dx = 1;
    int min_dy = (domain.domainDimension==2)?-1:0;
    int max_dy = (domain.domainDimension==2)? 1:0;
    int min_dz = (domain.domainDimension==3)?-1:0;
    int max_dz = (domain.domainDimension==3)? 1:0;

    // for(int l=0; l<domain.n_cells_total; ++l){
    //     int particle_idx = cells_arr[l];
    //     if (particle_idx==idx){
    //         if (particle_idx != -1){
    //             double Fx = 0.0;
    //             double Fy = 0.0;
    //             double Fz = 0.0;
    //             int i = (int)(l%domain.n_cells_x);
    //             int j = (int)(l/domain.n_cells_x)%domain.n_cells_x;
    //             int k = (int)(l/(domain.n_cells_x*domain.n_cells_x));

    for (int dz=min_dz; dz<=max_dz; ++dz){
        for (int dy=min_dy; dy<=max_dy; ++dy){
            for (int dx=min_dx; dx<=max_dx; ++dx){
                int neigh_i = i+dx;
                int neigh_j = j+dy;
                int neigh_k = k+dz;
                if (neigh_i<0 || neigh_i >= domain.n_cells_x) continue;
                if (neigh_j<0 || neigh_j >= domain.n_cells_x) continue;
                if (neigh_k<0 || neigh_k >= domain.n_cells_x) continue;
                int neigh_cell_idx = neigh_i + domain.n_cells_x*neigh_j + domain.n_cells_x*domain.n_cells_x*neigh_k;
                int neigh_cell_particle_idx = cells_arr[neigh_cell_idx];

                while(neigh_cell_particle_idx != -1){
                    if (neigh_cell_particle_idx != idx){
                        double dx_pos = px - position[3*neigh_cell_particle_idx+0];
                        double dy_pos = py - position[3*neigh_cell_particle_idx+1];
                        double dz_pos = pz - position[3*neigh_cell_particle_idx+2];

                        double xij_sqrd = dx_pos*dx_pos + dy_pos*dy_pos + dz_pos*dz_pos;

                        if (xij_sqrd<domain.rad_cutoff*domain.rad_cutoff){
                            double sigma_xij_powsix = sigma_pow_six/(xij_sqrd*xij_sqrd*xij_sqrd);
                            double multiplier = 24.0*domain.epsilon*sigma_xij_powsix*(2.0*sigma_xij_powsix-1)*(1.0/xij_sqrd);
                            Fx += multiplier*dx_pos;
                            if (domain.domainDimension == 2) Fy += multiplier*dy_pos;
                            if (domain.domainDimension == 3) Fz += multiplier*dz_pos;
                        }
                    }
                    neigh_cell_particle_idx = neigh_particle_arr[neigh_cell_particle_idx];
                }
            }
        }
    }
    acceleration[3*idx+0] = Fx/mass;
    acceleration[3*idx+1] = Fy/mass;
    acceleration[3*idx+2] = Fz/mass;
    velocity[3*idx+0] += acceleration[3*idx+0]*domain.deltaTime*0.5;
    velocity[3*idx+1] += acceleration[3*idx+1]*domain.deltaTime*0.5;
    velocity[3*idx+2] += acceleration[3*idx+2]*domain.deltaTime*0.5;
}
            
    //     }
        
    // } 
    // }
