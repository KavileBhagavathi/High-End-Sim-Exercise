#pragma once
#include <cmath>
#include <vector>

void createNeighbourList(std::vector<int>& cells, std::vector<int>& particles,
                         std::vector<std::vector<double>>& position, Domain& domain ){

    std::vector<int> cells_arr(domain.num_cells_x,-1);
    std::vector<int> particles_arr(domain.n_particles_total,-1);

    for (int i=0; i<particles.size(); ++i){
        int cell_x = (int)(position[i][0]/domain.rad_cutoff);
        int cell_y = (int)(position[i][1]/domain.rad_cutoff);
        int cell_z = (int)(position[i][2]/domain.rad_cutoff);
        int c = cell_x + domain.num_cells_x*cell_y + 
                            domain.num_cells_x*domain.num_cells_x*cell_z;
        particles[i] = cells[c];
        cells[c] = i;
    }
}

void LennardJones(Domain& domain, std::vector<int>& cells, std::vector<int>& particles,
                     std::vector<std::vector<double>>& position, std::vector<Particle>& p_arr){
    createNeighbourList(cells,particles,position,domain);
    for (int l=0; l<cells.size(); ++l){
        std::vector<double> Force(3,0); //initialized Fx, Fy, Fz with zero
        int particle_idx = cells[l];
        while(particle_idx!=-1){
            int next_particle_idx = particles[particle_idx];
            while (next_particle_idx != -1){
                double xij_sqrd = pow(position[particle_idx][0]-position[next_particle_idx][0],2)+
                                    pow(position[particle_idx][1]-position[next_particle_idx][1],2)+
                                    pow(position[particle_idx][2]-position[next_particle_idx][2],2);
                double sigma_xij_powsix = pow(domain.sigma,6)/ pow(xij_sqrd,3);
                std::vector<double> Xij (3,0);
                double multiplier = 24*domain.epsilon*sigma_xij_powsix*(2*sigma_xij_powsix-1)*(1.0/xij_sqrd);
                for (auto k=0; k<3;++k){
                    Xij[k] = position[particle_idx][k]-position[next_particle_idx][k];
                    Xij[k] *= multiplier;                       
                    Force[k] += Xij[k]; 
                    p_arr[particle_idx].p_acceleration[k] = Force[k]*(1.0/p_arr[particle_idx].p_mass);
                }
                next_particle_idx = particles[next_particle_idx];
            }
            particle_idx = particles[particle_idx];
        }    
    }
}