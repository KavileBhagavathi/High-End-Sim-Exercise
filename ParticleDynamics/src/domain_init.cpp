#include <cmath>
#include <iostream>
#include "domain.h"
void domainParamInitialization(Domain& domain){
    domain.dimension = 3;
    domain.n_particles_x = 10;
    domain.n_particles_y = domain.n_particles_x;
    domain.n_particles_z = domain.n_particles_x;
    domain.n_particles_total = domain.n_particles_x*domain.n_particles_y*domain.n_particles_z;
    domain.delta_x = 1.0;
    domain.delta_y = domain.delta_x;
    domain.delta_z = domain.delta_x;

    domain.sigma = 1.0; //reduced units
    domain.epsilon = 1.0;
    domain.rad_cutoff = 2.5*domain.sigma;
    domain.deltaTime = 0.005;
    domain.endTime = 1.0;
    domain.T = 2.68; //reduced units
    domain.gravityEnabled = true;
    domain.periodicBC = true;
    
    std::cout<<"Number of particles in x-axis: "<<domain.n_particles_x<<"\n";
    std::cout<<"Number of particles in y-axis: "<<domain.n_particles_y<<"\n";
    std::cout<<"Number of particles in z-axis: "<<domain.n_particles_z<<"\n";
    std::cout<<"Total number of particles: "<<domain.n_particles_total<<"\n";
}

void domainBoxInitialization(Domain& domain){
    domain.box_len_x = (double)(domain.n_particles_x)*domain.delta_x;
    domain.box_len_y = (double)(domain.n_particles_y)*domain.delta_y;
    domain.box_len_z = (double)(domain.n_particles_z)*domain.delta_z;

    domain.n_cells_x = (int)std::ceil(domain.box_len_x/domain.rad_cutoff);
    domain.n_cells_y = (int)std::ceil(domain.box_len_y/domain.rad_cutoff);
    domain.n_cells_z = (int)std::ceil(domain.box_len_z/domain.rad_cutoff);

    if (domain.n_cells_x<1) domain.n_cells_x = 1;
    if (domain.n_cells_y<1) domain.n_cells_y = 1;
    if (domain.n_cells_z<1) domain.n_cells_z = 1;

    domain.n_cells_total = domain.n_cells_x*domain.n_cells_y*domain.n_cells_z;
    std::cout << "Grid Dimensions (cells): [" << domain.n_cells_x << "x" 
              << domain.n_cells_y << "x" << domain.n_cells_z << "]\n";
    std::cout << "Total Cells: " << domain.n_cells_total << "\n";
    std::cout << "Box len: "<<domain.box_len_x<<"\n";
}