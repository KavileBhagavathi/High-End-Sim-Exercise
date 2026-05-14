#pragma once
struct Domain{
    int n_particles_x;
    int n_particles_y;
    int n_particles_z;
    int n_particles_total;
    double delta_x;
    double delta_y;
    double delta_z;
    double sigma;
    double epsilon;
    double rad_cutoff;
    int num_cells_x; //same across all axes
    int num_cells_total;
};