#pragma once
struct Domain{
    int domainDimension=1;
    int n_particles_x=10;
    int n_particles_y;
    int n_particles_z;
    int n_particles_total;
    double delta_x;
    // const double delta_y=delta_x;
    // const double delta_z=delta_x;
    const double sigma=1.0; //reduced units
    const double epsilon=1.0; //reduced units
    const double rad_cutoff=2.5*sigma; //2.5*sigma
    const float T = 2.68; //reduced temp 
    const double startTime = 0.0;
    const double endTime = 0.5;
    const double deltaTime = 0.005;
    
    int n_cells_x; //same across all axes
    int n_cells_total;
};
