#include <string>
struct Domain{
    int dimension;

    int n_particles_x;
    int n_particles_y;
    int n_particles_z;

    double box_len_x;
    double box_len_y;
    double box_len_z;

    int n_cells_x;
    int n_cells_y;
    int n_cells_z;
    int n_cells_total;

    double delta_x;
    double rad_cutoff;
    double sigma;
    double epsilon;
    double deltaTime;
    double endTime; 
    double T; //temperature

    std::string bc;
};