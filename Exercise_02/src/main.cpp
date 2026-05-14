#include <vector>
#include <random>
#include <iostream>
#include <cmath>
#include "particle.h"
#include "vtuWriter.h"
#include "domain.h"
#include "force.h"

std::vector<double> randVec(const double min, const double max){
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(min, max); // define the range
    //std::cout<<"Random num is: "<<distr(gen)<<"\n";
    std::vector<double> vec;
    for (int i=0; i<3; ++i){
        double randNum = (double)distr(gen)/1.0;
        vec.push_back(randNum);
        // std::cout<<vec[i]<<",";
    }
    return vec;
}

std::vector<std::vector<double>> generateRegularGrid(Domain& domain){
    std::vector<std::vector<double>> init_pos;
    for(int k=0; k<domain.n_particles_z; ++k){
        for(int j=0; j<domain.n_particles_y; ++j){
            for(int i=0; i<domain.n_particles_x; ++i){
                std::vector<double> pos{((double)i*domain.delta_x),
                                         ((double)j*domain.delta_y),
                                         ((double)k*domain.delta_z)};
                init_pos.push_back(pos);
            }
        }
    }
    return init_pos;
}

bool checkOverlap(const std::vector<Particle>& p_arr, 
                    const std::vector<double>& pos, const double radius){
    for (int i=0; i<p_arr.size(); ++i){
        double dist = pow(p_arr[i].p_position[0]-pos[0],2) +
                        pow(p_arr[i].p_position[1]-pos[1],2) +
                        pow(p_arr[i].p_position[2]-pos[2],2);
        if (dist < 4*radius*radius) return false;
    }
    return true;
}

std::vector<Particle> generateParticles(const int N, const std::vector<std::vector<double>>& init_pos){
    const double velMin = -1.0;
    double velMax = 1.0;
    double posMin = 0.0;
    double posMax = 100.0;
    const double mass = 2.0;
    const double radius = 1.0;
    const double delta = 5.0; //rcutoff
    
    bool firstparticle = true;
    std::vector<Particle> particle_arr;
    int iter = 0;
    while (iter<N){
        std::vector<double> vel;
        if (iter==0){
            vel = std::vector<double>(3,1);
        } 
        else {
            vel = std::vector<double>(3,0);
        }
        std::vector<double> pos = init_pos[iter];
        std::vector<double> acc (3,0);

        Particle p {mass,radius,pos,vel,acc};
        // if (firstparticle){
        //     particle_arr.push_back(p);
        //     firstparticle = false;
        //     continue;
        // }
        // if (checkOverlap(particle_arr,pos,radius)){
        //     particle_arr.push_back(p);
        //     active_particles ++;
        // }
        particle_arr.push_back(p);
        iter++;
    }
    std::cout<<"iter: "<<iter<<"\n";
    return particle_arr;
}





int main(){
    //initialize the pvd writer
    std::ofstream pvd("particles.pvd");
    pvd << "<?xml version=\"1.0\"?>\n";
    pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    pvd << "  <Collection>\n";
    //
    Domain domain;
    domain.n_particles_x = 10; //keep it min 10 
    domain.n_particles_y = 10;
    domain.n_particles_z = 10;
    domain.delta_x = 3.0;
    domain.delta_y = 3.0;    
    domain.delta_z = 3.0; 
    domain.sigma = 1.0; //subjected to change  
    domain.rad_cutoff = 2.5 * domain.sigma;
    domain.n_particles_total = domain.n_particles_x * domain.n_particles_y * domain.n_particles_z;
    double box_len = domain.n_particles_x*domain.delta_x; //assuming equidistant distribution in all directions
    
    std::vector<std::vector<double>> position = generateRegularGrid(domain);
    std::vector<Particle> p_arr = generateParticles(domain.n_particles_total,position);  
    double startTime = 0.0;
    double endTime = 20.0;
    double deltaTime = 1.0;
    int it = 0;
    double radius = p_arr[0].p_radius;
    domain.num_cells_x = (int)(box_len/domain.rad_cutoff);
    domain.num_cells_total = pow(domain.num_cells_x,3);

    std::cout<<"Num cells in x: "<<domain.num_cells_x<<"\n";

    std::vector<int> cells_arr(domain.num_cells_x,-1);
    std::vector<int> particles_arr(domain.n_particles_total,-1);

    std::string fname = "particle_"+std::to_string(it)+".vtu";
    writeVTU(fname,position,domain.n_particles_total,radius);
    while (it<=endTime){
        updatePositionVelocityFirstHalf(p_arr,deltaTime,position);
        // LennardJones(domain,cells_arr,particles_arr,position,p_arr);
        LennardJones(domain,position,p_arr);
        updateVelocitySecondHalf(p_arr,deltaTime,position);
        fname = "particle_"+std::to_string(it)+".vtu";
        writeVTU(fname,position,domain.n_particles_total,radius);
        //writePlane("plane.vtu");
        pvd<<"    <DataSet timestep=\"" << it << "\" file=\"" << fname << "\"/>\n";
        //pvd << "    <DataSet timestep=\"" << it << "\" file=\"plane.vtu\"/>\n";
        position.erase(position.begin(),position.end());
        //printf("Written %d\n",it);
        it += deltaTime;
    }
    pvd<<"  </Collection>\n";
    pvd<<"</VTKFile>\n";
    // std::vector<double> velocity = randNumGen(-100,100);
}