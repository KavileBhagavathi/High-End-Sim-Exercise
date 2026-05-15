#pragma once 
#include <vector>

struct ParticleSystem{

    double *p_mass; //in grams 
    double *p_radius; // in mm
    double *p_position_x,*p_position_y,*p_position_z; // in mm
    double *p_velocity_x,*p_velocity_y,*p_velocity_z; // in mm/s
    double *p_acceleration_x,*p_acceleration_y,*p_acceleration_z; // in mm/s^2
        
};

Particle* generateParticles(const Domain& domain){
    const double velMin = -0.001;
    double velMax = 0.01;
    const double mass = 2.0;
    const double radius = 0.5;
    
    Particle* particle_arr;

    for(int pk=0; pk<domain.n_particles_z; ++pk){
        for(int pj=0; pj<domain.n_particles_y; ++pj){
            for(int pi=0; pi<domain.n_particles_x; ++pi){
                // std::vector<double> vel (3,0);
                double* vel = [0.0,0.0,0.0];
                double* acc = [0.0,0.0,0.0];
                double* pos = [((double)pi*domain.delta_x),
                                         ((double)pj*domain.delta_y),
                                         ((double)pk*domain.delta_z)];
                Particle p {mass,radius,pos,vel,acc};
                particle_arr.push_back(p);
            }
        }
    }
    return particle_arr;
}


void updatePositionVelocityFirstHalf(std::vector<Particle>& particles, double deltaTime){
    int n = particles.size();
    for (int i=0; i<n; ++i){
        particles[i].p_position[0] += particles[i].p_velocity[0]*deltaTime 
                                + particles[i].p_acceleration[0]*deltaTime*deltaTime*0.5;

        particles[i].p_position[1] += particles[i].p_velocity[1]*deltaTime
                                + particles[i].p_acceleration[1]*deltaTime*deltaTime*0.5;

        particles[i].p_position[2] += particles[i].p_velocity[2]*deltaTime
                                + particles[i].p_acceleration[2]*deltaTime*deltaTime*0.5;

        // position[i] = particles[i].p_position;

        particles[i].p_velocity[0] += particles[i].p_acceleration[0]*deltaTime*0.5;
        particles[i].p_velocity[1] += particles[i].p_acceleration[1]*deltaTime*0.5; 
        particles[i].p_velocity[2] += particles[i].p_acceleration[2]*deltaTime*0.5; 
    }
}

void updateVelocitySecondHalf(std::vector<Particle>& particles, double deltaTime){
    int n = particles.size();
    for (int i=0; i<n; ++i){
        particles[i].p_velocity[0] += particles[i].p_acceleration[0]*deltaTime*0.5;
        particles[i].p_velocity[1] += particles[i].p_acceleration[1]*deltaTime*0.5; 
        particles[i].p_velocity[2] += particles[i].p_acceleration[2]*deltaTime*0.5; 
    }  
}