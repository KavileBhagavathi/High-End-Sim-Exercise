#pragma once 
#include <vector>

class Particle{
    public:
        double p_mass; //in grams 
        double p_radius; // in mm
        std::vector<double> p_position; // in mm
        std::vector<double> p_velocity; // in mm/s
        std::vector<double> p_acceleration; // in mm/s^2

        Particle(const double& mass, const double& radius,
                 const std::vector<double>& position,     
                 const std::vector<double>& velocity,
                 const std::vector<double>& acceleration) : p_mass{mass},p_radius{radius},p_position{position},p_velocity{velocity}, p_acceleration{acceleration}{}
        
};

void updatePositionVelocityFirstHalf(std::vector<Particle>& p_arr, double deltaTime,std::vector<std::vector<double>>& pos_major){
    int n = p_arr.size();
    for (int i=0; i<n; ++i){
        p_arr[i].p_position[0] += p_arr[i].p_velocity[0]*deltaTime 
                                + p_arr[i].p_acceleration[0]*deltaTime*deltaTime*0.5;

        p_arr[i].p_position[1] += p_arr[i].p_velocity[1]*deltaTime;
                                + p_arr[i].p_acceleration[1]*deltaTime*deltaTime*0.5;

        p_arr[i].p_position[2] += p_arr[i].p_velocity[2]*deltaTime;
                                + p_arr[i].p_acceleration[2]*deltaTime*deltaTime*0.5;

        pos_major.push_back(p_arr[i].p_position);

        p_arr[i].p_velocity[0] += p_arr[i].p_acceleration[0]*deltaTime*0.5;
        p_arr[i].p_velocity[1] += p_arr[i].p_acceleration[1]*deltaTime*0.5; 
        p_arr[i].p_velocity[2] += p_arr[i].p_acceleration[2]*deltaTime*0.5; 
    }
}

void updateVelocitySecondHalf(std::vector<Particle>& p_arr, double deltaTime, std::vector<std::vector<double>>& pos_major){
    int n = p_arr.size();
    for (int i=0; i<n; ++i){
        p_arr[i].p_velocity[0] += p_arr[i].p_acceleration[0]*deltaTime*0.5;
        p_arr[i].p_velocity[1] += p_arr[i].p_acceleration[1]*deltaTime*0.5; 
        p_arr[i].p_velocity[2] += p_arr[i].p_acceleration[2]*deltaTime*0.5; 
    }  
}