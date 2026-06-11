#include "domain.h"
#include "particle_system.h"


int main(int argc, char* argv[]){
    Domain domain;
    domainParamInitialization(domain);
    domainBoxInitialization(domain);

    ParticleSystem host_ps;
    ParticleSystem dev_ps;
    allocateHostMemory(domain, host_ps);
    allocateDeviceMemory(domain, dev_ps);

    freeHostMemory(host_ps);
    freeDeviceMemory(dev_ps);
}