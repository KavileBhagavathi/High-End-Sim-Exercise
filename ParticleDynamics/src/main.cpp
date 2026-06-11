#include "domain.h"
#include "particle_system.h"
#include "utils.h"

int main(int argc, char* argv[]){
    Domain domain;
    domainParamInitialization(domain);
    domainBoxInitialization(domain);

    ParticleSystem host_ps;
    ParticleSystem dev_ps;
    allocateHostMemory(domain, host_ps);
    allocateDeviceMemory(domain, dev_ps);

    initializeParticleSystem(domain,host_ps);

    int count = 0;
    outputFile(count,host_ps.pos,domain.n_particles_total,host_ps.radius);

    freeHostMemory(host_ps);
    freeDeviceMemory(dev_ps);

}