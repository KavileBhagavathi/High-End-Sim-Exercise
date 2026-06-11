#pragma once
#include "particle_props.h"
#include "domain.h"

void launchUpdatePosVelFirstHalfKernel(const Domain& domain, ParticleSystem& device_ps,int blocks, int threads);