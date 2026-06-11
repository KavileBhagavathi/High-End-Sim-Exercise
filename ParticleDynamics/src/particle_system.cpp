#inlcude

void initializeParticleSystem(ParticleSystem& host_ps, const Domain& domain){
    std::mt19937 generator(21);
    const double stddev = std::sqrt(domain.T);
    std::normal_distribution<double> guassianDistr(0.0,stddev);
}