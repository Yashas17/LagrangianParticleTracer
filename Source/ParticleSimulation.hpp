#pragma once

#include "Particles.hpp"
#include "Parameters.hpp"

class ParticleSimulation{
    private:
    Parameters& parameters_;
    std::list<Particles> particles_;

    public:
    ParticleSimulation(Parameters &parameters);

    void initializeParticles();
    void solveTimeStep();
    void plot(int timeSteps, RealType time);
};