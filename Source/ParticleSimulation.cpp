#include "ParticleSimulation.hpp"


ParticleSimulation::ParticleSimulation(Parameters &parameters)
: parameters_(parameters) {}

void ParticleSimulation::initializeParticles() {     
    const int particleCount = parameters_.particles.particleCount;
    const int dim = parameters_.geometry.dim;

    const RealType lengthY = parameters_.geometry.lengthY;
    const RealType lengthZ = parameters_.geometry.lengthZ;

    double spacingY = lengthY/(particleCount + 1);
    double spacingZ = lengthZ/(particleCount + 1);
    
    RealType y;
    RealType z = 0;

    for (int p = 0; p < parameters_.particles.particleCount; p++) {
        y = spacingY * (p+1);
        if (dim == 3)
            z = spacingZ * (p+1);
        particles_.push_back(Particles(0, y, z));
    }
}


void ParticleSimulation::solveTimeStep() {
    std::list<Particles>::iterator particle;

    for (particle = particles_.begin(); particle != particles_.end(); ++particle){
        particle->update();
    }

    // TODO: for loop to handle obstacles
    // TODO: for loop to handle going out of bounds (periodicity or dying out?)
}

