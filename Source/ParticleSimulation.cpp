#include "ParticleSimulation.hpp"


ParticleSimulation::ParticleSimulation(Parameters &parameters)
: parameters_(parameters) {}

void ParticleSimulation::initializeParticles() {  
    // Number of particles in each direction (2D: y - 3D: y & z)
    const int particleCount = parameters_.particles.particleCount;
    const int dim = parameters_.geometry.dim;

    const RealType lengthY = parameters_.geometry.lengthY;
    const RealType lengthZ = parameters_.geometry.lengthZ;

    double spacingY = lengthY/(particleCount + 1);
    double spacingZ = lengthZ/(particleCount + 1);
    
    RealType y;

    ASSERTION(dim == 2 || dim == 3);
    
    if (dim == 2) {
        for (int p = 0; p < particleCount; p++) {
            y = spacingY * (p+1);
            particles_.push_back(Particles(0, y, 0));
        }
    }
    else if(dim == 3) {
        RealType z = 0;

        for (int pz = 0; pz < particleCount; pz++) {
            for (int py = 0; py < particleCount; py++) {
                y = spacingY * (py+1);
                z = spacingZ * (pz+1);
                particles_.push_back(Particles(0, y, z));
            }
        }
    }
}


void ParticleSimulation::solveTimeStep() {
    std::list<Particles>::iterator particle;

    for (auto& particle : particles_){
        particle.update();
    }

    // TODO: for loop to handle obstacles
    // TODO: for loop to handle going out of bounds (periodicity or dying out?)
}

