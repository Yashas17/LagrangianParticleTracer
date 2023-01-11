#include "ParticleSimulation.hpp"


ParticleSimulation::ParticleSimulation(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField) {}

void ParticleSimulation::initializeParticles() {  
    // Number of particles in each direction (2D: y - 3D: y & z)
    const int particleCount = parameters_.particles.particleCount;
    const int dim = parameters_.geometry.dim;

    const RealType lengthY = parameters_.geometry.lengthY;
    const RealType lengthZ = parameters_.geometry.lengthZ;

    // Assuming uniform spacing between particles
    double spacingY = lengthY/(particleCount + 1);
    double spacingZ = lengthZ/(particleCount + 1);
    
    RealType y;

    ASSERTION(dim == 2 || dim == 3);
    
    if (dim == 2) {
        // Uniform distributution of the particles in the y-direction
        // Avoid having particles at exactly the bottom and top walls
        for (int p = 0; p < particleCount; p++) {
            y = spacingY * (p+1);
            particles_.push_back(Particles(0.0, y));
        }
    }
    else {
        RealType z;

        // Uniform distributution of the particles in the YZ-direction
        // Total of particleCount*particleCount particles in the YZ-plane
        // Avoid having particles at exactly the bottom, top, front, and back walls
        for (int pz = 0; pz < particleCount; pz++) {
            for (int py = 0; py < particleCount; py++) {
                y = spacingY * (py+1);
                z = spacingZ * (pz+1);
                particles_.push_back(Particles(0.0, y, z));
            }
        }
    }
  
}

void ParticleSimulation::solveTimeStep() {

    for (auto& particle : particles_){
        particle.update();
    }

  // TODO: for loop to handle obstacles
  // TODO: for loop to handle going out of bounds (periodicity or dying out?)
}
