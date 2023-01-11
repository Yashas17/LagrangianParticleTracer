#include "ParticleSimulation.hpp"

ParticleSimulation::ParticleSimulation(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField) {}

void ParticleSimulation::initializeParticles() {
  const int particleCount = parameters_.particles.particleCount;
  const int dim           = parameters_.geometry.dim;

  const RealType lengthY = parameters_.geometry.lengthY;
  const RealType lengthZ = parameters_.geometry.lengthZ;

  double spacingY = lengthY / (particleCount + 1);
  double spacingZ = lengthZ / (particleCount + 1);

  RealType y;
  RealType z = 0.0;

  for (int p = 0; p < parameters_.particles.particleCount; p++) {
    y = spacingY * (p + 1);
    if (dim == 3) {
      z = spacingZ * (p + 1);
      particles_.push_back(Particle::Particle(0.0, y, z, flowField_));
    } else {
      particles_.push_back(Particle::Particle(0.0, y, flowField_))
    }
  }
}

void ParticleSimulation::solveTimeStep() {

  for (auto& particle : particles_) {
    particle.update();
  }

  // TODO: for loop to handle obstacles
  // TODO: for loop to handle going out of bounds (periodicity or dying out?)
}
