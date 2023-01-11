#include "ParticleSimulation.hpp"

ParticleSimulation::ParticleSimulation(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField) {}

void ParticleSimulation::initializeParticles() {
  // Number of particles in each direction (2D: y - 3D: y & z)
  const int particleCount = parameters_.particles.particleCount;
  const int dim           = parameters_.geometry.dim;

  const RealType lengthY = parameters_.geometry.lengthY;
  const RealType lengthZ = parameters_.geometry.lengthZ;

  // Assuming uniform spacing between particles
  double spacingY = lengthY / (particleCount + 1);
  double spacingZ = lengthZ / (particleCount + 1);

  RealType           y;
  std::array<int, 3> index = {0, 0, 0};

  ASSERTION(dim == 2 || dim == 3);

  if (dim == 2) {
    // Uniform distributution of the particles in the y-direction
    // Avoid having particles at exactly the bottom and top walls
    int j = 0;
    for (int p = 0; p < particleCount; p++) {
      y = spacingY * (p + 1);
      while (parameters_.meshsize->getPosY(0, j) + parameters_.meshsize->getDy(0, j) <= y) {
        j++;
      }
      index[1] = j;
      particles_.push_back(Particle(0.0, y, index, flowField_, parameters_));
    }
  } else {
    RealType z;

    // Uniform distributution of the particles in the YZ-direction
    // Total of particleCount*particleCount particles in the YZ-plane
    // Avoid having particles at exactly the bottom, top, front, and back walls
    int k = 0;
    for (int pz = 0; pz < particleCount; pz++) {
      int j = 0;
      z     = spacingZ * (pz + 1);
      while (parameters_.meshsize->getPosZ(0, 0, k) + parameters_.meshsize->getDz(0, 0, k) <= z) {
        k++;
      }
      index[2] = k;
      for (int py = 0; py < particleCount; py++) {
        y = spacingY * (py + 1);
        while (parameters_.meshsize->getPosY(0, j, k) + parameters_.meshsize->getDy(0, j, k) <= y) {
          j++;
        }
        index[1] = j;
        particles_.push_back(Particle(0.0, y, z, index, flowField_, parameters_));
      }
    }
  }
}

void ParticleSimulation::solveTimestep() {

  for (auto& particle : particles_) {
    particle.update(parameters_.timestep.dt);
  }
  // TODO: for loop to handle obstacles
  // TODO: for loop to handle going out of bounds (periodicity or dying out?)
}

void ParticleSimulation::plot(int timeSteps, double time){}