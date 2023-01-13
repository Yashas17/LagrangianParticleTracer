#include "ParticleSimulation.hpp"

#include <fstream>

ParticleSimulation::ParticleSimulation(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField) {}

void ParticleSimulation::initializeParticles() {
  // Number of particles in each direction (2D: y - 3D: y & z)
  const int particleCount = parameters_.particles.particleCount;
  const int dim           = parameters_.geometry.dim;

  const RealType lengthY = (parameters_.bfStep.yRatio <= 0.0)
                             ? parameters_.geometry.lengthY
                             : parameters_.geometry.lengthY * (1 - parameters_.bfStep.yRatio);
  const RealType lengthZ = parameters_.geometry.lengthZ;

  // Assuming uniform spacing between particles
  double spacingY = lengthY / (particleCount + 1);
  double spacingZ = lengthZ / (particleCount + 1);

  RealType x = parameters_.meshsize->getDx(2, 0) / 2;
  RealType y = (parameters_.bfStep.yRatio >= 0) ? (parameters_.bfStep.yRatio * parameters_.geometry.lengthY) : 0;
  std::array<int, 3> index = {2, 0, 0};

  ASSERTION(dim == 2 || dim == 3);

  if (dim == 2) {
    // Uniform distributution of the particles in the y-direction
    // Avoid having particles at exactly the bottom and top walls
    int j = 0;
    for (int p = 0; p < particleCount; p++) {
      y += spacingY;
      while (parameters_.meshsize->getPosY(0, j) + parameters_.meshsize->getDy(0, j) <= y) {
        j++;
      }
      index[1] = j;
      particles_.push_back(Particle(x, y, index, flowField_, parameters_));
    }
  } else {
    RealType z;

    // Uniform distributution of the particles in the YZ-direction
    // Total of particleCount*particleCount particles in the YZ-plane
    // Avoid having particles at exactly the bottom, top, front, and back walls
    int j = 0;
    int k = 0;
    for (int py = 0; py < particleCount; py++) {
      k = 0;
      y += spacingY;
      z = 0;
      while (parameters_.meshsize->getPosY(0, j, k) + parameters_.meshsize->getDy(0, j, k) <= y) {
        j++;
      }
      index[1] = j;
      for (int pz = 0; pz < particleCount; pz++) {
        z += spacingZ;
        while (parameters_.meshsize->getPosZ(0, 0, k) + parameters_.meshsize->getDz(0, 0, k) <= z) {
          k++;
        }
        index[2] = k;
        particles_.push_back(Particle(x, y, z, index, flowField_, parameters_));
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

void ParticleSimulation::plot(int timeSteps, double time) {
  int i = 1;

  for (auto& particle : particles_) {
    std::cout
      << "Particle ID:" << i << "\t(x,y,z) = " << particle.getX() << "," << particle.getY() << "," << particle.getZ()
      << "\n";
    i++;
  }
}