#pragma once

#include "StdAfx.hpp"

#include "FlowField.hpp"
#include "Parameters.hpp"
#include "Particle.hpp"

class ParticleSimulation {
private:
  Parameters&          parameters_;
  FlowField&           flowField_;
  std::list<Particle> particles_;

public:
  ParticleSimulation(Parameters& parameters, FlowField& flowField);

  void initializeParticles();
  void solveTimestep();
  void plot(int timeSteps, RealType time);
};