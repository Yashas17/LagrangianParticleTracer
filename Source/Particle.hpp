#pragma once

#include "StdAfx.hpp"

#include "Definitions.hpp"
#include "FlowField.hpp"

class Particle {
private:
  RealType                x_, y_, z_; // coordinates
  std::array<RealType, 3> velocity_;
  FlowField&              flowField_;
  Parameters&             parameters_;
  std::array<int, 3>      index_; // index of cell particle is in

  void calculateVelocity();                                           // calculate particle velocity
  void applyBoundaryCondition();                                      // apply boundary condition on particle
  inline void wallCorrect(RealType& x, RealType xLimit, RealType& velocity); // helper function for boundary conditions

public:
  Particle(RealType x, RealType y, std::array<int, 3> index, FlowField& flowField_, Parameters& parameters);
  Particle(RealType x, RealType y, RealType z, std::array<int, 3> index, FlowField& flowField_, Parameters& parameters);
  RealType getX();
  RealType getY();
  RealType getZ();

  void update(RealType dt);
};