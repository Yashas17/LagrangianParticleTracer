#pragma once

#include "StdAfx.hpp"

#include "Definitions.hpp"
#include "FlowField.hpp"

class Particle {
private:
  RealType           x_, y_, z_; // coordinates
  FlowField&         flowField_;
  Parameters&        parameters_;
  std::array<int, 3> index_; // index of cell particle is in

  RealType* calculateVelocity(); // calculate particle velocity

public:
  Particle(RealType x, RealType y, std::array<int, 3> index, FlowField& flowField_, Parameters& parameters);
  Particle(RealType x, RealType y, RealType z, std::array<int, 3> index, FlowField& flowField_, Parameters& parameters);
  RealType getX();
  RealType getY();
  RealType getZ();

  void update(RealType dt);
};