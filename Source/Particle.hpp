#pragma once

#include "StdAfx.hpp"

#include "Definitions.hpp"
#include "FlowField.hpp"

class Particle {
private:
  RealType   x_, y_, z_ = -1; // coordinates
  FlowField& flowField_;

  std::array<int, 2>      getIndex(RealType x, RealType y);             // get cell index where the particle is located
  std::array<int, 3>      getIndex(RealType x, RealType y, RealType z); // get cell index where the particle is located
  std::array<RealType 2>  calculateVelocity(int i, int j);              // calculate particle velocity
  std::array<RealType, 3> calculateVelocity(int i, int j, int k);       // calculate particle velocity

public:
  Particle(RealType x, RealType y, FlowField& flowField_);
  Particle(RealType x, RealType y, RealType z, FlowField& flowField_);
  RealType getX();
  RealType getY();
  RealType getZ();
  
  void update();
}