#pragma once
#include "StdAfx.hpp"
#include "Definitions.hpp"

class Particles {
private:
  RealType x_, y_, z_; // coordinates

  std::array<int, 2>      getIndex(RealType x, RealType y);             // get cell index where the particle is located
  std::array<int, 3>      getIndex(RealType x, RealType y, RealType z); // get cell index where the particle is located
  std::array<RealType 2>  calculateVelocity(int i, int j);              // calculate particle velocity
  std::array<RealType, 3> calculateVelocity(int i, int j, int k);       // calculate particle velocity

public:
  Particles(RealType x, RealType y, RealType z);
  void update();
}