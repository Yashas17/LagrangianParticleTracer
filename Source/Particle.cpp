#include "Particle.hpp"

std::array<int, 2>      Particle::getIndex(RealType x, RealType y) {}
std::array<int, 3>      Particle::getIndex(RealType x, RealType y, RealType z) {}
std::array<RealType 2>  Particle::calculateVelocity(int i, int j) {}
std::array<RealType, 3> Particle::calculateVelocity(int i, int j, int k) {}

Particle::Particle(RealType x, RealType y, FlowField& flowField):
  x_(x),
  y_(y),
  flowField_(flowField) {}

Particle::Particle(RealType x, RealType y, RealType z, FlowField& flowField):
  x_(x),
  y_(y),
  z_(z),
  flowField_(flowField) {}

RealType Particle::getX(){return x_};
RealType Particle::getY(){return y_};
RealType Particle::getZ(){return z_};

void Particle::update() {
  // std::array<int,>
}