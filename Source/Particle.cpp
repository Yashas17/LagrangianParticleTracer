#include "Particle.hpp"

std::array<int, 2> Particle::getIndex(RealType x, RealType y) {}
std::array<int, 3> Particle::getIndex(RealType x, RealType y, RealType z) {}

RealType* Particle::calculateVelocity(int i, int j) { return flowField_.getVelocity.getVector(i, j); }

RealType* Particle::calculateVelocity(int i, int j, int k) { return flowField_.getVelocity.getVector(i, j, k); }

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

void Particle::update(RealType dt) {
  if (z_ >= 0) {
    std::array<int, 3> index       = getIndex(x_, y_, z_);
    RealType           velocity[3] = calculateVelocity(index[0], index[1], index[2]);
    x_ += dt * velocity[0];
    y_ += dt * velocity[1];
    z_ += dt * velocity[2];
  } else {
    std::array<int, 2> index       = getIndex(x_, y_);
    RealType           velocity[2] = calculateVelocity(index[0], index[1]);
    x_ += dt * velocity[0];
    y_ += dt * velocity[1];
  }
}