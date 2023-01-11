#include "Particle.hpp"

RealType* Particle::calculateVelocity() {
  if (parameters_.geometry.dim == 2) {
    return flowField_.getVelocity().getVector(index_[0], index_[1]);
  }
  return flowField_.getVelocity().getVector(index_[0], index_[1], index_[2]);
}

Particle::Particle(RealType x, RealType y, std::array<int, 3> index, FlowField& flowField, Parameters& parameters):
  x_(x),
  y_(y),
  index_(index),
  flowField_(flowField),
  parameters_(parameters) {}

Particle::Particle(
  RealType x, RealType y, RealType z, std::array<int, 3> index, FlowField& flowField, Parameters& parameters
):
  x_(x),
  y_(y),
  z_(z),
  index_(index),
  flowField_(flowField),
  parameters_(parameters) {}

RealType Particle::getX() { return x_; }
RealType Particle::getY() { return y_; }
RealType Particle::getZ() { return z_; }

void Particle::update(RealType dt) {
  if (parameters_.geometry.dim == 3) {

    RealType* velocity = calculateVelocity();
    x_ += dt * velocity[0];
    y_ += dt * velocity[1];
    z_ += dt * velocity[2];

    // Updating x index
    if (velocity[0] >= 0.0) {
      while (parameters_.meshsize->getPosX(index_[0], 0, 0) + parameters_.meshsize->getDx(index_[0], 0, 0) <= x_) {
        index_[0]++;
      }
    } else {
      while (parameters_.meshsize->getPosX(index_[0], 0, 0) >= x_) {
        index_[0]--;
      }
    }

    // Updating y index
    if (velocity[1] >= 0.0) {
      while (parameters_.meshsize->getPosY(0, index_[1], 0) + parameters_.meshsize->getDy(0, index_[1], 0) <= y_) {
        index_[1]++;
      }
    } else {
      while (parameters_.meshsize->getPosY(0, index_[1], 0) >= y_) {
        index_[1]--;
      }
    }

    // Updating z index
    if (velocity[2] >= 0.0) {
      while (parameters_.meshsize->getPosZ(0, 0, index_[2]) + parameters_.meshsize->getDz(0, 0, index_[2]) <= z_) {
        index_[2]++;
      }
    } else {
      while (parameters_.meshsize->getPosZ(0, 0, index_[2]) >= z_) {
        index_[2]--;
      }
    }

  } else {
    RealType* velocity = calculateVelocity();
    x_ += dt * velocity[0];
    y_ += dt * velocity[1];

    // Updating x index
    if (velocity[0] >= 0.0) {
      while (parameters_.meshsize->getPosX(index_[0], 0) + parameters_.meshsize->getDx(index_[0], 0) <= x_) {
        index_[0]++;
      }
    } else {
      while (parameters_.meshsize->getPosX(index_[0], 0) >= x_) {
        index_[0]--;
      }
    }

    // Updating y index
    if (velocity[1] >= 0.0) {
      while (parameters_.meshsize->getPosY(0, index_[1]) + parameters_.meshsize->getDy(0, index_[1]) <= y_) {
        index_[1]++;
      }
    } else {
      while (parameters_.meshsize->getPosY(0, index_[1]) >= y_) {
        index_[1]--;
      }
    }
  }
}