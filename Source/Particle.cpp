#include "Particle.hpp"

void Particle::calculateVelocity() {
  if (parameters_.geometry.dim == 2) {
    RealType* velocity1 = flowField_.getVelocity().getVector(index_[0], index_[1]);     //(i,j)
    RealType* velocity2 = flowField_.getVelocity().getVector(index_[0] - 1, index_[1]); //(i-1,j)
    RealType* velocity3 = flowField_.getVelocity().getVector(index_[0], index_[1] - 1); //(i,j-1)

    RealType dx   = parameters_.meshsize->getDx(index_[0], index_[1]);
    RealType dy   = parameters_.meshsize->getDy(index_[0], index_[1]);
    RealType posX = parameters_.meshsize->getPosX(index_[0], index_[1]);
    RealType posY = parameters_.meshsize->getPosY(index_[0], index_[1]);

    velocity_[0] = velocity1[0] * (x_ - posX) / dx + velocity2[0] * (posX + dx - x_) / dx;
    velocity_[1] = velocity1[1] * (y_ - posY) / dy + velocity3[1] * (posY + dy - y_) / dy;
  } else {
    RealType* velocity1 = flowField_.getVelocity().getVector(index_[0], index_[1], index_[2]);     //(i,j,k)
    RealType* velocity2 = flowField_.getVelocity().getVector(index_[0] - 1, index_[1], index_[2]); //(i-1,j,k)
    RealType* velocity3 = flowField_.getVelocity().getVector(index_[0], index_[1] - 1, index_[2]); //(i,j-1,k)
    RealType* velocity4 = flowField_.getVelocity().getVector(index_[0], index_[1], index_[2] - 1); //(i,j,k-1)

    RealType dx = parameters_.meshsize->getDx(index_[0], index_[1], index_[2]);
    RealType dy = parameters_.meshsize->getDy(index_[0], index_[1], index_[2]);
    RealType dz = parameters_.meshsize->getDz(index_[0], index_[1], index_[2]);

    RealType posX = parameters_.meshsize->getPosX(index_[0], index_[1], index_[2]);
    RealType posY = parameters_.meshsize->getPosY(index_[0], index_[1], index_[2]);
    RealType posZ = parameters_.meshsize->getPosZ(index_[0], index_[1], index_[2]);

    velocity_[0] = velocity1[0] * (x_ - posX) / dx + velocity2[0] * (posX + dx - x_) / dx;
    velocity_[1] = velocity1[1] * (y_ - posY) / dy + velocity3[1] * (posY + dy - y_) / dy;
    velocity_[2] = velocity1[2] * (z_ - posZ) / dz + velocity4[2] * (posZ + dz - z_) / dz;
  }
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

    calculateVelocity();
    x_ += dt * velocity_[0];
    y_ += dt * velocity_[1];
    z_ += dt * velocity_[2];

    // Updating x index
    if (velocity_[0] >= 0.0) {
      while (parameters_.meshsize->getPosX(index_[0], 0, 0) + parameters_.meshsize->getDx(index_[0], 0, 0) <= x_) {
        index_[0]++;
      }
    } else {
      while (parameters_.meshsize->getPosX(index_[0], 0, 0) >= x_) {
        index_[0]--;
      }
    }

    // Updating y index
    if (velocity_[1] >= 0.0) {
      while (parameters_.meshsize->getPosY(0, index_[1], 0) + parameters_.meshsize->getDy(0, index_[1], 0) <= y_) {
        index_[1]++;
      }
    } else {
      while (parameters_.meshsize->getPosY(0, index_[1], 0) >= y_) {
        index_[1]--;
      }
    }

    // Updating z index
    if (velocity_[2] >= 0.0) {
      while (parameters_.meshsize->getPosZ(0, 0, index_[2]) + parameters_.meshsize->getDz(0, 0, index_[2]) <= z_) {
        index_[2]++;
      }
    } else {
      while (parameters_.meshsize->getPosZ(0, 0, index_[2]) >= z_) {
        index_[2]--;
      }
    }

  } else {
    calculateVelocity();
    x_ += dt * velocity_[0];
    y_ += dt * velocity_[1];

    // Updating x index
    if (velocity_[0] >= 0.0) {
      while (parameters_.meshsize->getPosX(index_[0], 0) + parameters_.meshsize->getDx(index_[0], 0) <= x_) {
        index_[0]++;
      }
    } else {
      while (parameters_.meshsize->getPosX(index_[0], 0) >= x_) {
        index_[0]--;
      }
    }

    // Updating y index
    if (velocity_[1] >= 0.0) {
      while (parameters_.meshsize->getPosY(0, index_[1]) + parameters_.meshsize->getDy(0, index_[1]) <= y_) {
        index_[1]++;
      }
    } else {
      while (parameters_.meshsize->getPosY(0, index_[1]) >= y_) {
        index_[1]--;
      }
    }
  }

  applyBoundaryCondition();
}

void Particle::applyBoundaryCondition() {
  /******* Cavity Case *******/
  if (parameters_.simulation.scenario == "cavity") {

    // Check if particle is left of left boundary
    if (x_ < 0) {
      wallCorrect(x_, 0.0, velocity_[0]);
    }
    // Check if particle is right of right boundary
    else if (x_ > parameters_.geometry.lengthX) {
      wallCorrect(x_, parameters_.geometry.lengthX, velocity_[0]);
    }

    // Check if particle is below bottom boundary
    if (y_ < 0) {
      wallCorrect(y_, 0.0, velocity_[1]);
    }
    // Check if particle is above top boundary
    else if (y_ > parameters_.geometry.lengthY) {
      wallCorrect(y_, parameters_.geometry.lengthY, velocity_[1]);
    }

    if (parameters_.geometry.dim == 3) {
      if (z_ < 0) {
        wallCorrect(z_, 0.0, velocity_[2]);
      } else if (z_ > parameters_.geometry.lengthZ) {
        wallCorrect(z_, parameters_.geometry.lengthZ, velocity_[2]);
      }
    }

  }
  /******* Channel Case ********/
  else if (parameters_.bfStep.yRatio < 0) {
    // Check if particle is below bottom boundary
    if (y_ < 0) {
      wallCorrect(y_, 0.0, velocity_[1]);
    }
    // Check if particle is above top boundary
    else if (y_ > parameters_.geometry.lengthY) {
      wallCorrect(y_, parameters_.geometry.lengthY, velocity_[1]);
    }

    if (parameters_.geometry.dim == 3) {
      if (z_ < 0) {
        wallCorrect(z_, 0.0, velocity_[2]);
      } else if (z_ > parameters_.geometry.lengthZ) {
        wallCorrect(z_, parameters_.geometry.lengthZ, velocity_[2]);
      }
    }
  }
  /******* BFS Case ********/
  else {
    RealType stepLengthX = parameters_.geometry.lengthX * parameters_.bfStep.xRatio;
    RealType stepLengthY = parameters_.geometry.lengthY * parameters_.bfStep.yRatio;
    // Check if particle is inside the step
    if (x_ < stepLengthX && y_ < stepLengthY) {
      if (velocity_[0] < 0)
        wallCorrect(x_, stepLengthX, velocity_[0]);
      if (velocity_[1] < 0)
        wallCorrect(y_, stepLengthY, velocity_[1]);
    }
    // Check if particle is above top boundary
    else if (y_ > parameters_.geometry.lengthY) {
      wallCorrect(y_, parameters_.geometry.lengthY, velocity_[1]);
    }
    // Check if particle is below bottom boundary
    else if (y_ < 0) {
      wallCorrect(y_, 0.0, velocity_[1]);
    }
    if (parameters_.geometry.dim == 3) {
      if (z_ < 0) {
        wallCorrect(z_, 0.0, velocity_[2]);
      } else if (z_ > parameters_.geometry.lengthZ) {
        wallCorrect(z_, parameters_.geometry.lengthZ, velocity_[2]);
      }
    }
  }
}

inline void Particle::wallCorrect(RealType& x, RealType limitX, RealType& velocity) {
  x        = 2 * limitX - x;
  velocity = -velocity;
}