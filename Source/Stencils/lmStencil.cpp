#include "lmStencil.hpp"

#include <algorithm>
#include <cmath>

Stencils::lmStencil::lmStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::lmStencil::apply(FlowField& flowField, int i, int j) {

  RealType* const h  = flowField.getH().getScalar(i, j);
  RealType* const lm = flowField.getLm().getScalar(i, j);

  const RealType posX  = parameters_.meshsize->getPosX(i, j);
  RealType       Rex   = parameters_.flow.Re * posX;        // Flat Plate Reynolds Number
  RealType       delta = 0.382 * posX / std::pow(Rex, 0.2); // Boundary Layer Thickness

  // RealType       delta = 4.91 * posX / std::pow(Rex, 0.5); // Boundary Layer Thickness
  // RealType       delta = 0; // Boundary Layer Thickness

  lm = std::min(k_ * h, 0.09 * delta);
}

void Stencils::lmStencil::apply(FlowField& flowField, int i, int j, int k) {

  RealType* const h  = flowField.getH().getScalar(i, j, k);
  RealType* const lm = flowField.getLm().getScalar(i, j, k);

  const RealType posX  = parameters_.meshsize->getPosX(i, j, k);
  RealType       Rex   = parameters_.flow.Re * posX;        // Flat Plate Reynolds Number
  RealType       delta = 0.382 * posX / std::pow(Rex, 0.2); // Boundary Layer Thickness

  // RealType       delta = 4.91 * posX / std::pow(Rex, 0.5); // Boundary Layer Thickness
  // RealType       delta = 0; // Boundary Layer Thickness

  lm = std::min(k_ * h, 0.09 * delta);
}
