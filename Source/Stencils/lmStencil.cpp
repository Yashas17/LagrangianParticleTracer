#include "lmStencil.hpp"

#include <algorithm>
#include <cmath>

Stencils::lmStencil::lmStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::lmStencil::apply(FlowField& flowField, int i, int j) {

  RealType* const h  = flowField.getH().getScalar(i, j);
  RealType* const lm = flowField.getLm().getScalar(i, j);

  RealType delta = -1;

  if (parameters_.turbulence.boundaryLayer == "turbulent") {
    const RealType posX = parameters_.meshsize->getPosX(i, j);
    RealType       Rex  = parameters_.flow.Re * posX; // Flat Plate Reynolds Number
    delta               = 0.382 * posX / std::pow(Rex, 0.2);
  }
  elseif(parameters_.turbulence.boundaryLayer == "laminar") {
    const RealType posX = parameters_.meshsize->getPosX(i, j);
    RealType       Rex  = parameters_.flow.Re * posX; // Flat Plate Reynolds Number
    delta               = 4.91 * posX / std::pow(Rex, 0.5);
  }
  elseif(parameters_.turbulence.boundaryLayer == "none") { delta = 0; }
  else {
    throw std::runtime_error("Undefined Boundary Layer Function");
  }

  lm = std::min(k_ * h, 0.09 * delta);
}

void Stencils::lmStencil::apply(FlowField& flowField, int i, int j, int k) {

  RealType* const h  = flowField.getH().getScalar(i, j, k);
  RealType* const lm = flowField.getLm().getScalar(i, j, k);

  RealType delta = -1;

  if (parameters_.turbulence.boundaryLayer == "turbulent") {
    const RealType posX = parameters_.meshsize->getPosX(i, j, k);
    RealType       Rex  = parameters_.flow.Re * posX; // Flat Plate Reynolds Number
    delta               = 0.382 * posX / std::pow(Rex, 0.2);
  }
  elseif(parameters_.turbulence.boundaryLayer == "laminar") {
    const RealType posX = parameters_.meshsize->getPosX(i, j, k);
    RealType       Rex  = parameters_.flow.Re * posX; // Flat Plate Reynolds Number
    delta               = 4.91 * posX / std::pow(Rex, 0.5);
  }
  elseif(parameters_.turbulence.boundaryLayer == "none") { delta = 0; }
  else {
    throw std::runtime_error("Undefined Boundary Layer Function");
  }

  lm = std::min(k_ * h, 0.09 * delta);
}
