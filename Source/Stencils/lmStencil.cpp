#include "StdAfx.hpp"

#include "lmStencil.hpp"

Stencils::lmStencil::lmStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::lmStencil::apply(FlowField& flowField, int i, int j) {

  const RealType h  = flowField.getH().getScalar(i, j);
  RealType&       lm = flowField.getLm().getScalar(i, j);

  RealType delta = -1;

  if (parameters_.turbulence.boundaryLayer == 2) {
    const RealType posX = parameters_.meshsize->getPosX(i, j) + parameters_.meshsize->getDx(i,j)/2;
    RealType       Rex  = parameters_.flow.Re * posX; // Flat Plate Reynolds Number
    
    delta = (posX <= 0)  ? 0 : 0.382 * posX / std::pow(Rex, 0.2);
  
  } else if (parameters_.turbulence.boundaryLayer == 1) {
    const RealType posX = parameters_.meshsize->getPosX(i, j) + parameters_.meshsize->getDx(i,j)/2;
    RealType       Rex  = parameters_.flow.Re * posX; // Flat Plate Reynolds Number
    
    delta = (posX <= 0)  ? 0 : 4.91 * posX / std::pow(Rex, 0.5);

  } else if (parameters_.turbulence.boundaryLayer == 0) {
    delta = 0;
  } else {
    throw std::runtime_error("Undefined Boundary Layer Function");
  }

  lm = std::min(k_ * h, 0.09 * delta);
}

void Stencils::lmStencil::apply(FlowField& flowField, int i, int j, int k) {

  const RealType h  = flowField.getH().getScalar(i, j, k);
  RealType&       lm = flowField.getLm().getScalar(i, j, k);

  RealType delta = -1;

  if (parameters_.turbulence.boundaryLayer == 2) {
    const RealType posX = parameters_.meshsize->getPosX(i, j, k) + parameters_.meshsize->getDx(i,j,k)/2;
    RealType       Rex  = parameters_.flow.Re * posX; // Flat Plate Reynolds Number

    delta = (posX <= 0)  ? 0 : 0.382 * posX / std::pow(Rex, 0.2);

  } else if (parameters_.turbulence.boundaryLayer == 1) {
    const RealType posX = parameters_.meshsize->getPosX(i, j, k) + parameters_.meshsize->getDx(i,j,k)/2;
    RealType       Rex  = parameters_.flow.Re * posX; // Flat Plate Reynolds Number

    delta = (posX <= 0)  ? 0 : 4.91 * posX / std::pow(Rex, 0.5);

  } else if (parameters_.turbulence.boundaryLayer == 0) {
    delta = 0;
  } else {
    throw std::runtime_error("Undefined Boundary Layer Function");
  }

  lm = std::min(k_ * h, 0.09 * delta);
}
