#include "StdAfx.hpp"

#include "RHSStencil.hpp"

Stencils::RHSStencil::RHSStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j) {

  VectorField& FGH = flowField.getFGH();

  RealType dfdx = (FGH.getVector(i, j)[0] - FGH.getVector(i - 1, j)[0]) / parameters_.meshsize->getDx(i, j);
  RealType dgdy = (FGH.getVector(i, j)[1] - FGH.getVector(i, j - 1)[1]) / parameters_.meshsize->getDy(i, j);

  flowField.getRHS().getScalar(i, j) = (dfdx + dgdy) / parameters_.timestep.dt;
}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j, int k) {

  VectorField& FGH = flowField.getFGH();

  RealType dfdx = (FGH.getVector(i, j, k)[0] - FGH.getVector(i - 1, j, k)[0]) / parameters_.meshsize->getDx(i, j, k);
  RealType dgdy = (FGH.getVector(i, j, k)[1] - FGH.getVector(i, j - 1, k)[1]) / parameters_.meshsize->getDy(i, j, k);
  RealType dhdz = (FGH.getVector(i, j, k)[2] - FGH.getVector(i, j, k - 1)[2]) / parameters_.meshsize->getDz(i, j, k);

  flowField.getRHS().getScalar(i, j, k) = (dfdx + dgdy + dhdz) / parameters_.timestep.dt;
}