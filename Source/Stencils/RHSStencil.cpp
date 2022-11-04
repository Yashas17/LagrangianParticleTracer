#include "StdAfx.hpp"

#include "RHSStencil.hpp"

Stencils::RHSStencil::RHSStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j) {
  const RealType dt = parameters_.timestep.dt;
  const RealType dx = 0.5 * (parameters_.meshsize->getDx(i, j) + parameters_.meshsize->getDx(i - 1, j));
  const RealType dy = 0.5 * (parameters_.meshsize->getDy(i, j) + parameters_.meshsize->getDy(i, j - 1));

  VectorField& FGH = flowField.getFGH();

  RealType dudx = (FGH.getVector(i, j)[0] - FGH.getVector(i - 1, j)[0]) / dx;
  RealType dvdy = (FGH.getVector(i, j)[1] - FGH.getVector(i, j - 1)[1]) / dy;

  flowField.getRHS().getScalar(i, j) = (dudx + dvdy) / dt;
}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j, int k) {
  const RealType dt = parameters_.timestep.dt;
  const RealType dx = 0.5 * (parameters_.meshsize->getDx(i, j, k) + parameters_.meshsize->getDx(i - 1, j, k));
  const RealType dy = 0.5 * (parameters_.meshsize->getDy(i, j, k) + parameters_.meshsize->getDy(i, j - 1, k));
  const RealType dz = 0.5 * (parameters_.meshsize->getDz(i, j, k) + parameters_.meshsize->getDz(i, j, k - 1));

  VectorField& FGH = flowField.getFGH();

  RealType dudx = (FGH.getVector(i, j, k)[0] - FGH.getVector(i - 1, j, k)[0]) / dx;
  RealType dvdy = (FGH.getVector(i, j, k)[1] - FGH.getVector(i, j - 1, k)[1]) / dy;
  RealType dwdz = (FGH.getVector(i, j, k)[2] - FGH.getVector(i, j, k - 1)[2]) / dz;

  flowField.getRHS().getScalar(i, j, k) = (dudx + dvdy + dwdz) / dt;
}