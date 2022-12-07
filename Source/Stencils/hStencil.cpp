#include "StdAfx.hpp"

#include "hStencil.hpp"

#include <algorithm>

Stencils::hStencil::hStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  xLimit_(parameters.bfStep.xRatio * parameters.geometry.lengthX),
  yLimit_(parameters.bfStep.yRatio * parameters.geometry.lengthY) {}

void Stencils::hStencil::apply(FlowField& flowField, int i, int j) {
  RealType* const value = flowField.getH().getScalar(i, j);

  const RealType  posX   = parameters_.meshsize->getPosX(i, j);
  const RealType  posY   = parameters_.meshsize->getPosY(i, j);

  value = std::min(posX, posY);
}

void Stencils::hStencil::apply(FlowField& flowField, int i, int j, int k) {
  RealType* const value = flowField.getH().getScalar(i, j, k);

  const RealType  posX   = parameters_.meshsize->getPosX(i, j, k);
  const RealType  posY   = parameters_.meshsize->getPosY(i, j, k);
  const RealType  posZ   = parameters_.meshsize->getPosZ(i, j, k);

  value = std::min(posX, std::min(posY, posZ));
}
