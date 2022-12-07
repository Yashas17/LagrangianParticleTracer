#include "StdAfx.hpp"

#include "hStencil.hpp"

Stencils::hStencil::hStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  xLimit_(parameters.bfStep.xRatio * parameters.geometry.lengthX),
  yLimit_(parameters.bfStep.yRatio * parameters.geometry.lengthY) {}

  // TODO: Implement apply function 