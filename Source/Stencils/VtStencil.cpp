#include "StdAfx.hpp"

#include "VtStencil.hpp"

Stencils::vtStencil::vtStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::vtStencil::apply(FlowField& flowField, int i, int j) {
  RealType* const vt = flowField.getVt().getScalar(i, j);
  RealType* const lm = flowField.getLm().getScalar(i, j);
  loadLocalVelocity2D(flowField, localVelocity_, i, j);
  loadLocalMeshsize2D(parameters_, localMeshsize_, i, j);

  RealType S11 = dudx(localVelocity_, localMeshsize_);
  RealType S22 = dvdy(localVelocity_, localMeshsize_);
  RealType S12 = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));

  vt = lm * lm * sqrt(S11 * S11 + S22 * S22 + 2 * S12 * S12);
}

void Stencils::vtStencil::apply(FlowField& flowField, int i, int j, int k) {
  RealType* const vt = flowField.getVt().getScalar(i, j, k);
  RealType* const lm = flowField.getLm().getScalar(i, j, k);
  loadLocalVelocity2D(flowField, localVelocity_, i, j, k);
  loadLocalMeshsize2D(parameters_, localMeshsize_, i, j, k);

  RealType S11 = dudx(localVelocity_, localMeshsize_);
  RealType S22 = dvdy(localVelocity_, localMeshsize_);
  RealType S33 = dwdz(localVelocity_, localMeshsize_);
  RealType S12 = 0.5 * (dudy(localVelocity_, localMeshsize_) + dvdx(localVelocity_, localMeshsize_));
  RealType S13 = 0.5 * (dudz(localVelocity_, localMeshsize_) + dwdx(localVelocity_, localMeshsize_));
  RealType S23 = 0.5 * (dwdy(localVelocity_, localMeshsize_) + dvdz(localVelocity_, localMeshsize_));

  vt = lm * lm * sqrt(S11 * S11 + S22 * S22 + 2 * (S12 * S12 + S13 * S13 + S23 * S23));
}
