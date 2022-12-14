#include "StdAfx.hpp"

#include "VtBufferReadStencil.hpp"

#include <vector>

// 2d constructor
Stencils::VtBufferReadStencil::VtBufferReadStencil(
  const Parameters&       parameters,
  std::vector<RealType>&& leftBufferFrom,
  std::vector<RealType>&& rightBufferFrom,
  std::vector<RealType>&& bottomBufferFrom,
  std::vector<RealType>&& topBufferFrom
):
  BoundaryStencil<FlowField>(parameters) {
  // size includes the ghost cells
  cellsX = parameters.parallel.localSize[0] + 3;
  cellsY = parameters.parallel.localSize[1] + 3;

  leftBuffer   = std::move(leftBufferFrom);
  rightBuffer  = std::move(rightBufferFrom);
  bottomBuffer = std::move(bottomBufferFrom);
  topBuffer    = std::move(topBufferFrom);
}

// 3d constructor
Stencils::VtBufferReadStencil::VtBufferReadStencil(
  const Parameters&       parameters,
  std::vector<RealType>&& leftBufferFrom,
  std::vector<RealType>&& rightBufferFrom,
  std::vector<RealType>&& bottomBufferFrom,
  std::vector<RealType>&& topBufferFrom,
  std::vector<RealType>&& frontBufferFrom,
  std::vector<RealType>&& backBufferFrom
):
  BoundaryStencil<FlowField>(parameters) {
  // size includes the ghost cells
  cellsX = parameters.parallel.localSize[0] + 3;
  cellsY = parameters.parallel.localSize[1] + 3;
  cellsZ = parameters.parallel.localSize[2] + 3;

  leftBuffer   = std::move(leftBufferFrom);
  rightBuffer  = std::move(rightBufferFrom);
  bottomBuffer = std::move(bottomBufferFrom);
  topBuffer    = std::move(topBufferFrom);
  frontBuffer  = std::move(frontBufferFrom);
  backBuffer   = std::move(backBufferFrom);
}

// 2d cases
void Stencils::VtBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVt().getNx());
  assert(cellsY == flowField.getVt().getNy());
  flowField.getVt().getScalar(i, j) = leftBuffer[j]; // j changes
}

void Stencils::VtBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVt().getNx());
  assert(cellsY == flowField.getVt().getNy());
  flowField.getVt().getScalar(i, j) = rightBuffer[j]; // j changes
}

void Stencils::VtBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVt().getNx());
  assert(cellsY == flowField.getVt().getNy());
  flowField.getVt().getScalar(i, j) = bottomBuffer[i]; // i changes
}

void Stencils::VtBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVt().getNx());
  assert(cellsY == flowField.getVt().getNy());
  flowField.getVt().getScalar(i, j) = topBuffer[i]; // i changes
}

// 3d cases
void Stencils::VtBufferReadStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVt().getScalar(i, j, k) = leftBuffer[j * cellsZ + k]; // j,k changes
}

void Stencils::VtBufferReadStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVt().getScalar(i, j, k) = rightBuffer[j * cellsZ + k]; // j,k changes
}

void Stencils::VtBufferReadStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVt().getScalar(i, j, k) = bottomBuffer[i * cellsZ + k]; // i,k changes
}

void Stencils::VtBufferReadStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVt().getScalar(i, j, k) = topBuffer[i * cellsZ + k]; // i,k changes
}

void Stencils::VtBufferReadStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVt().getScalar(i, j, k) = frontBuffer[i * cellsY + j]; // i,j changes
}

void Stencils::VtBufferReadStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVt().getScalar(i, j, k) = backBuffer[i * cellsY + j]; // i,j changes
}