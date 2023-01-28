#include "StdAfx.hpp"

#include "VtBufferFillStencil.hpp"

#include <vector>

Stencils::VtBufferFillStencil::VtBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {

  // allocate buffer arrays
  if (parameters.geometry.dim == 2) {

    // size includes the ghost cells
    cellsX = parameters.parallel.localSize[0] + 3;
    cellsY = parameters.parallel.localSize[1] + 3;

    leftBuffer.resize(cellsY);
    rightBuffer.resize(cellsY);
    bottomBuffer.resize(cellsX);
    topBuffer.resize(cellsX);

  } else if (parameters.geometry.dim == 3) {

    // size includes the ghost cells
    cellsX = parameters.parallel.localSize[0] + 3;
    cellsY = parameters.parallel.localSize[1] + 3;
    cellsZ = parameters.parallel.localSize[2] + 3;

    leftBuffer.resize(cellsY * cellsZ);
    rightBuffer.resize(cellsY * cellsZ);
    bottomBuffer.resize(cellsX * cellsZ);
    topBuffer.resize(cellsX * cellsZ);
    frontBuffer.resize(cellsX * cellsY);
    backBuffer.resize(cellsX * cellsY);
  }
}

// 2d cases
void Stencils::VtBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVt().getNx());
  assert(cellsY == flowField.getVt().getNy());
  leftBuffer[j] = flowField.getVt().getScalar(i + 1, j); // j changes
}

void Stencils::VtBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVt().getNx());
  assert(cellsY == flowField.getVt().getNy());
  rightBuffer[j] = flowField.getVt().getScalar(i - 1, j); // j changes
}

void Stencils::VtBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVt().getNx());
  assert(cellsY == flowField.getVt().getNy());
  bottomBuffer[i] = flowField.getVt().getScalar(i, j + 1); // i changes
}

void Stencils::VtBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVt().getNx());
  assert(cellsY == flowField.getVt().getNy());
  topBuffer[i] = flowField.getVt().getScalar(i, j - 1); // i changes
}

// 3d cases
void Stencils::VtBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  leftBuffer[j * cellsZ + k] = flowField.getVt().getScalar(i + 1, j, k); // j,k changes
}

void Stencils::VtBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  rightBuffer[j * cellsZ + k] = flowField.getVt().getScalar(i - 1, j, k); // j,k changes
}

void Stencils::VtBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  bottomBuffer[i * cellsZ + k] = flowField.getVt().getScalar(i, j + 1, k); // i,k changes
}

void Stencils::VtBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  topBuffer[i * cellsZ + k] = flowField.getVt().getScalar(i, j - 1, k); // i,k changes
}

void Stencils::VtBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  frontBuffer[i * cellsY + j] = flowField.getVt().getScalar(i, j, k + 1); // i,j changes
}

void Stencils::VtBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  backBuffer[i * cellsY + j] = flowField.getVt().getScalar(i, j, k - 1); // i,j changes
}