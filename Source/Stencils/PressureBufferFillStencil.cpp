#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include <vector>

Stencils::PressureBufferFillStencil::PressureBufferFillStencil(const Parameters& parameters):
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
void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getPressure().getNx());
  assert(cellsY == flowField.getPressure().getNy());
  leftBuffer[j] = flowField.getPressure().getScalar(i + 1, j); // j changes
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getPressure().getNx());
  assert(cellsY == flowField.getPressure().getNy());
  rightBuffer[j] = flowField.getPressure().getScalar(i - 1, j); // j changes
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getPressure().getNx());
  assert(cellsY == flowField.getPressure().getNy());
  bottomBuffer[i] = flowField.getPressure().getScalar(i, j + 1); // i changes
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getPressure().getNx());
  assert(cellsY == flowField.getPressure().getNy());
  topBuffer[i] = flowField.getPressure().getScalar(i, j - 1); // i changes
}

// 3d cases
void Stencils::PressureBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  leftBuffer[j * cellsZ + k] = flowField.getPressure().getScalar(i + 1, j, k); // j,k changes
}

void Stencils::PressureBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  rightBuffer[j * cellsZ + k] = flowField.getPressure().getScalar(i - 1, j, k); // j,k changes
}

void Stencils::PressureBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  bottomBuffer[i * cellsZ + k] = flowField.getPressure().getScalar(i, j + 1, k); // i,k changes
}

void Stencils::PressureBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  topBuffer[i * cellsZ + k] = flowField.getPressure().getScalar(i, j - 1, k); // i,k changes
}

void Stencils::PressureBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  frontBuffer[i * cellsY + j] = flowField.getPressure().getScalar(i, j, k + 1); // i,j changes
}

void Stencils::PressureBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  backBuffer[i * cellsY + j] = flowField.getPressure().getScalar(i, j, k - 1); // i,j changes
}