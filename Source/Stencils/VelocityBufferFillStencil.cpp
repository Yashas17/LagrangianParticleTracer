#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

#include <vector>

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {
  // allocate buffer arrays
  if (parameters.geometry.dim == 2) {

    // size includes the ghost cells
    cellsX = parameters.parallel.localSize[0] + 3;
    cellsY = parameters.parallel.localSize[1] + 3;

    // 2 times for both velocity directions, flattened
    leftBuffer.resize(2 * cellsY);   // u,v
    rightBuffer.resize(3 * cellsY);  // u,v, another u
    bottomBuffer.resize(2 * cellsX); // u,v
    topBuffer.resize(3 * cellsX);    // u,v, another v

  } else if (parameters.geometry.dim == 3) {

    // size includes the ghost cells
    cellsX = parameters.parallel.localSize[0] + 3;
    cellsY = parameters.parallel.localSize[1] + 3;
    cellsZ = parameters.parallel.localSize[2] + 3;

    // 3 times for all velocity directions, flattened
    leftBuffer.resize(3 * cellsY * cellsZ);   // u,v,w
    rightBuffer.resize(4 * cellsY * cellsZ);  // u,v,w, another u
    bottomBuffer.resize(3 * cellsX * cellsZ); // u,v,w
    topBuffer.resize(4 * cellsX * cellsZ);    // u,v,w, another v
    frontBuffer.resize(3 * cellsX * cellsY);  // u,v,w
    backBuffer.resize(4 * cellsX * cellsY);   // u,v,w another w
  }
}

// 2d cases
void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVelocity().getNx());
  assert(cellsY == flowField.getVelocity().getNy());
  leftBuffer[2 * j]     = flowField.getVelocity().getVector(i + 1, j)[0]; // j changes
  leftBuffer[2 * j + 1] = flowField.getVelocity().getVector(i + 1, j)[1]; // j changes
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVelocity().getNx());
  assert(cellsY == flowField.getVelocity().getNy());
  rightBuffer[3 * j]     = flowField.getVelocity().getVector(i, j)[0];     // j changes
  rightBuffer[3 * j + 1] = flowField.getVelocity().getVector(i, j)[1];     // j changes
  rightBuffer[3 * j + 2] = flowField.getVelocity().getVector(i - 1, j)[0]; // another u
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVelocity().getNx());
  assert(cellsY == flowField.getVelocity().getNy());
  bottomBuffer[2 * i]     = flowField.getVelocity().getVector(i, j + 1)[0]; // i changes
  bottomBuffer[2 * i + 1] = flowField.getVelocity().getVector(i, j + 1)[1]; // i changes
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j) {
  assert(cellsX == flowField.getVelocity().getNx());
  assert(cellsY == flowField.getVelocity().getNy());
  topBuffer[3 * i]     = flowField.getVelocity().getVector(i, j)[0];     // i changes
  topBuffer[3 * i + 1] = flowField.getVelocity().getVector(i, j)[1];     // i changes
  topBuffer[3 * i + 2] = flowField.getVelocity().getVector(i, j - 1)[1]; // another v
}

// 3d cases
void Stencils::VelocityBufferFillStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  leftBuffer[3 * (j * cellsZ + k)]     = flowField.getVelocity().getVector(i + 1, j, k)[0]; // j,k changes
  leftBuffer[3 * (j * cellsZ + k) + 1] = flowField.getVelocity().getVector(i + 1, j, k)[1]; // j,k changes
  leftBuffer[3 * (j * cellsZ + k) + 2] = flowField.getVelocity().getVector(i + 1, j, k)[2]; // j,k changes
}

void Stencils::VelocityBufferFillStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  rightBuffer[4 * (j * cellsZ + k)]     = flowField.getVelocity().getVector(i, j, k)[0];     // j,k changes
  rightBuffer[4 * (j * cellsZ + k) + 1] = flowField.getVelocity().getVector(i, j, k)[1];     // j,k changes
  rightBuffer[4 * (j * cellsZ + k) + 2] = flowField.getVelocity().getVector(i, j, k)[2];     // j,k changes
  rightBuffer[4 * (j * cellsZ + k) + 3] = flowField.getVelocity().getVector(i - 1, j, k)[0]; // another u
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  bottomBuffer[3 * (i * cellsZ + k)]     = flowField.getVelocity().getVector(i, j + 1, k)[0]; // i,k changes
  bottomBuffer[3 * (i * cellsZ + k) + 1] = flowField.getVelocity().getVector(i, j + 1, k)[1]; // i,k changes
  bottomBuffer[3 * (i * cellsZ + k) + 2] = flowField.getVelocity().getVector(i, j + 1, k)[2]; // i,k changes
}

void Stencils::VelocityBufferFillStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  topBuffer[4 * (i * cellsZ + k)]     = flowField.getVelocity().getVector(i, j, k)[0];     // i,k changes
  topBuffer[4 * (i * cellsZ + k) + 1] = flowField.getVelocity().getVector(i, j, k)[1];     // i,k changes
  topBuffer[4 * (i * cellsZ + k) + 2] = flowField.getVelocity().getVector(i, j, k)[2];     // i,k changes
  topBuffer[4 * (i * cellsZ + k) + 3] = flowField.getVelocity().getVector(i, j - 1, k)[1]; // another v
}

void Stencils::VelocityBufferFillStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  frontBuffer[3 * (i * cellsY + j)]     = flowField.getVelocity().getVector(i, j, k + 1)[0]; // i,j changes
  frontBuffer[3 * (i * cellsY + j) + 1] = flowField.getVelocity().getVector(i, j, k + 1)[1]; // i,j changes
  frontBuffer[3 * (i * cellsY + j) + 2] = flowField.getVelocity().getVector(i, j, k + 1)[2]; // i,j changes
}

void Stencils::VelocityBufferFillStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  backBuffer[4 * (i * cellsY + j)]     = flowField.getVelocity().getVector(i, j, k)[0];     // i,j changes
  backBuffer[4 * (i * cellsY + j) + 1] = flowField.getVelocity().getVector(i, j, k)[1];     // i,j changes
  backBuffer[4 * (i * cellsY + j) + 2] = flowField.getVelocity().getVector(i, j, k)[2];     // i,j changes
  backBuffer[4 * (i * cellsY + j) + 3] = flowField.getVelocity().getVector(i, j, k - 1)[2]; // another w
}