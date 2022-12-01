#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"
#include <vector>

Stencils::VelocityBufferFillStencil::VelocityBufferFillStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {

    //allocate buffer arrays
    if(parameters.geometry.dim == 2){

      //size includes the ghost cells
      cellsX = parameters.parallel.localSize[0];
      cellsY = parameters.parallel.localSize[1];

      //2 times for both velocity directions, flattened
      std::vector<RealType> leftBuffer(2 * cellsY);
      std::vector<RealType> rightBuffer(2 * cellsY);
      std::vector<RealType> bottomBuffer(2 * cellsX);
      std::vector<RealType> topBuffer(2 * cellsX);

    } else if (parameters.geometry.dim == 3){

      //size includes the ghost cells
      cellsX = parameters.parallel.localSize[0];
      cellsY = parameters.parallel.localSize[1];
      cellsZ = parameters.parallel.localSize[2];

      //3 times for all velocity directions, flattened
      std::vector<RealType> leftBuffer(3 * cellsY*cellsZ);
      std::vector<RealType> rightBuffer(3 * cellsY*cellsZ);
      std::vector<RealType> bottomBuffer(3 * cellsX*cellsZ);
      std::vector<RealType> topBuffer(3 * cellsX*cellsZ);
      std::vector<RealType> frontBuffer(3 * cellsX*cellsY);
      std::vector<RealType> backBuffer(3 * cellsX*cellsY);

    }
    
  }

//2d cases
void Stencils::VelocityBufferFillStencil::applyLeftWall(
  FlowField& flowField, int i, int j
) {
  leftBuffer[2 * j] = flowField.getVelocity().getVector(i,j)[0]; //j changes
  leftBuffer[2 * j + 1] = flowField.getVelocity().getVector(i,j)[1]; //j changes
}

void Stencils::VelocityBufferFillStencil::applyRightWall(
  FlowField& flowField, int i, int j
) {
  rightBuffer[2 * j] = flowField.getVelocity().getVector(i,j)[0]; //j changes
  rightBuffer[2 * j + 1] = flowField.getVelocity().getVector(i,j)[1]; //j changes
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(
  FlowField& flowField, int i, int j
) {
  bottomBuffer[2 * i] = flowField.getVelocity().getVector(i,j)[0]; //i changes
  bottomBuffer[2 * i + 1] = flowField.getVelocity().getVector(i,j)[1]; //i changes
}

void Stencils::VelocityBufferFillStencil::applyTopWall(
  FlowField& flowField, int i, int j
) {
  topBuffer[2 * i] = flowField.getVelocity().getVector(i,j)[0]; //i changes
  topBuffer[2 * i + 1] = flowField.getVelocity().getVector(i,j)[1]; //i changes
}

//3d cases
void Stencils::VelocityBufferFillStencil::applyLeftWall(
  FlowField& flowField, int i, int j, int k
) {
  leftBuffer[3 * (j*cellsZ + k)] = flowField.getVelocity().getVector(i,j,k)[0]; //j,k changes
  leftBuffer[3 * (j*cellsZ + k) + 1] = flowField.getVelocity().getVector(i,j,k)[1]; //j,k changes
  leftBuffer[3 * (j*cellsZ + k) + 2] = flowField.getVelocity().getVector(i,j,k)[2]; //j,k changes
}

void Stencils::VelocityBufferFillStencil::applyRightWall(
  FlowField& flowField, int i, int j, int k
) {
  rightBuffer[3 * (j*cellsZ + k)] = flowField.getVelocity().getVector(i,j,k)[0]; //j,k changes
  rightBuffer[3 * (j*cellsZ + k) + 1] = flowField.getVelocity().getVector(i,j,k)[1]; //j,k changes
  rightBuffer[3 * (j*cellsZ + k) + 2] = flowField.getVelocity().getVector(i,j,k)[2]; //j,k changes
}

void Stencils::VelocityBufferFillStencil::applyBottomWall(
  FlowField& flowField, int i, int j, int k
) {
  bottomBuffer[3 * (i*cellsZ + k)] = flowField.getVelocity().getVector(i,j,k)[0]; //i,k changes
  bottomBuffer[3 * (i*cellsZ + k) + 1] = flowField.getVelocity().getVector(i,j,k)[1]; //i,k changes
  bottomBuffer[3 * (i*cellsZ + k) + 2] = flowField.getVelocity().getVector(i,j,k)[2]; //i,k changes
}

void Stencils::VelocityBufferFillStencil::applyTopWall(
  FlowField& flowField, int i, int j, int k
) {
  topBuffer[3 * (i*cellsZ + k)] = flowField.getVelocity().getVector(i,j,k)[0]; //i,k changes
  topBuffer[3 * (i*cellsZ + k) + 1] = flowField.getVelocity().getVector(i,j,k)[1]; //i,k changes
  topBuffer[3 * (i*cellsZ + k) + 2] = flowField.getVelocity().getVector(i,j,k)[2]; //i,k changes
}

void Stencils::VelocityBufferFillStencil::applyFrontWall(
  FlowField& flowField, int i, int j, int k
) {
  frontBuffer[3 * (i*cellsY + j)] = flowField.getVelocity().getVector(i,j,k)[0]; //i,j changes
  frontBuffer[3 * (i*cellsY + j) + 1] = flowField.getVelocity().getVector(i,j,k)[1]; //i,j changes
  frontBuffer[3 * (i*cellsY + j) + 2] = flowField.getVelocity().getVector(i,j,k)[2]; //i,j changes
}

void Stencils::VelocityBufferFillStencil::applyBackWall(
  FlowField& flowField, int i, int j, int k
) {
  backBuffer[3 * (i*cellsY + j)] = flowField.getVelocity().getVector(i,j,k)[0]; //i,j changes
  backBuffer[3 * (i*cellsY + j) + 1] = flowField.getVelocity().getVector(i,j,k)[1]; //i,j changes
  backBuffer[3 * (i*cellsY + j) + 2] = flowField.getVelocity().getVector(i,j,k)[2]; //i,j changes
}