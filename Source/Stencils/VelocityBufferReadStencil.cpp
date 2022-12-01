#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"
#include <vector>

//2d constructor
Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil
  (const Parameters& parameters,std::vector<RealType>& leftBufferFrom, 
                                std::vector<RealType>& rightBufferFrom,
                                std::vector<RealType>& bottomBufferFrom,
                                std::vector<RealType>& topBufferFrom):
  BoundaryStencil<FlowField>(parameters) {
    //size includes the ghost cells
    cellsX = parameters.parallel.localSize[0];
    cellsY = parameters.parallel.localSize[1];

    leftBuffer = std::move(leftBufferFrom);
    rightBuffer = std::move(rightBufferFrom);
    bottomBuffer = std::move(bottomBufferFrom);
    topBuffer = std::move(topBufferFrom);
  }

//3d constructor
Stencils::VelocityBufferReadStencil::VelocityBufferReadStencil
  (const Parameters& parameters,std::vector<RealType>& leftBufferFrom, 
                                std::vector<RealType>& rightBufferFrom,
                                std::vector<RealType>& bottomBufferFrom,
                                std::vector<RealType>& topBufferFrom,
                                std::vector<RealType>& frontBufferFrom,
                                std::vector<RealType>& backBufferFrom):
  BoundaryStencil<FlowField>(parameters) {
    //size includes the ghost cells
    cellsX = parameters.parallel.localSize[0];
    cellsY = parameters.parallel.localSize[1];
    cellsZ = parameters.parallel.localSize[2];

    leftBuffer = std::move(leftBufferFrom);
    rightBuffer = std::move(rightBufferFrom);
    bottomBuffer = std::move(bottomBufferFrom);
    topBuffer = std::move(topBufferFrom);
    frontBuffer = std::move(frontBufferFrom);
    backBuffer = std::move(backBufferFrom);
  }

//2d cases
void Stencils::VelocityBufferReadStencil::applyLeftWall(
  FlowField& flowField, int i, int j
) {
  flowField.getVelocity().getVector(i,j)[0] = leftBuffer[2 * j]; //j changes
  flowField.getVelocity().getVector(i,j)[1] = leftBuffer[2 * j + 1]; //j changes
}

void Stencils::VelocityBufferReadStencil::applyRightWall(
  FlowField& flowField, int i, int j
) {
  flowField.getVelocity().getVector(i,j)[0] = rightBuffer[2 * j]; //j changes
  flowField.getVelocity().getVector(i,j)[1] = rightBuffer[2 * j + 1]; //j changes
}

void Stencils::VelocityBufferReadStencil::applyBottomWall(
  FlowField& flowField, int i, int j
) {
  flowField.getVelocity().getVector(i,j)[0] = bottomBuffer[2 * i]; //i changes
  flowField.getVelocity().getVector(i,j)[1] = bottomBuffer[2 * i + 1]; //i changes
}

void Stencils::VelocityBufferReadStencil::applyTopWall(
  FlowField& flowField, int i, int j
) {
  flowField.getVelocity().getVector(i,j)[0] = topBuffer[2 * i]; //i changes
  flowField.getVelocity().getVector(i,j)[1] = topBuffer[2 * i + 1]; //i changes
}

//3d cases
void Stencils::VelocityBufferReadStencil::applyLeftWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getVelocity().getVector(i,j,k)[0] = leftBuffer[3 * (j*cellsZ + k)]; //j,k changes
  flowField.getVelocity().getVector(i,j,k)[1] = leftBuffer[3 * (j*cellsZ + k) + 1]; //j,k changes
  flowField.getVelocity().getVector(i,j,k)[2] = leftBuffer[3 * (j*cellsZ + k) + 2]; //j,k changes
}

void Stencils::VelocityBufferReadStencil::applyRightWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getVelocity().getVector(i,j,k)[0] = rightBuffer[3 * (j*cellsZ + k)]; //j,k changes
  flowField.getVelocity().getVector(i,j,k)[1] = rightBuffer[3 * (j*cellsZ + k) + 1]; //j,k changes
  flowField.getVelocity().getVector(i,j,k)[2] = rightBuffer[3 * (j*cellsZ + k) + 2]; //j,k changes
}

void Stencils::VelocityBufferReadStencil::applyBottomWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getVelocity().getVector(i,j,k)[0] = bottomBuffer[3 * (i*cellsZ + k)]; //i,k changes
  flowField.getVelocity().getVector(i,j,k)[1] = bottomBuffer[3 * (i*cellsZ + k) + 1]; //i,k changes
  flowField.getVelocity().getVector(i,j,k)[2] = bottomBuffer[3 * (i*cellsZ + k) + 2]; //i,k changes
}

void Stencils::VelocityBufferReadStencil::applyTopWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getVelocity().getVector(i,j,k)[0] = topBuffer[3 * (i*cellsZ + k)]; //i,k changes
  flowField.getVelocity().getVector(i,j,k)[1] = topBuffer[3 * (i*cellsZ + k) + 1]; //i,k changes
  flowField.getVelocity().getVector(i,j,k)[2] = topBuffer[3 * (i*cellsZ + k) + 2]; //i,k changes
}

void Stencils::VelocityBufferReadStencil::applyFrontWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getVelocity().getVector(i,j,k)[0] = frontBuffer[3 * (i*cellsY + j)]; //i,j changes
  flowField.getVelocity().getVector(i,j,k)[1] = frontBuffer[3 * (i*cellsY + j) + 1]; //i,j changes
  flowField.getVelocity().getVector(i,j,k)[2] = frontBuffer[3 * (i*cellsY + j) + 2]; //i,j changes
}

void Stencils::VelocityBufferReadStencil::applyBackWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getVelocity().getVector(i,j,k)[0] = backBuffer[3 * (i*cellsY + j)]; //i,j changes
  flowField.getVelocity().getVector(i,j,k)[1] = backBuffer[3 * (i*cellsY + j) + 1]; //i,j changes
  flowField.getVelocity().getVector(i,j,k)[2] = backBuffer[3 * (i*cellsY + j) + 2]; //i,j changes
}