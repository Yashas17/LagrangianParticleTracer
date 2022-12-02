#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"
#include <vector>

//2d constructor
Stencils::PressureBufferReadStencil::PressureBufferReadStencil
  (const Parameters& parameters,std::vector<RealType>& leftBufferFrom, 
                                std::vector<RealType>& rightBufferFrom,
                                std::vector<RealType>& bottomBufferFrom,
                                std::vector<RealType>& topBufferFrom):
  BoundaryStencil<FlowField>(parameters) {
    //size includes the ghost cells
    cellsX = parameters.parallel.localSize[0];
    cellsY = parameters.parallel.localSize[1];

    leftBuffer = leftBufferFrom;
    rightBuffer = rightBufferFrom;
    bottomBuffer = bottomBufferFrom;
    topBuffer = topBufferFrom;
  }

//3d constructor
Stencils::PressureBufferReadStencil::PressureBufferReadStencil
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

    leftBuffer = leftBufferFrom;
    rightBuffer = rightBufferFrom;
    bottomBuffer = bottomBufferFrom;
    topBuffer = topBufferFrom;
    frontBuffer = frontBufferFrom;
    backBuffer = backBufferFrom;
  }

//2d cases
void Stencils::PressureBufferReadStencil::applyLeftWall(
  FlowField& flowField, int i, int j
) {
  flowField.getPressure().getScalar(i - 1,j) = leftBuffer[j]; //j changes, i-1 for ghost layer
}

void Stencils::PressureBufferReadStencil::applyRightWall(
  FlowField& flowField, int i, int j
) {
  flowField.getPressure().getScalar(i + 1,j) = rightBuffer[j]; //j changes, i+1 for ghost
}

void Stencils::PressureBufferReadStencil::applyBottomWall(
  FlowField& flowField, int i, int j
) {
  flowField.getPressure().getScalar(i,j - 1) = bottomBuffer[i]; //i changes
}

void Stencils::PressureBufferReadStencil::applyTopWall(
  FlowField& flowField, int i, int j
) {
  flowField.getPressure().getScalar(i,j + 1) = topBuffer[i]; //i changes
}

//3d cases
void Stencils::PressureBufferReadStencil::applyLeftWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getPressure().getScalar(i - 1,j,k) = leftBuffer[j*cellsZ + k]; //j,k changes
}

void Stencils::PressureBufferReadStencil::applyRightWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getPressure().getScalar(i + 1,j,k) = rightBuffer[j*cellsZ + k]; //j,k changes
}

void Stencils::PressureBufferReadStencil::applyBottomWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getPressure().getScalar(i - 1,j,k) = bottomBuffer[i*cellsZ + k]; //i,k changes
}

void Stencils::PressureBufferReadStencil::applyTopWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getPressure().getScalar(i + 1,j,k) = topBuffer[i*cellsZ + k]; //i,k changes
}

void Stencils::PressureBufferReadStencil::applyFrontWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getPressure().getScalar(i - 1,j,k) = frontBuffer[i*cellsY + j]; //i,j changes
}

void Stencils::PressureBufferReadStencil::applyBackWall(
  FlowField& flowField, int i, int j, int k
) {
  flowField.getPressure().getScalar(i + 1,j,k) = backBuffer[i*cellsY + j]; //i,j changes
}