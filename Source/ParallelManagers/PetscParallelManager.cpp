#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"

#include "../Stencils/PressureBufferReadStencil.hpp"
#include "../Stencils/VelocityBufferReadStencil.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField),
  pressureBufferFillStencil_(parameters),
  pressureBufferFillIterator_(flowField, parameters, pressureBufferFillStencil_, 2, -1),
  //the read stencils will need the data to be initialized
  velocityBufferFillStencil_(parameters),
  velocityBufferFillIterator_(flowField, parameters, velocityBufferFillStencil_, 2, -1)
{
}

void ParallelManagers::PetscParallelManager::communicatePressure() {

  if(parameters_.geometry.dim == 2){   //2d case

    int cellsX = flowField_.getCellsX();
    int cellsY = flowField_.getCellsY();

    //Fill the buffers
    pressureBufferFillIterator_.iterate();

    //prepare receive buffers, send buffers are already on the stencil
    std::vector<RealType> leftBufferRecv(cellsY);
    std::vector<RealType> rightBufferRecv(cellsY);
    std::vector<RealType> bottomBufferRecv(cellsX);
    std::vector<RealType> topBufferRecv(cellsX);

    //send from left, receive on right
    MPI_Sendrecv(&pressureBufferFillStencil_.leftBuffer[0],
                  cellsY,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.leftNb,
                  0,
                  &rightBufferRecv[0],
                  cellsY,
                  MPI_DOUBLE,
                  parameters_.parallel.rightNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from right, receive on left
    MPI_Sendrecv(&pressureBufferFillStencil_.rightBuffer[0],
                  cellsY,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.rightNb,
                  0,
                  &leftBufferRecv[0],
                  cellsY,
                  MPI_DOUBLE,
                  parameters_.parallel.leftNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from bottom, receive on top
    MPI_Sendrecv(&pressureBufferFillStencil_.bottomBuffer[0],
                  cellsX,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.bottomNb,
                  0,
                  &topBufferRecv[0],
                  cellsX,
                  MPI_DOUBLE,
                  parameters_.parallel.topNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from top, receive on bottom
    MPI_Sendrecv(&pressureBufferFillStencil_.topBuffer[0],
                  cellsX,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.topNb,
                  0,
                  &bottomBufferRecv[0],
                  cellsX,
                  MPI_DOUBLE,
                  parameters_.parallel.bottomNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //Read the buffers into the flowField
    Stencils::PressureBufferReadStencil pressureBufferReadStencil_(parameters_, 
                                                                  leftBufferRecv,
                                                                  rightBufferRecv,
                                                                  bottomBufferRecv,
                                                                  topBufferRecv);
    ParallelBoundaryIterator<FlowField> pressureBufferReadIterator_(flowField_,
                                                                    parameters_,
                                                                    pressureBufferReadStencil_,
                                                                    2,
                                                                    -1);
    pressureBufferReadIterator_.iterate();

  } else if (parameters_.geometry.dim == 3){ //3d case

    int cellsX = flowField_.getCellsX();
    int cellsY = flowField_.getCellsY();
    int cellsZ = flowField_.getCellsZ();

    //Fill the buffers
    pressureBufferFillIterator_.iterate();

    //prepare receive buffers, send buffers are already on the stencil
    std::vector<RealType> leftBufferRecv(cellsY*cellsZ);
    std::vector<RealType> rightBufferRecv(cellsY*cellsZ);
    std::vector<RealType> bottomBufferRecv(cellsX*cellsZ);
    std::vector<RealType> topBufferRecv(cellsX*cellsZ);
    std::vector<RealType> frontBufferRecv(cellsX*cellsY);
    std::vector<RealType> backBufferRecv(cellsX*cellsY);

    //send from left, receive on right
    MPI_Sendrecv(&pressureBufferFillStencil_.leftBuffer[0],
                  cellsY*cellsZ,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.leftNb,
                  0,
                  &rightBufferRecv[0],
                  cellsY*cellsZ,
                  MPI_DOUBLE,
                  parameters_.parallel.rightNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from right, receive on left
    MPI_Sendrecv(&pressureBufferFillStencil_.rightBuffer[0],
                  cellsY*cellsZ,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.rightNb,
                  0,
                  &leftBufferRecv[0],
                  cellsY*cellsZ,
                  MPI_DOUBLE,
                  parameters_.parallel.leftNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from bottom, receive on top
    MPI_Sendrecv(&pressureBufferFillStencil_.bottomBuffer[0],
                  cellsX*cellsZ,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.bottomNb,
                  0,
                  &topBufferRecv[0],
                  cellsX*cellsZ,
                  MPI_DOUBLE,
                  parameters_.parallel.topNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from top, receive on bottom
    MPI_Sendrecv(&pressureBufferFillStencil_.topBuffer[0],
                  cellsX*cellsZ,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.topNb,
                  0,
                  &bottomBufferRecv[0],
                  cellsX*cellsZ,
                  MPI_DOUBLE,
                  parameters_.parallel.bottomNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from front, receive on back
    MPI_Sendrecv(&pressureBufferFillStencil_.frontBuffer[0],
                  cellsX*cellsY,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.frontNb,
                  0,
                  &backBufferRecv[0],
                  cellsX*cellsY,
                  MPI_DOUBLE,
                  parameters_.parallel.backNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from back, receive on front
    MPI_Sendrecv(&pressureBufferFillStencil_.backBuffer[0],
                  cellsX*cellsY,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.backNb,
                  0,
                  &frontBufferRecv[0],
                  cellsX*cellsY,
                  MPI_DOUBLE,
                  parameters_.parallel.frontNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //Read the buffers into the flowField
    Stencils::PressureBufferReadStencil pressureBufferReadStencil_(parameters_, 
                                                                  leftBufferRecv,
                                                                  rightBufferRecv,
                                                                  bottomBufferRecv,
                                                                  topBufferRecv,
                                                                  frontBufferRecv,
                                                                  backBufferRecv);
    ParallelBoundaryIterator<FlowField> pressureBufferReadIterator_(flowField_,
                                                                    parameters_,
                                                                    pressureBufferReadStencil_,
                                                                    2,
                                                                    -1);
    pressureBufferReadIterator_.iterate();

  }

}

void ParallelManagers::PetscParallelManager::communicateVelocities() {
  if(parameters_.geometry.dim == 2){   //2d case

    int cellsX = flowField_.getCellsX();
    int cellsY = flowField_.getCellsY();

    //Fill the buffers
    velocityBufferFillIterator_.iterate();

    //prepare receive buffers, send buffers are already on the stencil
    std::vector<RealType> leftBufferRecv(3 * cellsY);
    std::vector<RealType> rightBufferRecv(2 * cellsY);
    std::vector<RealType> bottomBufferRecv(3 * cellsX);
    std::vector<RealType> topBufferRecv(2 * cellsX);

    //send from left, receive on right
    MPI_Sendrecv(&velocityBufferFillStencil_.leftBuffer[0],
                  2 * cellsY,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.leftNb,
                  0,
                  &rightBufferRecv[0],
                  2 * cellsY,
                  MPI_DOUBLE,
                  parameters_.parallel.rightNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from right, receive on left
    MPI_Sendrecv(&velocityBufferFillStencil_.rightBuffer[0],
                  3 * cellsY,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.rightNb,
                  0,
                  &leftBufferRecv[0],
                  3 * cellsY,
                  MPI_DOUBLE,
                  parameters_.parallel.leftNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from bottom, receive on top
    MPI_Sendrecv(&velocityBufferFillStencil_.bottomBuffer[0],
                  2 * cellsX,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.bottomNb,
                  0,
                  &topBufferRecv[0],
                  2 * cellsX,
                  MPI_DOUBLE,
                  parameters_.parallel.topNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from top, receive on bottom
    MPI_Sendrecv(&velocityBufferFillStencil_.topBuffer[0],
                  3 * cellsX,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.topNb,
                  0,
                  &bottomBufferRecv[0],
                  3 * cellsX,
                  MPI_DOUBLE,
                  parameters_.parallel.bottomNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //Read the buffers into the flowField
    Stencils::VelocityBufferReadStencil velocityBufferReadStencil_(parameters_, 
                                                                  leftBufferRecv,
                                                                  rightBufferRecv,
                                                                  bottomBufferRecv,
                                                                  topBufferRecv);
    ParallelBoundaryIterator<FlowField> velocityBufferReadIterator_(flowField_,
                                                                    parameters_,
                                                                    velocityBufferReadStencil_,
                                                                    2,
                                                                    -1);
    velocityBufferReadIterator_.iterate();

  } else if (parameters_.geometry.dim == 3){ //3d case

    int cellsX = flowField_.getCellsX();
    int cellsY = flowField_.getCellsY();
    int cellsZ = flowField_.getCellsZ();

    //Fill the buffers
    velocityBufferFillIterator_.iterate();

    //prepare receive buffers, send buffers are already on the stencil
    std::vector<RealType> leftBufferRecv(4 * cellsY*cellsZ);
    std::vector<RealType> rightBufferRecv(3 * cellsY*cellsZ);
    std::vector<RealType> bottomBufferRecv(4 * cellsX*cellsZ);
    std::vector<RealType> topBufferRecv(3 * cellsX*cellsZ);
    std::vector<RealType> frontBufferRecv(4 * cellsX*cellsY);
    std::vector<RealType> backBufferRecv(3 * cellsX*cellsY);

    //send from left, receive on right
    MPI_Sendrecv(&velocityBufferFillStencil_.leftBuffer[0],
                  3 * cellsY*cellsZ,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.leftNb,
                  0,
                  &rightBufferRecv[0],
                  3 * cellsY*cellsZ,
                  MPI_DOUBLE,
                  parameters_.parallel.rightNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from right, receive on left
    MPI_Sendrecv(&velocityBufferFillStencil_.rightBuffer[0],
                  4 * cellsY*cellsZ,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.rightNb,
                  0,
                  &leftBufferRecv[0],
                  4 * cellsY*cellsZ,
                  MPI_DOUBLE,
                  parameters_.parallel.leftNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from bottom, receive on top
    MPI_Sendrecv(&velocityBufferFillStencil_.bottomBuffer[0],
                  3 * cellsX*cellsZ,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.bottomNb,
                  0,
                  &topBufferRecv[0],
                  3 * cellsX*cellsZ,
                  MPI_DOUBLE,
                  parameters_.parallel.topNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from top, receive on bottom
    MPI_Sendrecv(&velocityBufferFillStencil_.topBuffer[0],
                  4 * cellsX*cellsZ,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.topNb,
                  0,
                  &bottomBufferRecv[0],
                  4 * cellsX*cellsZ,
                  MPI_DOUBLE,
                  parameters_.parallel.bottomNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from front, receive on back
    MPI_Sendrecv(&velocityBufferFillStencil_.frontBuffer[0],
                  3 * cellsX*cellsY,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.frontNb,
                  0,
                  &backBufferRecv[0],
                  3 * cellsX*cellsY,
                  MPI_DOUBLE,
                  parameters_.parallel.backNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //send from back, receive on front
    MPI_Sendrecv(&velocityBufferFillStencil_.backBuffer[0],
                  4 * cellsX*cellsY,
                  MPI_DOUBLE, //TODO: also make it work for float
                  parameters_.parallel.backNb,
                  0,
                  &frontBufferRecv[0],
                  4 * cellsX*cellsY,
                  MPI_DOUBLE,
                  parameters_.parallel.frontNb,
                  0,
                  PETSC_COMM_WORLD,
                  MPI_STATUS_IGNORE);

    //Read the buffers into the flowField
    Stencils::VelocityBufferReadStencil velocityBufferReadStencil_(parameters_, 
                                                                  std::move(leftBufferRecv),
                                                                  std::move(rightBufferRecv),
                                                                  std::move(bottomBufferRecv),
                                                                  std::move(topBufferRecv),
                                                                  std::move(frontBufferRecv),
                                                                  std::move(backBufferRecv));
    ParallelBoundaryIterator<FlowField> velocityBufferReadIterator_(flowField_,
                                                                    parameters_,
                                                                    velocityBufferReadStencil_,
                                                                    2,
                                                                    -1);
    velocityBufferReadIterator_.iterate();

  }
}