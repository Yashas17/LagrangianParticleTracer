#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"

#include "../Stencils/PressureBufferReadStencil.hpp"
#include "../Stencils/VelocityBufferReadStencil.hpp"

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField),
  pressureBufferFillStencil_(parameters),
  pressureBufferFillIterator_(flowField, parameters, pressureBufferFillStencil_, 1, 0),
  // the read stencils will need the data to be initialized
  velocityBufferFillStencil_(parameters),
  velocityBufferFillIterator_(flowField, parameters, velocityBufferFillStencil_, 1, 0) {}

void ParallelManagers::PetscParallelManager::communicatePressure() {

  if (parameters_.geometry.dim == 2) { // 2d case

    int cellsX = flowField_.getPressure().getNx();
    int cellsY = flowField_.getPressure().getNy();

    // std::cout << parameters_.parallel.localSize[0] << ", " << cellsX << std::endl;
    // std::cout << parameters_.parallel.localSize[1] << ", " << cellsY << std::endl;

    // Fill the buffers
    pressureBufferFillIterator_.iterate();

    // prepare receive buffers, send buffers are already on the stencil
    std::vector<RealType> leftBufferRecv(cellsY);
    std::vector<RealType> rightBufferRecv(cellsY);
    std::vector<RealType> bottomBufferRecv(cellsX);
    std::vector<RealType> topBufferRecv(cellsX);

    if (pressureBufferFillStencil_.leftBuffer.size() != cellsY) {
      throw std::runtime_error(
        std::to_string(pressureBufferFillStencil_.leftBuffer.size()) + "!=" + std::to_string(cellsY)
      );
    }
    // send from left, receive on right
    MPI_Sendrecv(
      pressureBufferFillStencil_.leftBuffer.data(),
      cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      0,
      rightBufferRecv.data(),
      cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      0,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from right, receive on left
    MPI_Sendrecv(
      pressureBufferFillStencil_.rightBuffer.data(),
      cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      1,
      leftBufferRecv.data(),
      cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      1,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from bottom, receive on top
    MPI_Sendrecv(
      pressureBufferFillStencil_.bottomBuffer.data(),
      cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      2,
      topBufferRecv.data(),
      cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      2,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from top, receive on bottom
    MPI_Sendrecv(
      pressureBufferFillStencil_.topBuffer.data(),
      cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3,
      bottomBufferRecv.data(),
      cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // Read the buffers into the flowField
    Stencils::PressureBufferReadStencil pressureBufferReadStencil_(
      parameters_,
      std::move(leftBufferRecv),
      std::move(rightBufferRecv),
      std::move(bottomBufferRecv),
      std::move(topBufferRecv)
    );
    ParallelBoundaryIterator<FlowField> pressureBufferReadIterator_(
      flowField_, parameters_, pressureBufferReadStencil_, 1, 0
    );
    pressureBufferReadIterator_.iterate();

  } else if (parameters_.geometry.dim == 3) { // 3d case

    int cellsX = flowField_.getPressure().getNx();
    int cellsY = flowField_.getPressure().getNy();
    int cellsZ = flowField_.getPressure().getNz();

    // Fill the buffers
    pressureBufferFillIterator_.iterate();

    // prepare receive buffers, send buffers are already on the stencil
    std::vector<RealType> leftBufferRecv(cellsY * cellsZ);
    std::vector<RealType> rightBufferRecv(cellsY * cellsZ);
    std::vector<RealType> bottomBufferRecv(cellsX * cellsZ);
    std::vector<RealType> topBufferRecv(cellsX * cellsZ);
    std::vector<RealType> frontBufferRecv(cellsX * cellsY);
    std::vector<RealType> backBufferRecv(cellsX * cellsY);

    // send from left, receive on right
    MPI_Sendrecv(
      pressureBufferFillStencil_.leftBuffer.data(),
      cellsY * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      0,
      rightBufferRecv.data(),
      cellsY * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      0,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from right, receive on left
    MPI_Sendrecv(
      pressureBufferFillStencil_.rightBuffer.data(),
      cellsY * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      1,
      leftBufferRecv.data(),
      cellsY * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      1,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from bottom, receive on top
    MPI_Sendrecv(
      pressureBufferFillStencil_.bottomBuffer.data(),
      cellsX * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      2,
      topBufferRecv.data(),
      cellsX * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      2,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from top, receive on bottom
    MPI_Sendrecv(
      pressureBufferFillStencil_.topBuffer.data(),
      cellsX * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3,
      bottomBufferRecv.data(),
      cellsX * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from front, receive on back
    MPI_Sendrecv(
      pressureBufferFillStencil_.frontBuffer.data(),
      cellsX * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      4,
      backBufferRecv.data(),
      cellsX * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      4,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from back, receive on front
    MPI_Sendrecv(
      pressureBufferFillStencil_.backBuffer.data(),
      cellsX * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      5,
      frontBufferRecv.data(),
      cellsX * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      5,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // Read the buffers into the flowField
    Stencils::PressureBufferReadStencil pressureBufferReadStencil_(
      parameters_,
      std::move(leftBufferRecv),
      std::move(rightBufferRecv),
      std::move(bottomBufferRecv),
      std::move(topBufferRecv),
      std::move(frontBufferRecv),
      std::move(backBufferRecv)
    );
    ParallelBoundaryIterator<FlowField> pressureBufferReadIterator_(
      flowField_, parameters_, pressureBufferReadStencil_, 1, 0
    );
    pressureBufferReadIterator_.iterate();
  }
}

void ParallelManagers::PetscParallelManager::communicateVelocities() {
  if (parameters_.geometry.dim == 2) { // 2d case

    int cellsX = flowField_.getVelocity().getNx();
    int cellsY = flowField_.getVelocity().getNy();

    // Fill the buffers
    velocityBufferFillIterator_.iterate();

    // prepare receive buffers, send buffers are already on the stencil
    std::vector<RealType> leftBufferRecv(3 * cellsY);
    std::vector<RealType> rightBufferRecv(2 * cellsY);
    std::vector<RealType> bottomBufferRecv(3 * cellsX);
    std::vector<RealType> topBufferRecv(2 * cellsX);

    // send from left, receive on right
    MPI_Sendrecv(
      velocityBufferFillStencil_.leftBuffer.data(),
      2 * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      0,
      rightBufferRecv.data(),
      2 * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      0,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from right, receive on left
    MPI_Sendrecv(
      velocityBufferFillStencil_.rightBuffer.data(),
      3 * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      1,
      leftBufferRecv.data(),
      3 * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      1,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from bottom, receive on top
    MPI_Sendrecv(
      velocityBufferFillStencil_.bottomBuffer.data(),
      2 * cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      2,
      topBufferRecv.data(),
      2 * cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      2,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from top, receive on bottom
    MPI_Sendrecv(
      velocityBufferFillStencil_.topBuffer.data(),
      3 * cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3,
      bottomBufferRecv.data(),
      3 * cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // Read the buffers into the flowField
    Stencils::VelocityBufferReadStencil velocityBufferReadStencil_(
      parameters_,
      std::move(leftBufferRecv),
      std::move(rightBufferRecv),
      std::move(bottomBufferRecv),
      std::move(topBufferRecv)
    );
    ParallelBoundaryIterator<FlowField> velocityBufferReadIterator_(
      flowField_, parameters_, velocityBufferReadStencil_, 1, 0
    );
    velocityBufferReadIterator_.iterate();

  } else if (parameters_.geometry.dim == 3) { // 3d case

    int cellsX = flowField_.getVelocity().getNx();
    int cellsY = flowField_.getVelocity().getNy();
    int cellsZ = flowField_.getVelocity().getNz();

    // Fill the buffers
    velocityBufferFillIterator_.iterate();

    // prepare receive buffers, send buffers are already on the stencil
    std::vector<RealType> leftBufferRecv(4 * cellsY * cellsZ);
    std::vector<RealType> rightBufferRecv(3 * cellsY * cellsZ);
    std::vector<RealType> bottomBufferRecv(4 * cellsX * cellsZ);
    std::vector<RealType> topBufferRecv(3 * cellsX * cellsZ);
    std::vector<RealType> frontBufferRecv(4 * cellsX * cellsY);
    std::vector<RealType> backBufferRecv(3 * cellsX * cellsY);

    // send from left, receive on right
    MPI_Sendrecv(
      velocityBufferFillStencil_.leftBuffer.data(),
      3 * cellsY * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      0,
      rightBufferRecv.data(),
      3 * cellsY * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      0,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from right, receive on left
    MPI_Sendrecv(
      velocityBufferFillStencil_.rightBuffer.data(),
      4 * cellsY * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      1,
      leftBufferRecv.data(),
      4 * cellsY * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.leftNb,
      1,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from bottom, receive on top
    MPI_Sendrecv(
      velocityBufferFillStencil_.bottomBuffer.data(),
      3 * cellsX * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      2,
      topBufferRecv.data(),
      3 * cellsX * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      2,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from top, receive on bottom
    MPI_Sendrecv(
      velocityBufferFillStencil_.topBuffer.data(),
      4 * cellsX * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3,
      bottomBufferRecv.data(),
      4 * cellsX * cellsZ,
      MY_MPI_FLOAT,
      parameters_.parallel.bottomNb,
      3,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from front, receive on back
    MPI_Sendrecv(
      velocityBufferFillStencil_.frontBuffer.data(),
      3 * cellsX * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      4,
      backBufferRecv.data(),
      3 * cellsX * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      4,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // send from back, receive on front
    MPI_Sendrecv(
      velocityBufferFillStencil_.backBuffer.data(),
      4 * cellsX * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.backNb,
      5,
      frontBufferRecv.data(),
      4 * cellsX * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.frontNb,
      5,
      PETSC_COMM_WORLD,
      MPI_STATUS_IGNORE
    );

    // Read the buffers into the flowField
    Stencils::VelocityBufferReadStencil velocityBufferReadStencil_(
      parameters_,
      std::move(leftBufferRecv),
      std::move(rightBufferRecv),
      std::move(bottomBufferRecv),
      std::move(topBufferRecv),
      std::move(frontBufferRecv),
      std::move(backBufferRecv)
    );
    ParallelBoundaryIterator<FlowField> velocityBufferReadIterator_(
      flowField_, parameters_, velocityBufferReadStencil_, 1, 0
    );
    velocityBufferReadIterator_.iterate();
  }
}