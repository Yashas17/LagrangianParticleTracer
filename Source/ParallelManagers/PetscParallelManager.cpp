#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"

#include <fstream>

#include "../Stencils/PressureBufferReadStencil.hpp"
#include "../Stencils/VelocityBufferReadStencil.hpp"

void check_mpi_error(int err_code) {
  if (err_code != MPI_SUCCESS) {
    char* c   = new char[MPI_MAX_ERROR_STRING];
    int   len = 0;
    MPI_Error_string(err_code, c, &len);
    std::cout << err_code << std::endl;
    throw std::runtime_error(c);
  }
}

void dump_flowField(FlowField& flowField_, Parameters& parameters_) {
  static int        call = 0;
  int               rank = parameters_.parallel.rank;
  std::ofstream     of("rank" + std::to_string(rank) + ".txt", std::ios_base::app);
  std::stringstream ss;
  ss.setf(std::ios::fixed);
  ss.precision(4);
  ss << "Timestep " << (call >> 1) << ":\n";
  if (call % 2 == 0) {
    ss << "Before communicating the ghost layers\n";
  } else {
    ss << "After communicating the ghost layers\n";
  }
  ss << parameters_.parallel.rank << " has neighbors: "
     << "\tleft:" << parameters_.parallel.leftNb << "\tright:" << parameters_.parallel.rightNb
     << "\ttop:" << parameters_.parallel.topNb << "\tbottom:" << parameters_.parallel.bottomNb << std::endl;
  ss << "Pressure of rank " << rank << ":\n";
  for (int j = flowField_.getPressure().getNy() - 1; j >= 0; j--) {
    for (int i = 0; i < flowField_.getPressure().getNx(); i++) {
      if (flowField_.getPressure().getScalar(i, j) >= 0.0) {
        ss << " " << flowField_.getPressure().getScalar(i, j) << ", ";
      } else {
        ss << flowField_.getPressure().getScalar(i, j) << ", ";
      }
    }
    ss << "endrow\n";
  }
  ss << "Velocity U of rank " << rank << ":\n";
  for (int j = flowField_.getVelocity().getNy() - 1; j >= 0; j--) {
    for (int i = 0; i < flowField_.getVelocity().getNx(); i++) {
      if (flowField_.getVelocity().getVector(i, j)[0] >= 0.0) {
        ss << " " << flowField_.getVelocity().getVector(i, j)[0] << ", ";
      } else {
        ss << flowField_.getVelocity().getVector(i, j)[0] << ", ";
      }
    }
    ss << "endrow\n";
  }
  ss << "Velocity V of rank " << rank << ":\n";
  for (int j = flowField_.getVelocity().getNy() - 1; j >= 0; j--) {
    for (int i = 0; i < flowField_.getVelocity().getNx(); i++) {
      if (flowField_.getVelocity().getVector(i, j)[1] >= 0.0) {
        ss << " " << flowField_.getVelocity().getVector(i, j)[1] << ", ";
      } else {
        ss << flowField_.getVelocity().getVector(i, j)[1] << ", ";
      }
    }
    ss << "endrow\n";
  }
  of << ss.str();
  call += 1;
}

ParallelManagers::PetscParallelManager::PetscParallelManager(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField),
  pressureBufferFillStencil_(parameters),
  pressureBufferFillIterator_(flowField, parameters, pressureBufferFillStencil_, 1, -1),
  // the read stencils will need the data to be initialized
  velocityBufferFillStencil_(parameters),
  velocityBufferFillIterator_(flowField, parameters, velocityBufferFillStencil_, 1, -1) {}

void ParallelManagers::PetscParallelManager::communicatePressure() {

  if (parameters_.geometry.dim == 2) { // 2d case
    dump_flowField(flowField_, parameters_);

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
      flowField_, parameters_, pressureBufferReadStencil_, 1, -1
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
      flowField_, parameters_, pressureBufferReadStencil_, 1, -1
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
    std::vector<RealType> leftBufferRecv(4 * cellsY);
    std::vector<RealType> rightBufferRecv(2 * cellsY);
    std::vector<RealType> bottomBufferRecv(4 * cellsX);
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
      4 * cellsY,
      MY_MPI_FLOAT,
      parameters_.parallel.rightNb,
      1,
      leftBufferRecv.data(),
      4 * cellsY,
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
      4 * cellsX,
      MY_MPI_FLOAT,
      parameters_.parallel.topNb,
      3,
      bottomBufferRecv.data(),
      4 * cellsX,
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
      flowField_, parameters_, velocityBufferReadStencil_, 1, -1
    );
    velocityBufferReadIterator_.iterate();

    dump_flowField(flowField_, parameters_);
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
      flowField_, parameters_, velocityBufferReadStencil_, 1, -1
    );
    velocityBufferReadIterator_.iterate();
  }
}