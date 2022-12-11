#include "StdAfx.hpp"

#include "PetscParallelManager.hpp"

void check_mpi_error(int err_code) {
  if (err_code != MPI_SUCCESS) {
    char* c   = new char[MPI_MAX_ERROR_STRING];
    int   len = 0;
    MPI_Error_string(err_code, c, &len);
    std::cout << err_code << std::endl;
    throw std::runtime_error(c);
  }
}

ParallelManagers::PetscParallelManagerNonBlocking::PetscParallelManagerNonBlocking(
  FlowField& flowField, Parameters& parameters
):
  parameters_(parameters),
  flowField_(flowField),
  pressureBufferFillStencil_(parameters),
  pressureBufferReadStencil_(parameters),
  pressureBufferFillIterator_(flowField, parameters, pressureBufferFillStencil_),
  pressureBufferReadIterator_(flowField, parameters, pressureBufferReadStencil_),
  velocityBufferFillStencil_(parameters),
  velocityBufferReadStencil_(parameters),
  velocityBufferFillIterator_(flowField, parameters, velocityBufferFillStencil_),
  velocityBufferReadIterator_(flowField, parameters, velocityBufferReadStencil_) {}

ParallelManagers::PetscParallelManagerNonBlocking::~PetscParallelManagerNonBlocking() {}

void ParallelManagers::PetscParallelManagerNonBlocking::communicatePressure() {
  if (parameters_.geometry.dim == 2) {
    pressureBufferFillIterator_.iterate();

    int cellsY = flowField_.getPressure().getNy();
    int cellsX = flowField_.getPressure().getNx();

    std::vector<RealType> right_recv(cellsY, std::numeric_limits<RealType>::quiet_NaN());
    std::vector<RealType> left_recv(cellsY, std::numeric_limits<RealType>::quiet_NaN());
    std::vector<RealType> top_recv(cellsX, std::numeric_limits<RealType>::quiet_NaN());
    std::vector<RealType> bottom_recv(cellsX, std::numeric_limits<RealType>::quiet_NaN());

    MPI_Request requests[8];
    MPI_Status  array_of_statuses[8];
    std::fill(requests, requests + 8, MPI_REQUEST_NULL);

    std::cout
      << parameters_.parallel.rank << " has neighbors: "
      << "\tleft:" << parameters_.parallel.leftNb << "\tright:" << parameters_.parallel.rightNb
      << "\ttop:" << parameters_.parallel.topNb << "\tbottom:" << parameters_.parallel.bottomNb << std::endl;

    if (parameters_.parallel.rightNb >= 0) {
      int rres = MPI_Irecv(
        right_recv.data(), cellsY, MY_MPI_FLOAT, parameters_.parallel.rightNb, 1, PETSC_COMM_WORLD, &requests[4]
      );
      check_mpi_error(rres);
    }

    if (parameters_.parallel.leftNb >= 0) {
      int rres = MPI_Irecv(
        left_recv.data(), cellsY, MY_MPI_FLOAT, parameters_.parallel.leftNb, 3, PETSC_COMM_WORLD, &requests[5]
      );
      check_mpi_error(rres);
    }

    if (parameters_.parallel.topNb >= 0) {
      int rres = MPI_Irecv(
        top_recv.data(), cellsX, MY_MPI_FLOAT, parameters_.parallel.topNb, 2, PETSC_COMM_WORLD, &requests[6]
      );
      check_mpi_error(rres);
    }

    if (parameters_.parallel.bottomNb >= 0) {
      int rres = MPI_Irecv(
        bottom_recv.data(), cellsX, MY_MPI_FLOAT, parameters_.parallel.bottomNb, 4, PETSC_COMM_WORLD, &requests[7]
      );
      check_mpi_error(rres);
    }

    // MPI Communication begin
    if (parameters_.parallel.leftNb >= 0) {
      int ires = MPI_Isend(
        pressureBufferFillStencil_.left_buffer.data(),
        cellsY,
        MY_MPI_FLOAT,
        parameters_.parallel.leftNb,
        1,
        PETSC_COMM_WORLD,
        &requests[0]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.rightNb >= 0) {
      int ires = MPI_Isend(
        pressureBufferFillStencil_.right_buffer.data(),
        cellsY,
        MY_MPI_FLOAT,
        parameters_.parallel.rightNb,
        3,
        PETSC_COMM_WORLD,
        &requests[1]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.bottomNb >= 0) {
      int ires = MPI_Isend(
        pressureBufferFillStencil_.bottom_buffer.data(),
        cellsX,
        MY_MPI_FLOAT,
        parameters_.parallel.bottomNb,
        2,
        PETSC_COMM_WORLD,
        &requests[2]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.topNb >= 0) {
      int ires = MPI_Isend(
        pressureBufferFillStencil_.top_buffer.data(),
        cellsX,
        MY_MPI_FLOAT,
        parameters_.parallel.topNb,
        4,
        PETSC_COMM_WORLD,
        &requests[3]
      );
      check_mpi_error(ires);
    }

    int res = MPI_Waitall(8, requests, array_of_statuses);
    // MPI Communication end
    if (res != MPI_SUCCESS) {
      char* c   = new char[MPI_MAX_ERROR_STRING];
      int   len = 0;
      MPI_Error_string(res, c, &len);
      std::cout << res << std::endl;
      throw std::runtime_error(c);
    }

    pressureBufferReadStencil_.left_buffer   = std::move(left_recv);
    pressureBufferReadStencil_.right_buffer  = std::move(right_recv);
    pressureBufferReadStencil_.top_buffer    = std::move(top_recv);
    pressureBufferReadStencil_.bottom_buffer = std::move(bottom_recv);

    pressureBufferReadIterator_.iterate();
  }
}

void ParallelManagers::PetscParallelManagerNonBlocking::communicateVelocity() {
  if (parameters_.geometry.dim == 2) {
    velocityBufferFillIterator_.iterate();
    assert(flowField_.getVelocity().getNy() == flowField_.getPressure().getNy());

    int cellsY = flowField_.getVelocity().getNy();
    int cellsX = flowField_.getVelocity().getNx();

    std::vector<RealType> right_recv(2 * cellsY, std::numeric_limits<RealType>::quiet_NaN());
    std::vector<RealType> left_recv(3 * cellsY, std::numeric_limits<RealType>::quiet_NaN());
    std::vector<RealType> top_recv(2 * cellsX, std::numeric_limits<RealType>::quiet_NaN());
    std::vector<RealType> bottom_recv(3 * cellsX, std::numeric_limits<RealType>::quiet_NaN());

    MPI_Request requests[8];
    MPI_Status  array_of_statuses[8];
    std::fill(requests, requests + 8, MPI_REQUEST_NULL);

    if (parameters_.parallel.rightNb >= 0) {
      int ires = MPI_Irecv(
        right_recv.data(), 2 * cellsY, MY_MPI_FLOAT, parameters_.parallel.rightNb, 1, PETSC_COMM_WORLD, &requests[4]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.leftNb >= 0) {
      int ires = MPI_Irecv(
        left_recv.data(), 3 * cellsY, MY_MPI_FLOAT, parameters_.parallel.leftNb, 3, PETSC_COMM_WORLD, &requests[5]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.topNb >= 0) {
      int ires = MPI_Irecv(
        top_recv.data(), 2 * cellsX, MY_MPI_FLOAT, parameters_.parallel.topNb, 2, PETSC_COMM_WORLD, &requests[6]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.bottomNb >= 0) {
      int ires = MPI_Irecv(
        bottom_recv.data(), 3 * cellsX, MY_MPI_FLOAT, parameters_.parallel.bottomNb, 4, PETSC_COMM_WORLD, &requests[7]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.leftNb >= 0) {
      int ires = MPI_Isend(
        velocityBufferFillStencil_.left_buffer.data(),
        2 * cellsY,
        MY_MPI_FLOAT,
        parameters_.parallel.leftNb,
        1,
        PETSC_COMM_WORLD,
        &requests[0]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.rightNb >= 0) {
      int ires = MPI_Isend(
        velocityBufferFillStencil_.right_buffer.data(),
        3 * cellsY,
        MY_MPI_FLOAT,
        parameters_.parallel.rightNb,
        3,
        PETSC_COMM_WORLD,
        &requests[1]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.bottomNb >= 0) {
      int ires = MPI_Isend(
        velocityBufferFillStencil_.bottom_buffer.data(),
        2 * cellsX,
        MY_MPI_FLOAT,
        parameters_.parallel.bottomNb,
        2,
        PETSC_COMM_WORLD,
        &requests[2]
      );
      check_mpi_error(ires);
    }

    if (parameters_.parallel.topNb >= 0) {
      int ires = MPI_Isend(
        velocityBufferFillStencil_.top_buffer.data(),
        3 * cellsX,
        MY_MPI_FLOAT,
        parameters_.parallel.topNb,
        4,
        PETSC_COMM_WORLD,
        &requests[3]
      );
      check_mpi_error(ires);
    }

    int res = MPI_Waitall(8, requests, array_of_statuses);
    if (res != MPI_SUCCESS) {
      char* c   = new char[MPI_MAX_ERROR_STRING];
      int   len = 0;
      MPI_Error_string(res, c, &len);
      throw std::runtime_error(c);
    }

    velocityBufferReadStencil_.left_buffer   = std::move(left_recv);
    velocityBufferReadStencil_.right_buffer  = std::move(right_recv);
    velocityBufferReadStencil_.top_buffer    = std::move(top_recv);
    velocityBufferReadStencil_.bottom_buffer = std::move(bottom_recv);

    velocityBufferReadIterator_.iterate();
  }
}
