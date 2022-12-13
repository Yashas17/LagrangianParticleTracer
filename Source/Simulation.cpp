#include "StdAfx.hpp"

#include "Simulation.hpp"

#include "Solvers/PetscSolver.hpp"
#include "Solvers/SORSolver.hpp"

Simulation::Simulation(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField),
  maxUStencil_(parameters),
  maxUFieldIterator_(flowField_, parameters, maxUStencil_),
  maxUBoundaryIterator_(flowField_, parameters, maxUStencil_),
  globalBoundaryFactory_(parameters),
  wallVelocityIterator_(globalBoundaryFactory_.getGlobalBoundaryVelocityIterator(flowField_)),
  wallFGHIterator_(globalBoundaryFactory_.getGlobalBoundaryFGHIterator(flowField_)),
  fghStencil_(parameters),
  fghIterator_(flowField_, parameters, fghStencil_),
  velocityStencil_(parameters),
  obstacleStencil_(parameters),
  velocityIterator_(flowField_, parameters, velocityStencil_),
  obstacleIterator_(flowField_, parameters, obstacleStencil_),
  rhsStencil_(parameters),
  rhsIterator_(flowField_, parameters, rhsStencil_),
  petscParallelManager_(parameters, flowField_),
#ifdef ENABLE_PETSC
  solver_(std::make_unique<Solvers::PetscSolver>(flowField_, parameters))
#else
  solver_(std::make_unique<Solvers::SORSolver>(flowField_, parameters))
#endif
{
}

void Simulation::initializeFlowField() {
  if (parameters_.simulation.scenario == "taylor-green") {
    // Currently, a particular initialisation is only required for the taylor-green vortex.
    Stencils::InitTaylorGreenFlowFieldStencil stencil(parameters_);
    FieldIterator<FlowField>                  iterator(flowField_, parameters_, stencil);
    iterator.iterate();
  } else if (parameters_.simulation.scenario == "channel") {
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(flowField_, parameters_, stencil, 0, 1);
    iterator.iterate();
    wallVelocityIterator_.iterate();
  } else if (parameters_.simulation.scenario == "pressure-channel") {
    // Set pressure boundaries here for left wall
    const RealType value = parameters_.walls.scalarLeft;
    ScalarField&   rhs   = flowField_.getRHS();

    if (parameters_.geometry.dim == 2) {
      const int sizey = flowField_.getNy();
      for (int i = 0; i < sizey + 3; i++) {
        rhs.getScalar(0, i) = value;
      }
    } else {
      const int sizey = flowField_.getNy();
      const int sizez = flowField_.getNz();
      for (int i = 0; i < sizey + 3; i++) {
        for (int j = 0; j < sizez + 3; j++) {
          rhs.getScalar(0, i, j) = value;
        }
      }
    }

    // Do same procedure for domain flagging as for regular channel
    Stencils::BFStepInitStencil stencil(parameters_);
    FieldIterator<FlowField>    iterator(flowField_, parameters_, stencil, 0, 1);
    iterator.iterate();
  }

  solver_->reInitMatrix();
}

void Simulation::solveTimestep() {
  // Determine and set max. timestep which is allowed in this simulation
  setTimeStep();
  // Compute FGH
  fghIterator_.iterate();
  // Set global boundary values
  wallFGHIterator_.iterate();
  // Compute the right hand side (RHS)
  rhsIterator_.iterate();
  // Solve for pressure
  solver_->solve();

  // TODO WS2: communicate pressure values
  petscParallelManager_.communicatePressure();
  //double communication to ensure the corner values are communicated
  petscParallelManager_.communicatePressure();

  // Compute velocity
  velocityIterator_.iterate();
  obstacleIterator_.iterate();

  // TODO WS2: communicate velocity values
  petscParallelManager_.communicateVelocities();
  //double communication to ensure the corner values are communicated
  petscParallelManager_.communicateVelocities();

  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterate();
}

void Simulation::plotVTK(int timeStep, RealType simulationTime) {
#ifndef DISABLE_OUTPUT
  Stencils::VTKStencil     vtkStencil(parameters_);
  FieldIterator<FlowField> vtkIterator(flowField_, parameters_, vtkStencil, 1, 0);

  vtkIterator.iterate();
  vtkStencil.write(flowField_, timeStep, simulationTime);
#endif
}

void Simulation::setTimeStep() {
  RealType localMin, globalMin;
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  RealType factor = 1.0 / (parameters_.meshsize->getDxMin() * parameters_.meshsize->getDxMin())
                    + 1.0 / (parameters_.meshsize->getDyMin() * parameters_.meshsize->getDyMin());

  // Determine maximum velocity
  maxUStencil_.reset();
  maxUFieldIterator_.iterate();
  maxUBoundaryIterator_.iterate();
  if (parameters_.geometry.dim == 3) {
    factor += 1.0 / (parameters_.meshsize->getDzMin() * parameters_.meshsize->getDzMin());
    /**
     *We need the following 'if' condition to ensure that we do not encounter division by zero. In case the 'if'
     *condition is not satisfied, the value of dt will be same as previous time step or if it is the first time step,
     *the user defined time-step will be used. If user has not defined the timestep, dt will be set to 1 by default.
     **/
    if (maxUStencil_.getMaxValues()[2] != 0)
      parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[2]);
  } else {
    if (maxUStencil_.getMaxValues()[0] != 0)
      parameters_.timestep.dt = 1.0 / (maxUStencil_.getMaxValues()[0]);
  }

  // localMin = std::min(parameters_.timestep.dt, std::min(std::min(parameters_.flow.Re/(2 * factor), 1.0 /
  // maxUStencil_.getMaxValues()[0]), 1.0 / maxUStencil_.getMaxValues()[1]));
  if (maxUStencil_.getMaxValues()[0] != 0 && maxUStencil_.getMaxValues()[1] != 0) {
    localMin = std::min(
      parameters_.flow.Re / (2 * factor),
      std::min(
        parameters_.timestep.dt, std::min(1 / (maxUStencil_.getMaxValues()[0]), 1 / (maxUStencil_.getMaxValues()[1]))
      )
    );
  } else {
    localMin = std::min(parameters_.flow.Re / (2 * factor), parameters_.timestep.dt);
  }

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}
