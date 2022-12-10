#include "StdAfx.hpp"

#include "TurbulentSimulation.hpp"

#include "Solvers/PetscSolver.hpp"
#include "Solvers/SORSolver.hpp"

TurbulentSimulation::TurbulentSimulation(Parameters& parameters, FlowField& flowField):
  Simulation(parameters,flowField),
  fghTurbStencil_(parameters),
  fghTurbIterator_(flowField_, parameters, fghTurbStencil_),
  timeStepStencil_(parameters),
  timeStepIterator_(flowField_, parameters, timeStepStencil_),
  vtStencil_(parameters),
  vtIterator_(flowField_,parameters, vtStencil_)
{
}

void TurbulentSimulation::initializeFlowField() {
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

    // Initializing the "h" scalar field with the distances from the nearest walls
    Stencils::hStencil       hStencil(parameters_);
    FieldIterator<FlowField> hIterator(flowField_, parameters_, hStencil, 0, 1);
    hIterator.iterate();

    // Initialize mixing length scalar field
    Stencils::lmStencil      lmStencil(parameters_);
    FieldIterator<FlowField> lmIterator(flowField_, parameters_, lmStencil, 0, 1);
    lmIterator.iterate();

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

void TurbulentSimulation::solveTimestep() {
  // Determine and set max. timestep which is allowed in this simulation
  timeStepIterator_.iterate();
  setTimeStep();
  // Compute turbulent viscosity
  vtIterator_.iterate();
  // Compute FGH
  fghTurbIterator_.iterate();
  // Set global boundary values
  wallFGHIterator_.iterate();
  // Compute the right hand side (RHS)
  rhsIterator_.iterate();
  // Solve for pressure
  solver_->solve();
  // TODO WS2: communicate pressure values
  // Compute velocity
  velocityIterator_.iterate();
  obstacleIterator_.iterate();
  // TODO WS2: communicate velocity values
  // Iterate for velocities on the boundary
  wallVelocityIterator_.iterate();
}

void TurbulentSimulation::plotVTK(int timeStep, RealType simulationTime) {
  Stencils::VTKTurbStencil vtkStencil(parameters_);
  FieldIterator<FlowField> vtkIterator(flowField_, parameters_, vtkStencil, 1, 0);

  vtkIterator.iterate();
  vtkStencil.write(flowField_, timeStep, simulationTime);
}

void TurbulentSimulation::setTimeStep() {
  ASSERTION(parameters_.geometry.dim == 2 || parameters_.geometry.dim == 3);
  
  RealType localMin, globalMin;
  localMin = timeStepStencil_.dt;

  // Here, we select the type of operation before compiling. This allows to use the correct
  // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
  // machines.

  globalMin = MY_FLOAT_MAX;
  MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

  parameters_.timestep.dt = globalMin;
  parameters_.timestep.dt *= parameters_.timestep.tau;
}
