#pragma once

#include "Definitions.hpp"
#include "FlowField.hpp"
#include "GlobalBoundaryFactory.hpp"
#include "Iterators.hpp"
#include "Simulation.hpp"

#include "Solvers/LinearSolver.hpp"
#include "Stencils/BFInputStencils.hpp"
#include "Stencils/BFStepInitStencil.hpp"
#include "Stencils/FGHTurbStencil.hpp"
#include "Stencils/hStencil.hpp"
#include "Stencils/InitTaylorGreenFlowFieldStencil.hpp"
#include "Stencils/lmStencil.hpp"
#include "Stencils/MaxUStencil.hpp"
#include "Stencils/MovingWallStencils.hpp"
#include "Stencils/NeumannBoundaryStencils.hpp"
#include "Stencils/ObstacleStencil.hpp"
#include "Stencils/PeriodicBoundaryStencils.hpp"
#include "Stencils/RHSStencil.hpp"
#include "Stencils/TimeStepStencil.hpp"
#include "Stencils/VelocityStencil.hpp"
#include "Stencils/VTKTurbStencil.hpp"
#include "Stencils/VtStencil.hpp"

class TurbulentSimulation: public Simulation {
protected:
  Stencils::FGHTurbStencil fghTurbStencil_;
  FieldIterator<FlowField> fghTurbIterator_;

  Stencils::TimeStepStencil timeStepStencil_;
  FieldIterator<FlowField>  timeStepIterator_;

  Stencils::VtStencil      vtStencil_;
  FieldIterator<FlowField> vtIterator_;

public:
  TurbulentSimulation(Parameters& parameters, FlowField& flowField);
  ~TurbulentSimulation() = default;

  /** Initialises the flow field according to the scenario */
  void initializeFlowField() override;

  void solveTimestep() override;

  void setTimeStep() override;

  /** Plots the flow field */
  void plotVTK(int timeStep, RealType simulationTime) override;
};