#pragma once

#include "../Definitions.hpp"
#include "../FlowField.hpp"
#include "../Iterators.hpp"
#include "../Parameters.hpp"
#include "../Stencils/PressureBufferFillStencil.hpp"
#include "../Stencils/PressureBufferReadStencil.hpp"
#include "../Stencils/VelocityBufferFillStencil.hpp"
#include "../Stencils/VelocityBufferReadStencil.hpp"
#include "../Stencils/VtBufferFillStencil.hpp"
#include "../Stencils/VtBufferReadStencil.hpp"

namespace ParallelManagers {

  /** Class used to communicate pressure and velocity between MPI processes
   */
  class PetscParallelManager {
  private:
    Parameters& parameters_; //! Reference to the parameters

    FlowField& flowField_;

    Stencils::PressureBufferFillStencil pressureBufferFillStencil_;
    ParallelBoundaryIterator<FlowField> pressureBufferFillIterator_;

    Stencils::VelocityBufferFillStencil velocityBufferFillStencil_;
    ParallelBoundaryIterator<FlowField> velocityBufferFillIterator_;

    Stencils::VtBufferFillStencil       vtBufferFillStencil_;
    ParallelBoundaryIterator<FlowField> vtBufferFillIterator_;

  public:
    void communicatePressure();
    void communicateVelocities();
    void communicateVt();
    PetscParallelManager(Parameters& parameters, FlowField& flowField);
    ~PetscParallelManager() = default;
  };

} // namespace ParallelManagers
