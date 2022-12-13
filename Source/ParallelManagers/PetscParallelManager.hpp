#pragma once

#include "../Definitions.hpp"
#include "../FlowField.hpp"
#include "../Iterators.hpp"
#include "../Parameters.hpp"
#include "../Stencils/PressureBufferFillStencil.hpp"
#include "../Stencils/PressureBufferReadStencil.hpp"
#include "../Stencils/VelocityBufferFillStencil.hpp"
#include "../Stencils/VelocityBufferReadStencil.hpp"

namespace ParallelManagers {

  /** Class used to communicate pressure and velocity between MPI processes
   */
  class PetscParallelManager {
  private:
    Parameters& parameters_; //! Reference to the parameters

    FlowField& flowField_;

    Stencils::PressureBufferFillStencil pressureBufferFillStencil_;
    // Stencils::PressureBufferReadStencil pressureBufferReadStencil_;
    ParallelBoundaryIterator<FlowField> pressureBufferFillIterator_;
    // ParallelBoundaryIterator<FlowField> pressureBufferReadIterator_;

    Stencils::VelocityBufferFillStencil velocityBufferFillStencil_;
    // Stencils::VelocityBufferReadStencil velocityBufferReadStencil_;
    ParallelBoundaryIterator<FlowField> velocityBufferFillIterator_;
    // ParallelBoundaryIterator<FlowField> velocityBufferReadIterator_;

  public:
    void communicatePressure();
    void communicateVelocities();
    PetscParallelManager(Parameters& parameters, FlowField& flowField);
    ~PetscParallelManager() = default;
  };

} // namespace ParallelManagers
