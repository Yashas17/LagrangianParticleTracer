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
  class PetscParallelManagerNonBlocking {
  private:
    Parameters& parameters_; //! Reference to the parameters
    FlowField&  flowField_;

    Stencils::PressureBufferFillStencil pressureBufferFillStencil_;
    Stencils::PressureBufferReadStencil pressureBufferReadStencil_;
    GhostLayerIterator<FlowField>       pressureBufferFillIterator_;
    GhostLayerIterator<FlowField>       pressureBufferReadIterator_;

    Stencils::VelocityBufferFillStencil velocityBufferFillStencil_;
    Stencils::VelocityBufferReadStencil velocityBufferReadStencil_;
    GhostLayerIterator<FlowField>       velocityBufferFillIterator_;
    GhostLayerIterator<FlowField>       velocityBufferReadIterator_;

  public:
    void communicatePressure();
    void communicateVelocity();

    PetscParallelManagerNonBlocking(FlowField& flowField, Parameters& parameters);
    ~PetscParallelManagerNonBlocking();
  };

} // namespace ParallelManagers