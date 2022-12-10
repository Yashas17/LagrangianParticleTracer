#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /**Stencil to compute the time step for the turbulent simulation.
   */
  class TimeStepStencil: public FieldStencil<FlowField> {
  public:
    TimeStepStencil(const Parameters& parameters);
    ~TimeStepStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
