#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /**Stencil to compute RHS.
   */
  class RHSStencil: public FieldStencil<FlowField> {
  public:
    RHSStencil(const Parameters& parameters);
    ~RHSStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils
