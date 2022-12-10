#pragma once

#include <algorithm>
#include <cmath>

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"

namespace Stencils {

  /** Initialises the backward facing step scenario, i.e. sets the flag field.
   */
  class lmStencil: public FieldStencil<FlowField> {
  private:
    RealType k_ = 0.041;

  public:
    lmStencil(const Parameters& parameters);
    ~lmStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils