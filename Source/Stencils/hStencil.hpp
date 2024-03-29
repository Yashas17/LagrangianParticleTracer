#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"

namespace Stencils {

  /** Initialises the backward facing step scenario, i.e. sets the flag field.
   */
  class hStencil: public FieldStencil<FlowField> {
  private:
    const RealType xLimit_; //! Size of step in x-direction
    const RealType yLimit_; //! Same as for x

  public:
    hStencil(const Parameters& parameters);
    ~hStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils