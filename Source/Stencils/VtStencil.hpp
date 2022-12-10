#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"

namespace Stencils {

  /** Initialises the backward facing step scenario, i.e. sets the flag field.
   */
  class VtStencil: public FieldStencil<FlowField> {
  private:
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];

  public:
    vtStencil(const Parameters& parameters);
    ~vtStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils