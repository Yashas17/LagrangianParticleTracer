#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class FGHTurbStencil: public FieldStencil<FlowField> {
  private:
    // A local velocity variable that will be used to approximate derivatives. Size matches 3D
    // case, but can be used for 2D as well.
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];
    RealType localVTotal_[27];

  public:
    FGHTurbStencil(const Parameters& parameters);
    ~FGHTurbStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };
} // namespace Stencils