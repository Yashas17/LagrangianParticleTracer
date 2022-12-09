#pragma once

#include "GhostLayerStencil.hpp"
#include "Parameters.hpp"

#include "../FlowField.hpp"

namespace Stencils {

  /** Interface for operations on the global (or parallel) boundary
   */
  class PressureBufferReadStencil: public GhostLayerStencil<FlowField> {
  protected:
    const Parameters& parameters_;

  public:
    std::vector<RealType> left_buffer;
    std::vector<RealType> right_buffer;
    std::vector<RealType> bottom_buffer;
    std::vector<RealType> top_buffer;
    std::vector<RealType> front_buffer;
    std::vector<RealType> back_buffer;

    PressureBufferReadStencil(const Parameters& parameters);

    ~PressureBufferReadStencil() override;

    void applyLeftWall2D(FlowField& flowField) override;

    void applyRightWall2D(FlowField& flowField) override;

    void applyBottomWall2D(FlowField& flowField) override;

    void applyTopWall2D(FlowField& flowField) override;

    void applyLeftWall3D(FlowField& flowField) override;

    void applyRightWall3D(FlowField& flowField) override;

    void applyBottomWall3D(FlowField& flowField) override;

    void applyTopWall3D(FlowField& flowField) override;

    void applyFrontWall3D(FlowField& flowField) override;

    void applyBackWall3D(FlowField& flowField) override;
  };

} // namespace Stencils