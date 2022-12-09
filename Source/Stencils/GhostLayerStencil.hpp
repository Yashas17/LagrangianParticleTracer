#pragma once

#include "Parameters.hpp"

namespace Stencils {

  /** Interface for operations on the global (or parallel) boundary
   */
  template <class FlowFieldType>
  class GhostLayerStencil {
  protected:
    const Parameters& parameters_;

  public:
    GhostLayerStencil(const Parameters& parameters):
      parameters_(parameters) {}

    virtual ~GhostLayerStencil() = default;

    virtual void applyLeftWall2D(FlowFieldType& flowField) = 0;

    virtual void applyRightWall2D(FlowFieldType& flowField) = 0;

    virtual void applyBottomWall2D(FlowFieldType& flowField) = 0;

    virtual void applyTopWall2D(FlowFieldType& flowField) = 0;

    virtual void applyLeftWall3D(FlowFieldType& flowField) = 0;

    virtual void applyRightWall3D(FlowFieldType& flowField) = 0;

    virtual void applyBottomWall3D(FlowFieldType& flowField) = 0;

    virtual void applyTopWall3D(FlowFieldType& flowField) = 0;

    virtual void applyFrontWall3D(FlowFieldType& flowField) = 0;

    virtual void applyBackWall3D(FlowFieldType& flowField) = 0;
  };

} // namespace Stencils
