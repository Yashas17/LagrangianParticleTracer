#pragma once

#include "Parameters.hpp"

namespace Stencils {

  /** Interface for operations on the global (or parallel)
   */
  template <class FlowFieldType>
  class GhostLayerStencil {
  protected:
    const Parameters& parameters_;

  public:
    GhostLayerStencil(const Parameters& parameters):
      parameters_(parameters) {}

    virtual ~GhostLayerStencil() = default;

    /** Represents an operation in the left wall of a 2D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyLeftWall2D(FlowFieldType& flowField) = 0;

    /** Represents an operation in the right wall of a 2D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyRightWall2D(FlowFieldType& flowField) = 0;

    /** Represents an operation in the bottom wall of a 2D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyBottomWall2D(FlowFieldType& flowField) = 0;

    /** Represents an operation in the top wall of a 2D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyTopWall2D(FlowFieldType& flowField) = 0;

    /** Represents an operation in the left wall of a 3D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyLeftWall2D(FlowFieldType& flowField) = 0;

    /** Represents an operation in the right wall of a 3D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyRightWall3D(FlowFieldType& flowField) = 0;

    /** Represents an operation in the bottom wall of a 3D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyBottomWall3D(FlowFieldType& flowField) = 0;

    /** Represents an operation in the top wall of a 3D domain.
     *
     * @param flowField State of the flow field
    virtual void applyTopWall(FlowFieldType& flowField) = 0;

    /** Represents an operation in the front wall of a 3D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyFrontWall3D(FlowFieldType& flowField) = 0;

    /** Represents an operation in the back wall of a 3D domain.
     *
     * @param flowField State of the flow field
     */
    virtual void applyBackWall3D(FlowFieldType& flowField) = 0;
  };

} // namespace Stencils
